# Prior compiler -----------------------------------------------------------

# This file translates symbolic SD/scale-parameterized prior rules into the
# resolved values stored in fit$prior and the metadata shown by prior_summary().
# JAGS precision conversion remains a separate final code-generation step.

# Serialize one numeric prior value for fit$prior and generated code.
# Uses the shortest exact representation, independent of options(digits).
format_prior_number = function(x) {
  if (length(x) != 1 || !is.numeric(x))
    stop_github("Expected one numeric prior value.")
  if (is.na(x))
    return(if (is.nan(x)) "NaN" else "NA")
  if (is.infinite(x))
    return(ifelse(x > 0, "Inf", "-Inf"))
  x = as.double(x)
  out = sprintf("%.15g", x)
  if (as.numeric(out) != x)
    out = sprintf("%.17g", x)
  out
}


# Deparse a prior expression without rounding numbers (deparse() uses 15 digits).
deparse_prior = function(expr) {
  out = deparse1(expr, backtick = TRUE)
  if (nzchar(out) && !identical(str2lang(out), expr))
    out = deparse1(expr, backtick = TRUE, control = c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17"))
  out
}


# Round numbers in prior strings for display in prior_summary(). Values used in the model remain exact.
format_prior_display = function(x, digits = 5) {
  stringr::str_replace_all(x, "(?<![\\w.])\\d+\\.?\\d*(?:[eE][+-]?\\d+)?", function(num) {
    vapply(as.numeric(num), function(value) format_prior_number(signif(value, digits)), character(1))
  })
}


# Quote non-syntactic data names for readable symbolic prior rules.
display_data_name = function(x) {
  if (make.names(x) == x && !grepl("^[.][0-9]", x)) x else paste0("`", x, "`")
}


# Collect data and model summaries used to resolve symbolic prior templates.
prior_context = function(data, segments, design_specs = list()) {
  par_x = segments$x[1]
  y_col = segments$y[1]
  x = data[[par_x]]

  # mu offsets are available as `offset` or `offset_<segment>` in prior expressions.
  # `offsets` names the symbol for each segment, e.g., c(`1` = "offset_1", `2` = "offset_2").
  mu_offsets = Filter(function(design_spec) isTRUE(design_spec$has_offset) && identical(design_spec$dpar, "mu") && any(design_spec$offset_data != 0), design_specs)
  offsets = character()
  for (design_spec in mu_offsets) {
    data[[paste0("offset_", design_spec$segment)]] = design_spec$offset_data
    data$offset = design_spec$offset_data
    offsets[[as.character(design_spec$segment)]] = if (length(mu_offsets) == 1) "offset" else paste0("offset_", design_spec$segment)
  }

  list(
    data = data,
    par_x = par_x,
    y_col = y_col,
    x_display = display_data_name(par_x),
    y_display = display_data_name(y_col),
    x_min = min(x, na.rm = TRUE),
    y_sd = suppressWarnings(stats::sd(as.numeric(data[[y_col]]), na.rm = TRUE)),
    offsets = offsets,
    n_cp = nrow(segments) - 1,
    n_segments = nrow(segments)
  )
}


# Substitute the actual predictor and response names into a prior template.
instantiate_prior_template = function(x, context) {
  # Whole names only, e.g., not the ".x" in "stateb.x_1"
  x = stringr::str_replace_all(x, "(?<![[:alnum:]._])\\.x(?![[:alnum:]._])", context$x_display)
  stringr::str_replace_all(x, "(?<![[:alnum:]._])\\.y(?![[:alnum:]._])", context$y_display)
}


# Parse a simple prior call into its function name and argument strings.
parse_prior_call = function(x) {
  if (is.null(x) || length(x) != 1 || is.na(x) || !is.character(x))
    return(NULL)
  expr = tryCatch(str2lang(trimws(x)), error = function(e) NULL)
  if (is.null(expr) || !is.call(expr) || !is.symbol(expr[[1]]))
    return(NULL)
  name = as.character(expr[[1]])
  if (!grepl("^[A-Za-z.][A-Za-z0-9._]*$", name))
    return(NULL)

  # Reject named arguments
  call_args = as.list(expr)[-1]
  arg_names = names(call_args)
  if (!is.null(arg_names) && any(nzchar(arg_names))) {
    target = if (name == "T") "truncation" else "distributions"
    stop("Named arguments are not supported in prior ", target, ": '", trimws(x), "'.")
  }
  args = vapply(call_args, deparse_prior, character(1), USE.NAMES = FALSE)
  list(name = name, args = args)
}


# Parse a prior into its distribution `dist`, argument strings `args`, and truncation bounds
# `trunc` ("" where unbounded), e.g., "dnorm(0, 1) T(0, )" into list(dist = "dnorm",
# args = c("0", "1"), trunc = c("0", "")). Constants, aliases, and expressions have
# dist = NA and the expression as their only argument. Priors are parsed once and written
# as strings with format_prior() (or prior_to_jags() for JAGS code).
parse_prior = function(x) {
  if (is.numeric(x))
    return(list(dist = NA_character_, args = format_prior_number(x), trunc = c("", "")))

  trunc = c("", "")
  match = regexpr("(?<=\\))\\s*T\\s*\\(", x, perl = TRUE)
  if (match[1] != -1) {
    call = parse_prior_call(substr(x, match[1], nchar(x)))
    if (is.null(call) || call$name != "T" || length(call$args) != 2)
      stop("Could not parse the truncation in prior '", x, "'. Use T(lower, upper), e.g., T(0, ).", call. = FALSE)
    trunc = call$args
    x = substr(x, 1, match[1] - 1L)
  }

  call = parse_prior_call(x)
  if (!is.null(call) && grepl("^d[A-Za-z]", call$name))
    return(list(dist = call$name, args = call$args, trunc = trunc))
  if (any(nzchar(trunc)))
    stop("Only distributions can be truncated. Got '", x, "' with T().", call. = FALSE)
  expr = tryCatch(str2lang(x), error = function(e) stop("Could not parse prior expression '", x, "'.", call. = FALSE))
  list(dist = NA_character_, args = deparse_prior(expr), trunc = c("", ""))
}


# Write a parsed prior (see parse_prior()) as a string, e.g., "dnorm(0, 1) T(0, )".
format_prior = function(parsed) {
  if (is.na(parsed$dist))
    return(parsed$args)
  out = paste0(parsed$dist, "(", paste(parsed$args, collapse = ", "), ")")
  if (any(nzchar(parsed$trunc)))
    out = paste0(out, " T(", paste(parsed$trunc, collapse = ", "), ")")
  out
}


# The lower and upper bounds of a parsed prior: the arguments of dunif() or T() ("" where unbounded).
get_prior_bounds = function(parsed) {
  if (!identical(parsed$dist, "dunif"))
    return(parsed$trunc)
  if (any(nzchar(parsed$trunc)))
    stop("JAGS cannot truncate dunif(). Set its lower and upper bounds as its arguments instead, e.g., \"dunif(20, 50)\".", call. = FALSE)
  parsed$args
}


set_prior_bounds = function(parsed, bounds) {
  if (identical(parsed$dist, "dunif")) parsed$args = bounds else parsed$trunc = bounds
  parsed
}


# Apply `shift()` to the location and the bounds of a parsed prior, e.g., to write the
# prior on x + c from the prior on x. Both arguments of dunif() are bounds.
shift_prior = function(parsed, shift) {
  location = if (parsed$dist == "dunif") 1:2 else 1
  parsed$args[location] = vapply(parsed$args[location], shift, character(1), USE.NAMES = FALSE)
  parsed$trunc = vapply(parsed$trunc, function(bound) if (nzchar(bound)) shift(bound) else "", character(1), USE.NAMES = FALSE)
  parsed
}


# Check whether an expression contains only approved data and math operations.
allowed_data_expression = function(expr, data_names) {
  if (is.numeric(expr))
    return(TRUE)
  if (is.name(expr))
    return(as.character(expr) %in% c(data_names, "Inf", "pi"))
  if (!is.call(expr))
    return(FALSE)
  func_name = as.character(expr[[1]])
  if (func_name %notin% c(
    "+", "-", "*", "/", "^", "(", "log", "exp", "sqrt", "abs",
    "min", "max", "pmin", "pmax", "round", "ifelse", "<="
  )) return(FALSE)
  all(vapply(as.list(expr)[-1], allowed_data_expression, logical(1), data_names = data_names))
}


# Resolve prior-expression symbols, summaries, and arithmetic where possible.
resolve_prior_ast = function(expr, context) {
  if (is.numeric(expr))
    return(as.numeric(expr))
  if (is.name(expr)) {
    legacy_value = legacy_prior_value(as.character(expr), context)
    return(if (is.null(legacy_value)) expr else unname(legacy_value))
  }
  if (!is.call(expr))
    return(expr)

  func_name = as.character(expr[[1]])
  args = as.list(expr)[-1]
  data_names = names(context$data)

  if (func_name %in% c("n_cp", "n_segments")) {
    if (length(args) != 0)
      stop("`", func_name, "()` does not take arguments in prior expressions.")
    return(if (func_name == "n_cp") context$n_cp else context$n_segments)
  }

  # In change point priors, min(x) and max(x) are the JAGS nodes cp_0 and cp_K
  if (isTRUE(context$cp_limits) && func_name %in% c("min", "max") && identical(args, list(as.name(context$par_x))))
    return(as.name(paste0("cp_", if (func_name == "min") 0 else context$n_segments)))

  summaries = c("min", "max", "mean", "median", "sd", "mad")
  if (func_name %in% summaries && length(args) == 1 && allowed_data_expression(args[[1]], data_names)) {
    data_env = list2env(as.list(context$data), parent = baseenv())
    values = eval(args[[1]], envir = data_env)
    return(as.numeric(switch(
      func_name,
      min = min(values, na.rm = TRUE),
      max = max(values, na.rm = TRUE),
      mean = mean(values, na.rm = TRUE),
      median = stats::median(values, na.rm = TRUE),
      sd = stats::sd(values, na.rm = TRUE),
      mad = stats::mad(values, na.rm = TRUE)
    )))
  }

  resolved_args = lapply(args, resolve_prior_ast, context = context)
  resolved = as.call(c(list(expr[[1]]), resolved_args))
  if (allowed_data_expression(resolved, character()))
    return(as.numeric(eval(resolved, envir = baseenv())))
  resolved
}


# Parse and resolve one symbolic prior expression to displayable code.
resolve_prior_expression = function(x, context) {
  expr = tryCatch(parse(text = x, keep.source = FALSE)[[1]], error = function(e) NULL)
  if (is.null(expr))
    stop("Could not parse prior expression '", x, "'.")
  resolved = resolve_prior_ast(expr, context)
  if (is.numeric(resolved))
    return(format_prior_number(resolved))
  deparse_prior(resolved)
}


# Resolve every argument and bound in one parsed prior.
resolve_prior = function(parsed, context) {
  parsed$args = vapply(parsed$args, resolve_prior_expression, character(1), context = context, USE.NAMES = FALSE)
  parsed$trunc = vapply(parsed$trunc, function(bound) {
    if (nzchar(bound)) resolve_prior_expression(bound, context) else ""
  }, character(1), USE.NAMES = FALSE)
  parsed
}


# Quantile of a prior in mcp's parameterization (scale, not precision) with numeric
# arguments `args` and truncation `bounds`. Returns NULL for unsupported distributions.
prior_quantile = function(dist, args, bounds, prob) {
  q_p = switch(dist,
    dnorm = list(function(p) stats::qnorm(p, args[1], args[2]), function(x) stats::pnorm(x, args[1], args[2])),
    dt = list(function(p) args[1] + args[2] * stats::qt(p, args[3]), function(x) stats::pt((x - args[1]) / args[2], args[3])),
    dbeta = list(function(p) stats::qbeta(p, args[1], args[2]), function(x) stats::pbeta(x, args[1], args[2])),
    dgamma = list(function(p) stats::qgamma(p, args[1], args[2]), function(x) stats::pgamma(x, args[1], args[2])),
    dunif = list(function(p) stats::qunif(p, args[1], args[2]), function(x) stats::punif(x, args[1], args[2]))
  )
  if (is.null(q_p))
    return(NULL)
  prob_bounds = q_p[[2]](bounds)
  q_p[[1]](prob_bounds[1] + prob * diff(prob_bounds))
}


# The lower and upper bounds of a parsed prior for prior_summary(): the arguments of dunif()
# or T(), or the support of the distribution.
describe_prior_bounds = function(parsed, context = list()) {
  if (is.na(parsed$dist))
    return("none")
  bounds = get_prior_bounds(parsed)
  bounds[!nzchar(bounds)] = c("-Inf", "Inf")[!nzchar(bounds)]
  if (!identical(bounds, c("-Inf", "Inf")))
    return(paste0("[", bounds[1], ", ", bounds[2], "]"))
  if (parsed$dist == "dirichlet")
    return(paste0("[", paste(vapply(range(context$data[[context$par_x]], na.rm = TRUE), format_prior_number, character(1)), collapse = ", "), "]"))
  if (parsed$dist == "dbeta") return("[0, 1]")
  if (parsed$dist %in% c("dgamma", "dinvgamma", "dexp", "dlnorm", "dweib")) return("[0, Inf]")
  "none"
}


# Classify a parsed prior as a distribution, alias, expression, or constant.
prior_kind = function(parsed, all_names) {
  if (!is.na(parsed$dist))
    return("distribution")
  if (!is.na(suppressWarnings(as.numeric(parsed$args))))
    return("constant")
  if (parsed$args %in% all_names)
    return("alias")
  "expression"
}


# Compile one prior into a row of the prior table:
# - value, parsed: resolved for JAGS, keeping cp_0 and cp_K for the ends of the observed range in change point priors.
# - prior, bounds: resolved for prior_summary(), with numbers for cp_0 and cp_K.
# - rule: the specified prior, e.g., with `sd(y)`.
compile_prior_record = function(parameter, parsed, all_names, context, source,
                                 description = NULL) {
  kind = prior_kind(parsed, all_names)
  display_original = legacy_translate_prior(format_prior(parsed), context)
  context$cp_limits = kind == "distribution" && parameter %in% context$cp_pars
  resolved = resolve_prior(parsed, context)
  value = format_prior(resolved)
  # Show the exact ends of the observed range, cp_0 and cp_K, as numbers
  display_context = context
  display_context$cp_limits = FALSE
  display = resolve_prior(parsed, display_context)
  bounds = describe_prior_bounds(display, context)
  if (context$cp_limits && parameter %in% context$cp_names)
    bounds = paste0(bounds, ", ordered")

  if (kind == "expression") {
    expr = tryCatch(parse(text = value, keep.source = FALSE)[[1]], error = function(e) NULL)
    unknown = if (is.null(expr)) character() else setdiff(all.vars(expr), all_names)
    if (length(unknown) > 0) {
      stop(
        "Unknown model parameter(s) in prior expression for '", parameter, "': ",
        and_collapse(unknown), "."
      )
    }
  }

  if (is.null(description) || is.na(description)) {
    description = switch(
      kind,
      distribution = if (source == "default") "Family-defined default prior" else "User-specified prior",
      alias = paste0("Same value as ", trimws(display_original)),
      expression = "Deterministic expression",
      constant = paste0("Fixed at ", format_prior_display(value))
    )
  }

  tibble::tibble(
    parameter = parameter,
    prior = format_prior_display(format_prior(display)),
    bounds = format_prior_display(bounds),
    rule = format_prior_display(display_original),
    description = description,
    source = source,
    kind = kind,
    value = value,
    parsed = list(resolved)
  )
}


# Compile all symbolic prior specifications into the final prior table.
compile_prior_specs = function(prior_specs, all_names, context) {
  required = c("parameter", "parsed", "description", "source")
  assert_data_cols(prior_specs, required)
  records = lapply(seq_len(nrow(prior_specs)), function(i) {
    compile_prior_record(
      prior_specs$parameter[i], prior_specs$parsed[[i]], all_names, context, prior_specs$source[i],
      prior_specs$description[i]
    )
  })
  dplyr::bind_rows(records)
}
