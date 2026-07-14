# Prior compiler -----------------------------------------------------------

# This file translates symbolic SD/scale-parameterized prior rules into the
# resolved values stored in fit$prior and the metadata shown by prior_summary().
# JAGS precision conversion remains a separate final code-generation step.

format_prior_number = function(x) {
  if (length(x) != 1 || !is.numeric(x))
    stop_github("Expected one numeric prior value.")
  if (is.na(x))
    return(if (is.nan(x)) "NaN" else "NA")
  if (is.infinite(x))
    return(ifelse(x > 0, "Inf", "-Inf"))
  format(signif(x, 15), scientific = FALSE, trim = TRUE)
}


display_data_name = function(x) {
  if (make.names(x) == x && !grepl("^[.][0-9]", x)) x else paste0("`", x, "`")
}


prior_context = function(data, ST) {
  x_name = ST$x[1]
  y_name = ST$y[1]
  x = data[[x_name]]
  list(
    data = data,
    x_name = x_name,
    y_name = y_name,
    x_display = display_data_name(x_name),
    y_display = display_data_name(y_name),
    x_min = min(x, na.rm = TRUE),
    x_max = max(x, na.rm = TRUE),
    x_span = diff(range(x, na.rm = TRUE)),
    n_cp = nrow(ST) - 1,
    n_segments = nrow(ST),
    segment_width = diff(range(x, na.rm = TRUE)) / nrow(ST)
  )
}


instantiate_prior_template = function(x, context) {
  x = stringr::str_replace_all(x, stringr::fixed(".x"), context$x_display)
  stringr::str_replace_all(x, stringr::fixed(".y"), context$y_display)
}


split_top_level = function(x, separator = ",") {
  if (!nzchar(trimws(x)))
    return(character())
  chars = strsplit(x, "", fixed = TRUE)[[1]]
  depth = 0L
  starts = 1L
  pieces = character()
  for (i in seq_along(chars)) {
    if (chars[i] == "(") depth = depth + 1L
    if (chars[i] == ")") depth = depth - 1L
    if (chars[i] == separator && depth == 0L) {
      pieces = c(pieces, substr(x, starts, i - 1L))
      starts = i + 1L
    }
  }
  c(pieces, substr(x, starts, nchar(x))) %>% stringr::str_trim()
}


split_prior_truncation = function(x) {
  match = regexpr("(?<=\\))\\s*T\\s*\\(", x, perl = TRUE)
  if (match[1] == -1)
    return(list(distribution = trimws(x), truncation = NULL))
  list(
    distribution = trimws(substr(x, 1, match[1] - 1L)),
    truncation = trimws(substr(x, match[1], nchar(x)))
  )
}


parse_prior_call = function(x) {
  x = trimws(x)
  if (!grepl("^[A-Za-z.][A-Za-z0-9._]*\\s*\\(", x))
    return(NULL)
  open = regexpr("\\(", x)[1]
  if (open < 2 || substr(x, nchar(x), nchar(x)) != ")")
    return(NULL)
  list(
    name = trimws(substr(x, 1, open - 1L)),
    args = split_top_level(substr(x, open + 1L, nchar(x) - 1L))
  )
}


allowed_data_expression = function(expr, data_names) {
  if (is.numeric(expr))
    return(TRUE)
  if (is.name(expr))
    return(as.character(expr) %in% c(data_names, "Inf", "pi"))
  if (!is.call(expr))
    return(FALSE)
  fun = as.character(expr[[1]])
  if (fun %notin% c(
    "+", "-", "*", "/", "^", "(", "log", "exp", "sqrt", "abs",
    "min", "max", "pmin", "pmax", "round"
  )) return(FALSE)
  all(vapply(as.list(expr)[-1], allowed_data_expression, logical(1), data_names = data_names))
}


resolve_prior_ast = function(expr, context) {
  if (is.numeric(expr))
    return(as.numeric(expr))
  if (is.name(expr)) {
    legacy_value = resolve_legacy_prior_name(as.character(expr), context)
    return(if (is.null(legacy_value)) expr else unname(legacy_value))
  }
  if (!is.call(expr))
    return(expr)

  fun = as.character(expr[[1]])
  args = as.list(expr)[-1]
  data_names = names(context$data)

  if (fun %in% c("n_cp", "n_segments")) {
    if (length(args) != 0)
      stop("`", fun, "()` does not take arguments in prior expressions.")
    return(if (fun == "n_cp") context$n_cp else context$n_segments)
  }

  summaries = c("min", "max", "mean", "median", "sd", "mad", "segment_width")
  if (fun %in% summaries && length(args) == 1 && allowed_data_expression(args[[1]], data_names)) {
    data_env = list2env(as.list(context$data), parent = baseenv())
    values = eval(args[[1]], envir = data_env)
    return(as.numeric(switch(
      fun,
      min = min(values, na.rm = TRUE),
      max = max(values, na.rm = TRUE),
      mean = mean(values, na.rm = TRUE),
      median = stats::median(values, na.rm = TRUE),
      sd = stats::sd(values, na.rm = TRUE),
      mad = stats::mad(values, na.rm = TRUE),
      segment_width = diff(range(values, na.rm = TRUE)) / context$n_segments
    )))
  }

  resolved_args = lapply(args, resolve_prior_ast, context = context)
  resolved = as.call(c(list(expr[[1]]), resolved_args))
  if (allowed_data_expression(resolved, character()))
    return(as.numeric(eval(resolved, envir = baseenv())))
  resolved
}


resolve_prior_expression = function(x, context) {
  expr = tryCatch(parse(text = x, keep.source = FALSE)[[1]], error = function(e) NULL)
  if (is.null(expr))
    stop("Could not parse prior expression '", x, "'.")
  resolved = resolve_prior_ast(expr, context)
  if (is.numeric(resolved))
    return(format_prior_number(resolved))
  paste(deparse(resolved, width.cutoff = 500L), collapse = "")
}


resolve_prior_value = function(x, context) {
  if (is.numeric(x))
    return(format_prior_number(x))
  parts = split_prior_truncation(x)
  call = parse_prior_call(parts$distribution)
  if (!is.null(call)) {
    args = vapply(call$args, resolve_prior_expression, character(1), context = context)
    value = paste0(call$name, "(", paste(args, collapse = ", "), ")")
  } else {
    value = resolve_prior_expression(parts$distribution, context)
  }
  if (!is.null(parts$truncation)) {
    trunc = parse_prior_call(trimws(parts$truncation))
    if (is.null(trunc) || trunc$name != "T")
      stop("Could not parse prior truncation '", parts$truncation, "'.")
    trunc_args = vapply(trunc$args, function(arg) {
      if (!nzchar(arg)) "" else resolve_prior_expression(arg, context)
    }, character(1))
    value = paste0(value, " T(", paste(trunc_args, collapse = ", "), ")")
  }
  value
}


human_prior_call = function(x) {
  parts = split_prior_truncation(x)
  call = parse_prior_call(parts$distribution)
  if (is.null(call))
    return(parts$distribution)
  a = call$args
  switch(
    call$name,
    dt = if (length(a) == 3) paste0("student_t(df = ", a[3], ", location = ", a[1], ", scale = ", a[2], ")") else parts$distribution,
    dnorm = if (length(a) == 2) paste0("normal(mean = ", a[1], ", sd = ", a[2], ")") else parts$distribution,
    dunif = if (length(a) == 2) paste0("uniform(min = ", a[1], ", max = ", a[2], ")") else parts$distribution,
    dbeta = if (length(a) == 2) paste0("beta(shape1 = ", a[1], ", shape2 = ", a[2], ")") else parts$distribution,
    dgamma = if (length(a) == 2) paste0("gamma(shape = ", a[1], ", rate = ", a[2], ")") else parts$distribution,
    dloginvgamma = if (length(a) == 2) paste0("log_inverse_gamma(shape = ", a[1], ", scale = ", a[2], ")") else parts$distribution,
    dirichlet = if (length(a) == 1) paste0("dirichlet(alpha = ", a[1], ")") else parts$distribution,
    parts$distribution
  )
}


prior_bounds = function(x) {
  parts = split_prior_truncation(x)
  if (!is.null(parts$truncation)) {
    trunc = parse_prior_call(trimws(parts$truncation))
    if (!is.null(trunc) && trunc$name == "T") {
      lower = if (length(trunc$args) >= 1 && nzchar(trunc$args[1])) trunc$args[1] else "-Inf"
      upper = if (length(trunc$args) >= 2 && nzchar(trunc$args[2])) trunc$args[2] else "Inf"
      if (lower == "-Inf" && upper == "Inf") return("none")
      return(paste0("[", lower, ", ", upper, "]"))
    }
  }
  call = parse_prior_call(parts$distribution)
  if (is.null(call))
    return("none")
  if (call$name == "dunif" && length(call$args) == 2)
    return(paste0("[", call$args[1], ", ", call$args[2], "]"))
  if (call$name == "dbeta") return("[0, 1]")
  if (call$name %in% c("dgamma", "dexp", "dlnorm", "dweib")) return("[0, Inf]")
  "none"
}


prior_kind = function(value, all_names) {
  if (is.numeric(value) || (is.character(value) && grepl("^[-+]?[0-9.eE]+$", trimws(value))))
    return("constant")
  parts = split_prior_truncation(value)
  call = parse_prior_call(parts$distribution)
  if (!is.null(call) && (grepl("^d[A-Za-z]", call$name) || call$name == "dirichlet"))
    return("distribution")
  if (trimws(value) %in% all_names)
    return("alias")
  "expression"
}


compile_prior_record = function(parameter, code, all_names, context, source,
                                 description = NULL) {
  kind = prior_kind(code, all_names)
  original = if (is.numeric(code)) format_prior_number(code) else code
  display_original = translate_legacy_prior(original, context)
  resolved = resolve_prior_value(code, context)

  if (kind == "expression") {
    expr = tryCatch(parse(text = resolved, keep.source = FALSE)[[1]], error = function(e) NULL)
    unknown = if (is.null(expr)) character() else setdiff(all.vars(expr), all_names)
    if (length(unknown) > 0) {
      stop(
        "Unknown model parameter(s) in prior expression for '", parameter, "': ",
        and_collapse(unknown), "."
      )
    }
  }

  rule = if (kind == "distribution") human_prior_call(display_original) else display_original
  if (is.null(description) || is.na(description)) {
    description = switch(
      kind,
      distribution = if (source == "default") "Family-defined default prior" else "User-specified prior",
      alias = paste0("Same value as ", trimws(display_original)),
      expression = "Deterministic expression",
      constant = paste0("Fixed at ", resolved)
    )
  }

  tibble::tibble(
    parameter = parameter,
    prior = if (kind == "distribution") human_prior_call(resolved) else resolved,
    bounds = prior_bounds(display_original),
    rule = rule,
    description = description,
    source = source,
    kind = kind,
    value = resolved
  )
}


compile_prior_specs = function(specs, all_names, context) {
  required = c("parameter", "code", "description", "source")
  assert_data_cols(specs, required)
  records = lapply(seq_len(nrow(specs)), function(i) {
    code = specs$code[[i]]
    if (specs$source[i] == "default")
      code = instantiate_prior_template(code, context)
    compile_prior_record(
      specs$parameter[i], code, all_names, context, specs$source[i],
      specs$description[i]
    )
  })
  dplyr::bind_rows(records)
}
