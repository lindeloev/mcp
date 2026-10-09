# ABOUT: Formula-structure helpers. These parse mcp's segment-formula syntax
# (tildes, response, change point, and right-hand-side sections) without
# building any design matrices.
# ------------------------------------------------------------

# Takes any formula-like input (formula or string) and returns a formula
to_formula = function(form) {
  checkmate::assert(
    checkmate::check_character(form, min.len = 1, max.len = 3),
    checkmate::check_formula(form),
    .var.name = "form"
  )
  if (is.character(form)) {
    # Add tilde
    if (!stringr::str_detect(form, "^(\\s|)~")) {
      form = paste0("~", form)
    }
    form = stats::as.formula(form)
  }

  form
}


# Convert formula to string using deparse1() to guarantee a single string.
# One-sided formulas get a space after the tilde ("~ 0 + x"), like two-sided ones.
formula_to_char = function(form) {
  checkmate::assert_formula(form)
  sub("^~", "~ ", deparse1(form))
}


# Returns the right-hand-side of a formula, excluding the change point in cp()
get_rhs = function(form) {
  checkmate::assert_formula(form)
  rhs = split_cp(form[[length(form)]])$rhs
  stats::as.formula(call("~", rhs), env = environment(form))
}


# Check whether an expression is a formula call, e.g., `y ~ x`
is_tilde_call = function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("~"))
}


# Check whether an expression is a call to cp()
is_cp_call = function(expr) {
  is.call(expr) && identical(deparse1(expr[[1]]), "cp")
}


# Recursively check whether an expression contains any call to cp()
contains_cp_call = function(expr) {
  if (is_cp_call(expr))
    return(TRUE)
  if (!is.call(expr))
    return(FALSE)

  any(vapply(as.list(expr)[-1], contains_cp_call, logical(1)))
}


# Split a segment RHS expression into the cp() content (NULL if absent) and the remaining terms
# E.g., `cp(1 + (1 | id)) + 0 + x` --> list(cp = `1 + (1 | id)`, rhs = `0 + x`).
split_cp = function(rhs) {
  leaves = unpack_additive(rhs)
  is_cp = vapply(leaves, is_cp_call, logical(1))
  if (!any(is_cp))
    return(list(cp = NULL, rhs = rhs))
  if (sum(is_cp) > 1)
    stop("Only one `cp()` allowed in each segment formula. Got `", deparse1(rhs), "`.", call. = FALSE)

  cp = leaves[is_cp][[1]]
  if (length(cp) != 2)
    stop("`cp()` takes exactly one formula, e.g., `cp(1 + (1 | id))`. Got `", deparse1(cp), "`.", call. = FALSE)

  rest = leaves[!is_cp]
  list(
    cp = cp[[2]],
    rhs = if (length(rest) == 0) 1 else Reduce(function(a, b) call("+", a, b), rest)
  )
}




# Returns all variables in the predictor parts of an mcpmodel
get_rhs_vars = function(model) {
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")

  vars = model %>%
    lapply(get_rhs) %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()

  unique(c(vars, get_arma_series(model)))
}


# Returns grouping-factor variables in predictor group-level terms
get_rhs_group_vars = function(model) {
  find_groups = function(expr) {
    if (!is.call(expr))
      return(character())
    if (as.character(expr[[1]])[1] %in% c("|", "||"))
      return(all.vars(expr[[3]]))
    nested_groups = lapply(rlang::call_args(expr), find_groups)
    unique(unlist(nested_groups))
  }

  model %>%
    lapply(get_rhs) %>%
    lapply(function(form) find_groups(form[[2]])) %>%
    unlist() %>%
    unique()
}


# Returns variables appearing inside offset() calls
get_rhs_offset_vars = function(model) {
  find_offsets = function(expr) {
    if (!is.call(expr))
      return(character())
    if (deparse1(expr[[1]]) %in% c("offset", "stats::offset"))
      return(all.vars(expr))
    unique(unlist(lapply(rlang::call_args(expr), find_offsets)))
  }

  unique(unlist(lapply(model, function(m) find_offsets(get_rhs(m)[[2]]))))
}


# Returns variables appearing in RHS outside offset() calls
get_rhs_non_offset_vars = function(model) {
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")

  find_non_offset = function(expr) {
    if (!is.call(expr))
      return(all.vars(expr))
    if (deparse1(expr[[1]]) %in% c("offset", "stats::offset"))
      return(character())
    if (deparse1(expr[[1]]) %in% c("|", "||"))
      return(find_non_offset(expr[[2]]))
    unique(unlist(lapply(rlang::call_args(expr), find_non_offset)))
  }

  vars = model %>%
    lapply(get_rhs) %>%
    lapply(function(form) find_non_offset(form[[2]])) %>%
    unlist() %>%
    unique()

  unique(c(vars, get_arma_series(model)))
}

# Returns all variables in the predictor parts of an mcpmodel
get_model_vars = function(model) {
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")

  vars = model %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()

  unique(c(vars, get_arma_series(model)))
}

#' Remove varying or population terms from a formula
#'
#' WARNING: removes response side from the formula
#'
#' @param form A formula
#' @param remove Either "varying" or "population". These are removed.
#' @return A formula
#' @noRd
remove_terms = function(form, remove) {
  checkmate::assert_formula(form)
  remove = rlang::arg_match0(remove, c("varying", "population"))

  # Find terms with "|"
  attrs = attributes(stats::terms(form))
  term.labels = attrs$term.labels
  varying_bool = stringr::str_detect(term.labels, "\\|")

  # Add parenthesis back to them
  term.labels[varying_bool] = paste0("(", term.labels[varying_bool], ")")

  # Remove non-matching types
  if (remove == "varying") {
    term.labels = term.labels[!varying_bool]
    # base::terms() omits offset() from term.labels; re-attach it for population formulas
    if (!is.null(attrs$offset)) {
      offset_terms = vapply(attrs$offset, function(i) deparse1(attrs$variables[[i + 1]]), character(1))
      term.labels = c(term.labels, offset_terms)
    }
    term.labels = c(attrs$intercept, term.labels)  # Add intercept indicator
  } else if (remove == "population") {
    term.labels = term.labels[varying_bool]
  }

  # Build formula from terms and return
  if (length(term.labels) == 0) {
    return(NULL)
  } else {
    formula_terms = paste0(term.labels, collapse = " + ")
    formula_str = paste0("~", formula_terms)
    # Rebuild in the caller's environment so local transformations still resolve.
    return(stats::as.formula(formula_str, env = environment(form)))
  }
}


#' Get formula inside a wrapper
#'
#' @param term E.g., "mu(1 + x)", "sigma(0 + I(x^2))", etc.
#' @param env Environment in which the term was originally defined.
#' @return char formula with the content inside the brackets.
#' @noRd
get_term_content = function(term, env = parent.frame()) {
  # Handle cases of no input or several inputs
  if (length(term) == 0) {
    return(NA)
  } else if (length(term) > 1) {
    stop("Only one ", term, " allowed in each formula.")
  } else if (is.na(term)) {
    return(NA)
  } else if (length(term) == 1) {
    # Get formula inside wrapper
    content_start = stringr::str_locate(term, "\\(") + 1  # Location of first character in contents
    content_end = stringr::str_length(term) - 1  # Location of last character in contents
    content = substr(term, content_start, content_end)

    # To formula
    if (content == "")
      stop("Empty terms not allowed in the formulas. Found '", term, "'.")
    # Wrapper terms are strings, so carry their source formula environment explicitly.
    form = stats::as.formula(paste0("~", content), env = env)
    return(form)
  }
}


# Unpack an additive expression into a list of leaf terms
unpack_additive = function(expr) {
  if (is.call(expr) && identical(deparse1(expr[[1]]), "+")) {
    c(unpack_additive(expr[[2]]), unpack_additive(expr[[3]]))
  } else {
    list(expr)
  }
}


# Canonicalize a segment RHS formula so that bare mu terms are wrapped in mu(...)
# E.g., `~ 1 + x + sigma(1 + x)` --> `~ mu(1 + x) + sigma(1 + x)`.
# This allows treating mu like any other wrapper, simplifying downstream parsing.
canonicalize_rhs = function(form_rhs, family) {
  checkmate::assert_formula(form_rhs)
  checkmate::assert_true(is.mcpfamily(family), .var.name = "family")
  env = environment(form_rhs)

  # Unpack top-level additive leaves and their call names
  leaves = unpack_additive(form_rhs[[2]])
  leaf_names = vapply(leaves, function(x) if (is.call(x)) deparse1(x[[1]]) else "", character(1))

  # Separate component wrappers, explicit mu() calls, and bare terms
  dpar_wrappers = setdiff(unique(c(known_dpar_wrappers(), family$dpar_specs$dpar)), "mu")
  known_wrappers = c(dpar_wrappers, "ar", "ma")

  wrapper_leaves = leaves[leaf_names %in% known_wrappers]
  explicit_mu = leaves[leaf_names == "mu"]
  bare_leaves = leaves[leaf_names %notin% c(known_wrappers, "mu")]

  if (length(explicit_mu) > 1)
    stop("Only one `mu()` allowed in each formula.", call. = FALSE)

  # Merge explicit mu() contents with any bare terms
  inner_mu = if (length(explicit_mu) == 1 && length(explicit_mu[[1]]) > 1) {
    unpack_additive(explicit_mu[[1]][[2]])
  } else {
    list()
  }
  mu_leaves = c(inner_mu, bare_leaves)

  # Helper to chain additive expressions
  sum_exprs = function(exprs) Reduce(function(a, b) call("+", a, b), exprs)

  # Wrap mu terms into mu(...) (defaulting to mu(1) if empty)
  mu_expr = if (length(mu_leaves) == 0) call("mu", 1) else call("mu", sum_exprs(mu_leaves))

  # Return canonical one-sided formula
  stats::as.formula(call("~", sum_exprs(c(list(mu_expr), wrapper_leaves))), env = env)
}


# Canonicalize a segment formula into a uniform AST structure
# Returns a list: response (expr), cp (expr or NULL for segment 1), rhs (expr), and form (canonical formula)
# E.g., `~ cp((1 | id)) + 0 + x` gives cp = `(1 | id)` and rhs = `0 + x`.
canonicalize_segment = function(form, i, default_response = NULL) {
  form = to_formula(form)
  if (i > 1 && length(form) == 3) {
    if (identical(form[[2]], default_response))
      stop("Error in segment ", i, ": segments 2+ use the response from segment 1. Remove `", deparse1(form[[2]]), " ~`, e.g., `~ 0 + x`.", call. = FALSE)
    form = translate_cp_lhs(form, i, default_response)
  }

  split = split_cp(form[[length(form)]])
  if (contains_cp_call(split$rhs))
    stop("Error in segment ", i, ": `cp()` must be a separate term in the formula, e.g., `~ cp((1 | id)) + 0 + x`. Got `", deparse1(form), "`.", call. = FALSE)

  # Segment 1 must define the response and cannot contain a change point
  if (i == 1) {
    if (length(form) == 2)
      stop("No response variable in segment 1.", call. = FALSE)
    if (is_tilde_call(form[[2]]))
      stop("The first segment must have exactly one tilde. Got two.", call. = FALSE)
    if (!is.null(split$cp))
      stop("Segment 1 has no change point, so `cp()` is only allowed in segments 2+.", call. = FALSE)

    return(list(response = form[[2]], cp = NULL, rhs = split$rhs, form = form))
  }

  # Segment > 1: the response is inherited from segment 1 and the change point is in cp(), defaulting to cp(1)
  cp = if (is.null(split$cp)) 1 else split$cp
  list(response = default_response, cp = cp, rhs = split$rhs, form = form)
}


#' Unpacks y variable name
#'
#' @inheritParams mcp
#' @param form_y Character representation of formula
#' @param i Segment number
#' @return A one-row tibble with the response and auxiliary-data columns.
#' @noRd
unpack_y = function(form_y, i, family, env = parent.frame()) {
  if (is.language(form_y))
    form_y = deparse1(form_y)

  family_aux = names(family$response$auxiliary)
  aux_names = unique(c("trials", "weights", family_aux))
  response = stats::setNames(as.list(rep(NA_character_, length(aux_names) + 1)), c("y", aux_names))

  # If NA and not segment 1, just return empty
  if (is.na(form_y)) {
    if (i == 1)
      stop("A response must be defined in segment 1, e.g., 'y ~ 1'")

    return(tibble::as_tibble(response))
  }


  # Split by |
  y_split = strsplit(form_y, "\\|")[[1]]
  if (length(y_split) > 2)
    stop("There can only be zero or one pipe in response. Got '", form_y, "' in segment ", i)

  # RESPONSE
  lhs = y_split[1]
  y_col = attr(stats::terms(to_formula(lhs)), "term.labels")
  if (length(y_col) != 1)
    stop("There should be exactly one response variable. Got ", length(y_col), " in segment ", i)
  response$y = y_col

  term_labels = character()
  if (length(y_split) == 2) {
    rhs = y_split[2]
    term_labels = attr(stats::terms(to_formula(rhs)), "term.labels")
    ok_terms = if (length(family_aux) == 0) rep(FALSE, length(term_labels)) else
      vapply(term_labels, function(term) any(stringr::str_detect(term, paste0("^", family_aux, "\\("))), logical(1))
    if (!all(ok_terms))
      stop(
        "Only ", if (length(family_aux) == 0) "no terms are" else and_collapse(paste0("`", family_aux, "()`")),
        " allowed after the pipe for family = ", family$family, "(). Got '", rhs, "'."
      )
  }

  for (name in family_aux) {
    term_index = stringr::str_detect(term_labels, paste0("^", name, "\\("))
    got_term = any(term_index)
    if (family$response$auxiliary[[name]]$required && !got_term)
      stop("Error in response of segment ", i, ": need a valid ", name, "() specification.")
    if (!got_term)
      next
    if (sum(term_index) > 1)
      stop("Only one ", name, "() term is allowed in segment ", i, ".")

    term = term_labels[term_index]
    content = get_term_content(term, env)
    column = attr(stats::terms(content), "term.labels")
    if (length(column) != 1)
      stop("There must be exactly one term inside ", name, "(). Got ", term, " in segment ", i)
    response[[name]] = column
  }

  tibble::as_tibble(response)
}


#' Takes a cp formula (as a string) and returns its properties
#'
#' @param form_cp Segment formula as string.
#' @param i segment number
#' @return A one-row tibble with columns:
#'   * `cp_varying`: Logical scalar. Is there a group-level intercept on the change point?
#'   * `cp_group_col`: char or NA. Which data column defines the grouping factor?
#' @noRd
unpack_cp = function(form_cp, i, env = parent.frame()) {
  # Segment 1 has no change point
  if (is.null(form_cp)) {
    return(tibble::tibble(
      cp_varying = FALSE,
      cp_group_col = NA
    ))
  }

  form_cp = if (rlang::is_formula(form_cp)) form_cp else stats::as.formula(call("~", form_cp), env = env)

  # Group-level effects
  form_varying = remove_terms(form_cp, "population")

  if (!is.null(form_varying)) {
    varying_terms = attr(stats::terms(form_varying), "term.labels")
    if (length(varying_terms) > 1)
      stop("Error in segment ", i, " (change point): only one group-level effect is allowed. Found ", deparse1(form_cp))

    varying_parts = strsplit(gsub(" ", "", varying_terms), "\\|\\|?")[[1]]  # (1 || id) is (1 | id) for intercepts
    if (!varying_parts[1] == "1")
      stop("Error in segment ", i, " (change point): Only plain intercepts are allowed in group-level effects, e.g., `cp((1 | id))`.")

    if (!grepl("^[A-Za-z._0-9]+$", varying_parts[2]))
      stop("Error in segment ", i, " (change point): invalid grouping-variable format in group-level effect. Got: ", varying_parts[2])
  }

  # Population-level effects
  attrs = attributes(stats::terms(remove_terms(form_cp, "varying")))
  if (length(attrs$term.labels) > 0)
    stop("Error in segment ", i, " (change point): Only intercepts are allowed in population-level effects, e.g., `cp(1 + (1 | id))`.")

  if (attrs$intercept == 0)
    stop("Error in segment ", i, " (change point): change points always include a population-level intercept. Remove `0 +`, e.g., `cp((1 | id))`.")

  tibble::tibble(
    cp_varying = !is.null(form_varying),
    cp_group_col = if (is.null(form_varying)) NA else varying_parts[2]
  )
}
