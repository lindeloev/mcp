# ABOUT: Formula-structure helpers. These parse mcp's segment-formula syntax
# (tildes, response, change point, and right-hand-side sections) without
# building any design matrices.
# ------------------------------------------------------------

#' Takes any formula-like input and returns a formula
#' @aliases to_formula
#' @keywords internal
#' @noRd
#' @param form Formula or character (with or without initial tilde/"~")
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
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


#' Converts formula to string
#'
#' Note: this uses base R and circumvents the length-limitation of `deparse()`
#' and `format()`.
#'
#' @aliases formula_to_char
#' @keywords internal
#' @noRd
#' @param form Any valid formula with any number of tildes.
#' @return A character.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
formula_to_char = function(form) {
  checkmate::assert_formula(form)
  form_char = as.character(form)
  if (length(form_char) == 2 & form_char[1] == "~") {
    return(paste0(form_char, collapse = " "))
  } else if (length(form_char == 3) & form_char[1] == "~") {
    return(paste0(form_char[c(2, 3)], collapse = " ~ "))
  } else {
    stop_github("Could not decode formula ", deparse(form, width.cutoff = 500))
  }
}


#' Returns the right-hand-side of a formula
#'
#' @aliases get_rhs
#' @keywords internal
#' @noRd
#' @param form Formula, e.g. `~x`, `y ~ x` or `y ~ z ~ x`
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs = function(form) {
  checkmate::assert_formula(form)
  if (length(form) == 2) {
    return(form)
  } else if (length(form) == 3) {
    return(form[-2])
  }
}


#' Returns all variables in the predictor parts of an mcpmodel
#'
#' @aliases get_rhs_vars
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return Character vector with unique term names
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs_vars = function(model) {
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")

  model %>%
    lapply(get_rhs) %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()
}


#' Returns grouping-factor variables in predictor group-level terms
#'
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return Character vector of grouping-factor column names.
get_rhs_group_vars = function(model) {
  find_groups = function(expr) {
    if (!is.call(expr))
      return(character())
    if (as.character(expr[[1]]) %in% c("|", "||"))
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

#' Returns all variables in the predictor parts of an mcpmodel
#'
#' @aliases get_model_vars
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return Character vector with unique term names
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_model_vars = function(model) {
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")

  model %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()
}


#' Remove varying or population terms from a formula
#'
#' WARNING: removes response side from the formula
#'
#' @aliases remove_terms
#' @keywords internal
#' @noRd
#' @param form A formula
#' @param remove Either "varying" or "population". These are removed.
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
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
    return(stats::as.formula(formula_str, env=globalenv()))
  }
}


#' Get formula inside a wrapper
#'
#' @aliases get_term_content
#' @keywords internal
#' @noRd
#' @param term E.g., "mu(1 + x)", "sigma(0 + I(x^2))", etc.
#' @return char formula with the content inside the brackets.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_term_content = function(term) {
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
    form = stats::as.formula(paste0("~", content), env = globalenv())
    return(form)
  }
}


#' Takes a formula and returns a string representation of y, cp, and rhs
#' @aliases unpack_tildes
#' @keywords internal
#' @noRd
#' @param form A formula
#' @param i The segment number
#' @return A one-row tibble with columns:
#'   * `form`: String. The full formula for this segment.
#'   * `form_y`: String. The expression for y (without tilde)
#'   * `form_cp`: String. The formula for the change point.
#'   * `form_rhs`: String. The predictor formula. Only used to build the
#'     formula representation in `summary.mcpfit()`.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_tildes = function(form, i) {
  has_LHS = attributes(stats::terms(form))$response == 1
  form_str = formula_to_char(form)
  if (has_LHS == FALSE && i == 1) {
    stop("No response variable in segment 1.")
  } else if (has_LHS == FALSE && i > 1) {
    # If no LHS, add a change point "intercept"
    form_str = paste("1", form_str)
  }

  # List of strings for each section
  chunks = stringr::str_trim(strsplit(form_str, "~")[[1]])

  if (length(chunks) == 2) {
    # Only one tilde. This is the first segment or y is implicit from earlier segment(s)
    return(tibble::tibble(
      form = form_str,
      form_y = ifelse(i == 1, chunks[1], NA),
      form_cp = ifelse(i == 1, NA, paste0(" ~ ", chunks[1])),
      form_rhs = chunks[2]
    ))
  } else if (length(chunks) == 3) {
    if (i == 1)
      stop("The first segment must have exactly one tilde. Got two.")

    return(tibble::tibble(
      form = form_str,
      form_y = chunks[1],
      form_cp = paste0(" ~ ", chunks[2]),
      form_rhs = chunks[3]
    ))
  } else {
    stop("Error in segment ", i, ": Got none or more than two ~ in a segment formula.")
  }
}


#' Unpacks y variable name
#'
#' @aliases unpack_y
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param form_y Character representation of formula
#' @param i Segment number
#' @return A one-row tibble with the response and auxiliary-data columns.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_y = function(form_y, i, family) {
  declared = names(family$response$auxiliary)
  aux_names = unique(c("trials", "weights", declared))
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
    ok_terms = if (length(declared) == 0) rep(FALSE, length(term_labels)) else
      vapply(term_labels, function(term) any(stringr::str_detect(term, paste0("^", declared, "\\("))), logical(1))
    if (!all(ok_terms))
      stop(
        "Only ", if (length(declared) == 0) "no terms are" else and_collapse(paste0("`", declared, "()`")),
        " allowed after the pipe for family = ", family$family, "(). Got '", rhs, "'."
      )
  }

  for (name in declared) {
    term_index = stringr::str_detect(term_labels, paste0("^", name, "\\("))
    got_term = any(term_index)
    if (family$response$auxiliary[[name]]$required && !got_term)
      stop("Error in response of segment ", i, ": need a valid ", name, "() specification.")
    if (!got_term)
      next
    if (sum(term_index) > 1)
      stop("Only one ", name, "() term is allowed in segment ", i, ".")

    term = term_labels[term_index]
    content = get_term_content(term)
    column = attr(stats::terms(content), "term.labels")
    if (length(column) != 1)
      stop("There must be exactly one term inside ", name, "(). Got ", term, " in segment ", i)
    response[[name]] = column
  }

  tibble::as_tibble(response)
}


#' Takes a cp formula (as a string) and returns its properties
#'
#' @aliases unpack_cp
#' @keywords internal
#' @noRd
#' @param form_cp Segment formula as string.
#' @param i segment number
#' @return A one-row tibble with columns:
#'   * `cp_int`: bool. Whether there is an intercept change in the change point.
#'   * `cp_varying`: bool or NA. Is there a group-level intercept on the change point?
#'   * `cp_group_col`: char or NA. Which data column defines the grouping factor?
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_cp = function(form_cp, i) {
  if (is.na(form_cp)) {
    return(tibble::tibble(
      cp_int = FALSE,
      cp_varying = FALSE,
      cp_group_col = NA
    ))
  }
  form_cp = stats::as.formula(form_cp)

  # Group-level effects
  form_varying = remove_terms(form_cp, "population")

  if (!is.null(form_varying)) {
    varying_terms = attr(stats::terms(form_varying), "term.labels")
    if (length(varying_terms) > 1)
      stop("Error in segment", i, " (change point): only one group-level effect is allowed. Found ", form_cp)

    varying_parts = strsplit(gsub(" ", "", varying_terms), "\\|")[[1]]
    if (!varying_parts[1] == "1")
      stop("Error in segment ", i, " (change point): Only plain intercepts are allowed in group-level effects, e.g., (1|id).", i)

    if (!grepl("^[A-Za-z._0-9]+$", varying_parts[2]))
      stop("Error in segment ", i, " (change point): invalid grouping-variable format in group-level effect. Got: ", varying_parts[2])
  }

  # Population-level effects
  attrs = attributes(stats::terms(remove_terms(form_cp, "varying")))
  if (length(attrs$term.labels) > 0)
    stop("Error in segment ", i, " (change point): Only intercepts (1) are allowed in population-level effects.")

  if (is.null(form_varying) && attrs$intercept == 0)
    stop("Error in segment ", i, " (change point): no population-level intercept or group-level effect. You can do e.g., ~ 1 or ~ (1 |id).")

  # Return as list.
  if (!is.null(form_varying)) {
    # If there is a group-level effect
    return(tibble::tibble(
      cp_int = attrs$intercept == 1,
      cp_varying = ifelse(varying_parts[1] == "1", TRUE, NA),  # placeholder for later
      cp_group_col = varying_parts[2]
    ))
  } else {
    # If there is no group-level effect
    return(tibble::tibble(
      cp_int = attrs$intercept == 1,
      cp_varying = FALSE,
      cp_group_col = NA
    ))
  }
}
