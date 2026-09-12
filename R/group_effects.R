# ABOUT: Group-level ("varying"/"random") predictor effects: detecting
# `(... | group)` terms in a predictor formula and parsing them into
# group-effect coefficient definitions.
# ------------------------------------------------------------

# Is a formula term a group-level term?
is_group_term = function(term) {
  expr = str2lang(term)
  is.call(expr) && as.character(expr[[1]])[1] %in% c("|", "||")
}


# Extract group-level terms from a one-sided formula
get_group_terms = function(form) {
  checkmate::assert_formula(form)
  terms = attr(stats::terms(form), "term.labels")
  terms[vapply(terms, is_group_term, logical(1))]
}


#' Parse one predictor group-level term
#'
#' @keywords internal
#' @noRd
#' @param term A term such as `"1 | id"`, `"factor || id"`, or `"0 | id"`.
#' @inheritParams get_predictors_dpar
#' @return A tibble with one row per group-level coefficient, or one inactive
#'   row for an explicit turn-off.
parse_predictor_group_term = function(
  term, segment, dpar, data, par_x, check_rank = TRUE, env = parent.frame()
) {
  expr = str2lang(term)
  operator = as.character(expr[[1]])
  if (!is.call(expr) || operator %notin% c("|", "||") || length(expr) != 3)
    stop_github("Expected a group-level term, got '", term, "'.")

  group_expr = expr[[3]]
  if (!is.symbol(group_expr))
    stop(
      "Grouping factors in predictor group-level terms must be plain data-column names. Got `",
      paste(deparse(group_expr), collapse = ""), "` in segment ", segment, "."
    )
  group_col = as.character(group_expr)
  if (group_col %notin% names(data))
    stop("Grouping factor '", group_col, "' was not found in data.")
  if (anyNA(data[[group_col]]))
    stop("Grouping factor '", group_col, "' contains missing values.")
  if (is.numeric(data[[group_col]]) && any(data[[group_col]] != floor(data[[group_col]])))
    stop("Grouping factor '", group_col, "' must be integer, character, or factor.")

  # Group-level coefficient terms inherit the environment of their source formula.
  coefficient_form = stats::as.formula(call("~", expr[[2]]), env = env)
  attrs = attributes(stats::terms(coefficient_form))
  if (operator == "|" && length(attrs$term.labels) > 0) {
    stop(
      "Predictor group-level slopes and factor terms currently require `||`, ",
      "for example `(1 + x || id)`. Correlated `|` terms are not yet supported. ",
      "Got `(", term, ")` in segment ", segment, "."
    )
  }
  coefficient = get_predictors_dpar(
    data, coefficient_form, segment, dpar, par_x,
    order = NULL, check_rank = check_rank,
    design_id = paste("group", dpar, group_col, segment, sep = ":")
  )
  assert_unique_predictor_names(coefficient)
  active = nrow(coefficient) > 0

  if (!active) {
    return(tibble::tibble(
      dpar = dpar,
      segment = segment,
      group_col = group_col,
      group_term = term,
      active = FALSE,
      correlated = operator == "|",
      population_name = NA_character_,
      name = NA_character_,
      sd_name = NA_character_,
      par_type = NA_character_,
      matrix_name = NA_character_,
      display_name = NA_character_,
      order = NA_integer_,
      x_factor = NA_character_,
      design_id = NA_character_,
      design_col = NA_integer_,
      design_spec = list(NULL),
      matrix_data = list(NULL)
    ))
  }

  coefficient %>%
    dplyr::transmute(
      dpar = .data$dpar,
      segment = .data$segment,
      group_col = .env$group_col,
      group_term = .env$term,
      active = TRUE,
      correlated = .env$operator == "|",
      population_name = .data$code_name,
      name = paste0(.data$code_name, "_", .env$group_col),
      sd_name = paste0(.data$code_name, "_", .env$group_col, "_sd"),
      par_type = .data$par_type,
      matrix_name = .data$matrix_name,
      display_name = .data$display_name,
      order = as.integer(.data$order),
      x_factor = .data$x_factor,
      design_id = .data$design_id,
      design_col = .data$design_col,
      design_spec = .data$design_spec,
      matrix_data = .data$matrix_data
    )
}


#' Get predictor group-level definitions for one segment
#'
#' @keywords internal
#' @noRd
#' @inheritParams get_predictors_segment
#' @return A tibble containing active definitions and explicit turn-offs.
get_predictor_group_definitions_segment = function(
  form_rhs, segment, family, data, par_x, check_rank = TRUE
) {
  form_rhs = stats::as.formula(form_rhs)
  form_env = environment(form_rhs)
  term_labels = attr(stats::terms(form_rhs), "term.labels")
  top_level = term_labels[vapply(term_labels, is_group_term, logical(1))]

  definitions = lapply(
    top_level,
    parse_predictor_group_term,
    segment = segment,
    dpar = "mu",
    data = data,
    par_x = par_x,
    check_rank = check_rank,
    env = form_env
  )

  for (dpar in setdiff(family$dpar_specs$dpar, "mu")) {
    dpar_terms = term_labels[stringr::str_detect(term_labels, paste0("^", dpar, "\\("))]
    if (length(dpar_terms) == 0)
      next
    dpar_form = get_term_content(dpar_terms, form_env)
    definitions = c(
      definitions,
      lapply(
        get_group_terms(dpar_form),
        parse_predictor_group_term,
        segment = segment,
        dpar = dpar,
        data = data,
        par_x = par_x,
        check_rank = check_rank,
        env = form_env
      )
    )
  }

  definitions = dplyr::bind_rows(definitions)
  if (nrow(definitions) == 0)
    return(definitions)

  group_terms = definitions %>%
    dplyr::distinct(
      .data$dpar, .data$group_col, .data$segment, .data$group_term
    )
  keys = paste(group_terms$dpar, group_terms$group_col, group_terms$segment)
  if (anyDuplicated(keys)) {
    duplicated_keys = unique(keys[duplicated(keys) | duplicated(keys, fromLast = TRUE)])
    stop(
      "Only one predictor group-level term per distributional parameter and grouping factor is allowed in a segment. Found ",
      and_collapse(duplicated_keys), " in segment ", segment, "."
    )
  }
  definitions
}


# Get information about group-level parameters
#
# Returns parameters, data columns, and effect metadata given parameter
# name(s), model part(s), or column(s).
#
# - pars: `NULL`/`FALSE` for nothing. `TRUE` for all. A character vector
#   containing `"cp"`, `"predictor"`, or exact group-level parameter names.
# - cols: `NULL`/`FALSE` for nothing. `TRUE` for all. A vector of grouping
#   column names for specifics. Usually provided via `facet_by` elsewhere.
# Returns: A list. See details.
#   
# Returns a list with
# * `pars`: Character vector of parameter names, or `NULL` if empty.
# * `cols`: Character vector of data column names, or `NULL` if empty.
# * `indices`: Logical vector indexing the group-effects table.
# * `effects`: The selected rows of the group-effects table.
unpack_group_effects = function(fit, pars = NULL, cols = NULL) {
  checkmate::assert_multi_class(pars, c("logical", "character"), null.ok = TRUE)
  checkmate::assert_multi_class(cols, c("logical", "character"), null.ok = TRUE)
  if (is.logical(pars))
    checkmate::assert_flag(pars)
  if (is.logical(cols))
    checkmate::assert_flag(cols)
  group_effects = get_fit_model_tables(fit)$group_effects
  use_group = rep(FALSE, nrow(group_effects))

  if (!is.null(pars) && !is.null(cols)) {
    stop("One of `pars` and `cols` must be NULL.")
  }


  if (!is.null(pars)) {
    if (all(pars == FALSE)) {
      # Select no group-level effects
      use_group[] = FALSE
    } else if (all(pars == TRUE)) {
      # Select all group-level effects
      use_group[] = TRUE
    } else if (is.character(pars)) {
      allowed = c("cp", "predictor", group_effects$name)
      unknown = pars[pars %notin% allowed]
      if (length(unknown) > 0)
        stop(
          "Unknown group-effect selection: ", and_collapse(unknown), ". ",
          "Use TRUE, FALSE, \"cp\", \"predictor\", or a group-effect name."
        )
      use_group = group_effects$part %in% pars | group_effects$name %in% pars
    }
  } else if (!is.null(cols)) {
    if (all(cols == TRUE)) {
      use_group[] = TRUE
    } else if (!all(cols == FALSE)) {
      use_group = group_effects$group_col %in% cols
    }
  }

  # Return
  list(
    pars = logical0_to_null(group_effects$name[use_group]),
    cols = logical0_to_null(group_effects$group_col[use_group]),
    indices = use_group,
    effects = group_effects[use_group, , drop = FALSE]
  )
}
