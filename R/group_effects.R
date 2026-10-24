# ABOUT: Group-level ("random") effects in the predictor part: detecting
# `(... | group)` terms and parsing them into group-level coefficient definitions.
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


# Key that is identical for equivalent group-level terms, e.g., `(x + 1 || id)` and `(1 + x || id)`
group_block_key = function(term) {
  expr = str2lang(term)
  attrs = attributes(stats::terms(stats::as.formula(call("~", expr[[2]]))))
  paste(deparse1(expr[[3]]), attrs$intercept,
    paste(sort(attrs$term.labels), collapse = "+"),
    if (length(attrs$term.labels)) as.character(expr[[1]]) else "", sep = ":")
}


#' Parse one predictor group-level term
#'
#' @param term A term such as `"1 | id"`, `"factor || id"`, or `"0 | id"`.
#' @inheritParams get_coefs_dpar
#' @return A tibble with one row per group-level coefficient, or one row with
#'   `has_coefficients = FALSE` for an explicit turn-off.
#' @noRd
parse_predictor_group_term = function(
  term, segment, dpar, data, par_x, env = parent.frame()
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
  coefficient = get_coefs_dpar(
    data, coefficient_form, segment, dpar, par_x,
    order = NULL,
    design_id = paste("group", dpar, group_col, segment, sep = ":")
  )
  assert_unique_coef_names(coefficient)
  has_coefficients = nrow(coefficient) > 0

  # A turn-off like `(0 | id)`: one row without coefficients
  if (!has_coefficients) {
    return(tibble::tibble(
      dpar = dpar,
      segment = segment,
      group_col = group_col,
      has_coefficients = FALSE,
      population_parameter = NA_character_,
      name = NA_character_,
      sd_name = NA_character_,
      par_type = NA_character_,
      matrix_name = NA_character_,
      order = NA_integer_,
      x_degree = NA_integer_,
      design_id = NA_character_,
      design_col = NA_integer_,
      design_spec = list(NULL),
      matrix_data = list(NULL)
    ))
  }

  coefficient |>
    dplyr::transmute(
      dpar = .data$dpar,
      segment = .data$segment,
      group_col = .env$group_col,
      has_coefficients = TRUE,
      population_parameter = .data$name,
      name = paste0(.data$population_parameter, "_", .env$group_col),
      sd_name = paste0(.data$population_parameter, "_", .env$group_col, "_sd"),
      par_type = .data$par_type,
      matrix_name = .data$matrix_name,
      order = as.integer(.data$order),
      x_degree = .data$x_degree,
      design_id = .data$design_id,
      design_col = .data$design_col,
      design_spec = .data$design_spec,
      matrix_data = .data$matrix_data
    )
}


# Predictor group-level definitions for one segment, including explicit turn-offs
get_predictor_group_definitions_segment = function(
  components, segment, family, data, par_x, previous = NULL
) {
  # Parse group-level terms for all distributional parameters (including mu)
  definitions = list()
  for (dpar in intersect(family$dpar_specs$dpar, names(components))) {
    form_env = environment(components[[dpar]])
    standardized = standardize_same_component(components[[dpar]], segment)
    selected = standardized$selected
    selected = selected[vapply(selected$term_key, is_group_term, logical(1)), , drop = FALSE]

    # Canonical block keys unify equivalent formulas before enforcing one block per group
    terms = get_group_terms(standardized$form)
    blocks = tibble::tibble(term = terms,
      key = vapply(terms, group_block_key, character(1)),
      group = vapply(terms, function(term) deparse1(str2lang(term)[[3]]), character(1)),
      source = selected$source[match(terms, selected$term_key)]) |>
      dplyr::distinct(.data$key, .data$source, .keep_all = TRUE)
    if (anyDuplicated(blocks$group))
      stop("Only one predictor group-level term per distributional parameter and grouping factor is allowed in a segment, including terms in `same()`.", call. = FALSE)

    # Group-level terms are atomic blocks; same() takes their source coding and all coefficients together
    for (i in seq_len(nrow(blocks))) {
      term = blocks$term[i]
      source_segment = blocks$source[i]
      if (!is.na(source_segment)) {
        # `previous` has no columns before the first group-level term
        source = if (NROW(previous) > 0) previous[previous$dpar == dpar & previous$segment == source_segment &
          previous$has_coefficients & previous$group_key == blocks$key[i], , drop = FALSE]
        if (NROW(source) == 0)
          stop("`same((", term, "))` found no identical group-level term for ", dpar, " in segment ", source_segment,
            ". `same()` must repeat a complete group-level term from its source segment.", call. = FALSE)
        block = source
        block$segment = segment
      } else {
        block = parse_predictor_group_term(term, segment, dpar, data, par_x, form_env)
        block$definition_segment = segment
      }
      block$group_key = blocks$key[i]
      definitions[[length(definitions) + 1L]] = block
    }
  }

  dplyr::bind_rows(definitions)
}


# Select group-level parameters by name/part (`pars`) or grouping column (`cols`, e.g., from `facet_by`).
# Each is `NULL`/`FALSE` for none or `TRUE` for all. `pars` also takes `"cp"`, `"predictor"`,
# parameter names, and `*` wildcards.
# Returns list(pars, cols, groups), where `groups` is `pars` split by grouping column.
unpack_group_effects = function(fit, pars = NULL, cols = NULL) {
  checkmate::assert_multi_class(pars, c("logical", "character"), null.ok = TRUE)
  checkmate::assert_multi_class(cols, c("logical", "character"), null.ok = TRUE)
  if (is.logical(pars))
    checkmate::assert_flag(pars)
  if (is.logical(cols))
    checkmate::assert_flag(cols)
  group_pars = mcp_pars(fit, scope = "group")
  use_group = rep(FALSE, nrow(group_pars))

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
      # `*` is a wildcard, e.g., "cp_*" or "state*_1_id"
      is_glob = grepl("*", pars, fixed = TRUE)
      glob_matches = lapply(pars[is_glob], function(glob) grepl(utils::glob2rx(glob), group_pars$parameter))
      known = pars %in% c("cp", "predictor", group_pars$parameter)
      known[is_glob] = vapply(glob_matches, any, logical(1))
      unknown = pars[!known]
      if (length(unknown) > 0)
        stop(
          "Unknown group-effect selection: ", and_collapse(unknown), ". ",
          "Use TRUE, FALSE, \"cp\", \"predictor\", group-level parameter names, or `*` wildcards. ",
          "Available group-level parameters: ", paste0(group_pars$parameter, collapse = ", "), "."
        )
      use_group = group_pars$part %in% pars | group_pars$parameter %in% pars
      use_group = Reduce(`|`, glob_matches, init = use_group)
    }
  } else if (!is.null(cols)) {
    if (all(cols == TRUE)) {
      use_group[] = TRUE
    } else if (!all(cols == FALSE)) {
      use_group = group_pars$group_col %in% cols
    }
  }

  # Return
  list(
    pars = logical0_to_null(group_pars$parameter[use_group]),
    cols = logical0_to_null(group_pars$group_col[use_group]),
    groups = split(group_pars$parameter[use_group], group_pars$group_col[use_group])
  )
}


# Check the hierarchy before drawing a group-level deviation. Change-point deviations are
# drawn from their default bounded normal, so user-specified priors are not supported.
assert_generatable_group_prior = function(fit, name, sd_name) {
  prior_table = fit$.internal$prior_table
  row = match(name, prior_table$parameter)
  prior = prior_table$value[row]
  if (name %in% get_fit_model_tables(fit)$cps$group_name) {
    if (prior_table$source[row] == "user")
      stop(
        "Cannot generate group-level change point '", name, "' from the user-specified prior '", prior_table$rule[row],
        "'. `fit$simulate()` supports only the default prior here. Remove '", name, "' from `mcp(prior = ...)` to simulate.",
        call. = FALSE
      )
    return(invisible(NULL))
  }
  expected = paste0("dnorm(0, ", sd_name, ")")
  if (!identical(prior, expected))
    stop(
      "Cannot generate group-level deviation '", name, "' from prior '", prior,
      "'. This generator supports only '", expected, "' here.",
      call. = FALSE
    )
}
