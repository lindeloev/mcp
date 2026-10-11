# ABOUT: AR/MA-specific code: formula parsing, GARMA recurrence generation,
# R-side simulation/evaluation, and root-condition (stationarity/invertibility)
# checks. In prose, "AR/MA" refers to model terms and effects, while "GARMA"
# refers to the generalized link-scale recurrence.
# ------------------------------------------------------------

#' Unpack an ar() or ma() term
#'
#' @param expr A call such as `ar(2)`, `ma(1, 1 + x)`, or `ar(1, threshold = 0.2, series = id)`.
#' @return A list with `$order`, `$formula` (the expression inside the term; `1`
#'   if no formula is given and `0` if `order = 0`), `$threshold`, and `$series`.
#' @noRd
unpack_arma = function(expr) {
  term = deparse1(expr)
  component = rlang::arg_match0(as.character(expr[[1]]), c("ar", "ma"))

  component_args = as.list(expr)[-1]
  component_arg_names = names(component_args)
  if (is.null(component_arg_names))
    component_arg_names = rep("", length(component_args))

  if (length(component_args) == 0 || component_arg_names[1] %notin% c("", "order"))
    stop("The first argument to ", component, "() must be its order.")

  threshold_index = which(component_arg_names == "threshold")
  if (length(threshold_index) > 1)
    stop("Only one `threshold` value is allowed in ", component, "().")

  series_index = which(component_arg_names == "series")
  if (length(series_index) > 1)
    stop("Only one `series` value is allowed in ", component, "().")

  formula_index = setdiff(seq_along(component_args), c(1, threshold_index, series_index))
  if (length(formula_index) > 1 || any(component_arg_names[formula_index] %notin% c("", "formula")))
    stop(component, "() accepts only `order`, an optional formula, `threshold`, and `series`.")

  # GET ORDER
  order = suppressWarnings(as.numeric(deparse1(component_args[[1]])))

  # Check the order
  if (is.na(order))
    stop("Wrong specification of order in '", term, "'. Must be ", component, "(order) or ", component, "(order, formula) where order is a non-negative integer.")
  checkmate::assert_int(order, lower = 0, .var.name = term)

  # GET FORMULA AND THRESHOLD AND SERIES
  if (order == 0) {
    if (length(formula_index) == 1)
      stop("Formula cannot be specified when order = 0 in '", term, "'.")
    formula = 0
  } else if (length(formula_index) == 1) {
    formula = component_args[[formula_index]]
  } else {
    # If there is no formula, use an intercept-only component formula.
    formula = 1
  }

  if (length(threshold_index) == 1) {
    threshold = suppressWarnings(as.numeric(deparse1(component_args[[threshold_index]])))
    if (length(threshold) != 1 || is.na(threshold) || !is.finite(threshold) || threshold <= 0 || threshold >= 1)
      stop("`threshold` in ", component, "() must be one number between 0 and 1.")
  } else {
    # Defer the default until AR and MA have been parsed together for a segment.
    threshold = NA_real_
  }

  if (length(series_index) == 1) {
    series_arg = component_args[[series_index]]
    if (is.symbol(series_arg) || is.character(series_arg)) {
      series = as.character(series_arg)
    } else {
      stop("`series` in ", component, "() must identify a data column name.")
    }
    if (length(series) != 1 || !nzchar(series))
      stop("`series` in ", component, "() must identify a data column name.")
  } else {
    series = NA_character_
  }

  # Return
  list(
    order = order,
    formula = formula,
    threshold = threshold,
    series = series
  )
}


# Returns the requested AR/MA order, or NA if the term is absent
get_arma_order = function(predictors, term) {
  term = rlang::arg_match0(term, c("ar", "ma"))
  orders = predictors$order[predictors$dpar == term]
  if (length(orders) == 0) NA else max(orders, na.rm = TRUE)
}


# The ar()/ma() terms included in each segment, given `get_segment_components()` for each segment.
# Positive-order ~ 0 formulas have no coefficient rows but are joined to the earlier
# intercept and x-terms.
get_included_arma = function(components) {
  rows = list()
  for (segment in seq_along(components)) {
    for (component in intersect(c("ar", "ma"), names(components[[segment]])))
      rows[[length(rows) + 1L]] = tibble::tibble(dpar = component, segment = segment, order = components[[segment]][[component]]$order)
  }
  dplyr::bind_rows(rows)
}


# Validate an optional independent-series column
assert_arma_series = function(data, series) {
  checkmate::assert_string(series, null.ok = TRUE)
  if (is.null(series))
    return(invisible(NULL))

  assert_data_cols(data, series, fail_funcs = c(is.na))
  if (!is.atomic(data[[series]]))
    stop("`series` must identify an atomic data column.")
  runs = rle(match(data[[series]], unique(data[[series]])))$values
  if (anyDuplicated(runs))
    stop("Rows belonging to each `series` must be contiguous. Sort the data by `series` and then by time.")

  invisible(NULL)
}


# Validate GARMA observation thresholds against family-appropriate bounds
assert_arma_thresholds = function(family, threshold, data = list()) {
  threshold = unique(stats::na.omit(threshold))
  if (length(threshold) == 0)
    return(invisible(TRUE))

  if (!is.null(family$garma$validate_threshold)) {
    family$garma$validate_threshold(threshold, data)
  } else {
    if (any(threshold <= 0 | threshold >= 1))
      stop("`threshold` in ar() / ma() must be strictly between 0 and 1.")
  }

  invisible(TRUE)
}


# Check whether this model has any AR/MA terms
is_arma = function(fit) {
  any(get_fit_model_tables(fit)$predictors$dpar %in% c("ar", "ma"))
}


# Check whether missing responses enter a later observed GARMA history
has_missing_garma_history = function(fit, data = fit$data) {
  # Is not AR/MA
  if (!is_arma(fit))
    return(FALSE)

  # No missing values in columns
  columns = mcp_columns(fit)
  y = data[[columns$response]]
  if (is.null(y) || !anyNA(y))
    return(FALSE)

  # Missing values in the history of later observed values (per series)
  # (e.g., c(1, NA, 3) --> TRUE; c(1, 3, NA, NA) --> FALSE
  series = if (is.null(columns$series) || is.null(data[[columns$series]])) rep(1, length(y)) else data[[columns$series]]
  garma_missing_vec = vapply(split(y, series), function(series_y) {
    missing_rows = which(is.na(series_y))
    observed_rows = which(!is.na(series_y))

    length(missing_rows) > 0 &&
      length(observed_rows) > 0 &&
      min(missing_rows) < max(observed_rows)
  }, logical(1))
  any(garma_missing_vec)
}


# Reject observed-data likelihoods with unresolved GARMA histories
# - arma: Whether AR/MA effects are being evaluated (logical).
# - what: String describing the operation, used in the error message.
assert_loglik_garma_history = function(fit, data = fit$data, arma = TRUE,
                                        what = "Log-likelihood evaluation") {
  if (!arma || !has_missing_garma_history(fit, data))
    return(invisible(NULL))

  stop(
    what, " is unavailable with `arma = TRUE` because a missing response occurs ",
    "in the history of a later observed response. Correct evaluation requires ",
    "integrating over the missing GARMA history, which mcp does not currently ",
    "implement. `arma = FALSE` evaluates the model without its AR/MA contribution ",
    "and should only be used if that is intentional.",
    call. = FALSE
  )
}


# Match retained JAGS imputations to expanded posterior evaluation rows
get_imputed_response_draws = function(fit, draws) {
  # Initialize one optional imputation for each expanded evaluation row
  imputed = fit$.internal$imputed_response
  rows = fit$.internal$imputed_response_rows
  values = rep(NA_real_, nrow(draws))

  # Return missing values when the fit contains no retained imputations
  if (is.null(imputed) || length(rows) == 0)
    return(values)

  # Convert stored chains and align them with the requested posterior draws
  imputed = posterior::as_draws_df(posterior::as_draws_array(imputed))
  draw_index = match(draws$.draw, imputed$.draw)
  response = mcp_columns(fit)$response

  # Copy each monitored response node to its matching data row and draw
  for (row in rows) {
    use = draws$.data_row == row
    node = paste0(response, "[", row, "]")
    values[use] = imputed[[node]][draw_index[use]]
  }

  values
}


# Warn about model-checking limitations for AR/MA models
# - fit: An mcpfit object.
# - arma: Whether AR and MA effects are used in the evaluation (logical).
# - check: Currently only "information_criterion".
warn_arma_check = function(fit, arma, check) {
  if (!arma || !is_arma(fit))
    return(invisible(NULL))

  warning(
    switch(
      check,
      information_criterion = paste(
        "Observationwise PSIS-LOO/WAIC is problematic for AR/MA models because both treat",
        "individual conditional likelihood terms as validation units. In PSIS-LOO, a held-out",
        "response also remains in the conditioning history of later terms. Prefer leave-future-out",
        "or blocked cross-validation. These methods are not currently implemented in mcp."
      ),
      stop_github("Unknown AR/MA check type: ", check)
    ),
    call. = FALSE
  )
  invisible(NULL)
}


# Check AR stationarity or MA invertibility row by row
# - values: Evaluated model parameters, including ar1_, ma1_, etc.
# - component: Either "ar" or "ma".
arma_root_violations = function(values, component) {
  component = rlang::arg_match0(component, c("ar", "ma"))
  pattern = paste0("^", component, "([0-9]+)_$")
  coefficient_names = grep(pattern, names(values), value = TRUE)
  if (length(coefficient_names) == 0)
    return(rep(FALSE, nrow(values)))
  orders = as.integer(sub(pattern, "\\1", coefficient_names))
  coefficients = as.matrix(values[, coefficient_names[order(orders)], drop = FALSE])
  if (ncol(coefficients) == 1)
    return(!is.finite(coefficients[, 1]) | abs(coefficients[, 1]) >= 1)

  apply(coefficients, 1, function(x) {
    polynomial = c(1, if (component == "ar") -x else x)
    any(!is.finite(x)) || any(Mod(polyroot(polynomial)) <= 1)
  })
}


#' Warn when a posterior AR/MA root smoke test finds violations
#'
#' @param fit An `mcpfit` object with posterior draws.
#' @param ndraws,nrows Maximum numbers of draws and observed rows to check.
#' @param diagnostics A resolved diagnostics configuration.
#' @return `NULL`, invisibly.
#' @noRd
warn_arma_fit = function(fit, ndraws = 500, nrows = 100, diagnostics = list()) {
  diagnostics = resolve_diagnostics(diagnostics)
  model_tables = get_fit_model_tables(fit)
  model_predictors = model_tables$predictors
  model_group_effects = model_tables$group_effects
  components = intersect(c("ar", "ma"), unique(model_predictors$dpar))
  enabled = vapply(components, function(component) {
    !is.null(diagnostics[[component]])
  }, logical(1))
  if (!any(enabled))
    return(invisible(NULL))

  rows = unique(round(seq(1, nrow(fit$data), length.out = min(nrows, nrow(fit$data)))))
  newdata = fit$data[rows, , drop = FALSE]
  newdata$.data_row = seq_len(nrow(newdata))
  group_info = unpack_group_effects(fit, pars = TRUE)
  draws = as.matrix(.subset2(fit, "mcmc_post"))
  # Spread the check over all retained post-warmup draws and chains.
  keep = unique(round(seq(1, nrow(draws), length.out = min(ndraws, nrow(draws)))))
  smoke_fit = fit
  smoke_fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws[keep, , drop = FALSE]))
  draws = mcp_draws(
    smoke_fit, population = TRUE, group = length(group_info$cols) > 0
  )
  predictor_data = add_rhs_predictors(newdata, fit)
  if (length(group_info$cols) > 0) {
    draws_predictors = dplyr::left_join(
      predictor_data, draws, by = unique(group_info$cols),
      relationship = "many-to-many"
    )
  } else {
    draws_predictors = tidyr::expand_grid(draws, predictor_data)
  }

  values = evaluate_model_dpars(
    fit, as.list(draws_predictors),
    paste0(".pred_", get_predictor_design_names(model_predictors, model_group_effects))
  )
  probabilities = vapply(components, function(component) {
    violations = arma_root_violations(values, component)
    max(tapply(violations, draws_predictors$.data_row, mean))
  }, numeric(1))
  bad = vapply(names(probabilities), function(component) {
    threshold = diagnostics[[component]]
    !is.null(threshold) && probabilities[[component]] > threshold
  }, logical(1))
  if (!any(bad))
    return(invisible(NULL))

  details = paste0(
    toupper(names(probabilities)[bad]), ": ",
    round(100 * probabilities[bad]), "%"
  )
  warning(
    "Posterior AR/MA root smoke test found violations at observed predictor values ",
    "(maximum checked-draw violation rate: ",
    paste(details, collapse = "; "), "). ",
    "For time-varying coefficients this is a local check, not proof of global ",
    "stationarity or invertibility.",
    call. = FALSE
  )
  invisible(NULL)
}


# Warn when fresh-series AR/MA coefficients violate root conditions
# - values: Evaluated model parameters, including ar1_, ma1_, etc.
warn_arma_simulation = function(values) {
  bad = vapply(c("ar", "ma"), function(component) {
    any(arma_root_violations(values, component))
  }, logical(1))
  if (!any(bad))
    return(invisible(NULL))

  warning(
    "Generating a fresh series with locally non-",
    if (all(bad)) "stationary AR and non-invertible MA" else if (bad["ar"]) {
      "stationary AR"
    } else {
      "invertible MA"
    },
    " coefficients.",
    call. = FALSE
  )
  invisible(NULL)
}


#' Get JAGS code for GARMA residual recursion
#'
#' @param ar_order,ma_order Positive integer or `NA` when absent.
#' @param x_name Character. Name of some vector that has the length of the dataset.
#' @param series Whether generated lags must remain within `series_id_`.
#' @return Character JAGS code
#' @noRd
get_garma_jagscode = function(ar_order, ma_order, x_name, series = FALSE) {
  ar_order = ifelse(is.na(ar_order), 0, ar_order)
  ma_order = ifelse(is.na(ma_order), 0, ma_order)
  checkmate::assert_int(ar_order, lower = 0)
  checkmate::assert_int(ma_order, lower = 0)
  checkmate::assert_string(x_name)
  checkmate::assert_flag(series)
  max_order = max(ar_order, ma_order)
  if (max_order == 0)
    stop_github("get_garma_jagscode() requires a positive AR or MA order.")

  get_terms = function(row, available_ar, available_ma) {
    terms = character()
    term_lags = integer()
    if (available_ar > 0) {
      lags = seq_len(available_ar)
      terms = c(terms, paste0("ar", lags, "_[", row, "] * resid_abs_[", row, " - ", lags, "]"))
      term_lags = c(term_lags, lags)
    }
    if (available_ma > 0) {
      lags = seq_len(available_ma)
      terms = c(terms, paste0("ma", lags, "_[", row, "] * resid_ma_[", row, " - ", lags, "]"))
      term_lags = c(term_lags, lags)
    }
    if (series)
      terms = paste0(
        "equals(series_id_[", row, "], series_id_[", row, " - ",
        term_lags, "]) * ", terms
      )
    paste0(terms, collapse = " +\n              ")
  }

  jagscode = "
  # Apply GARMA recursion to link-scale residuals
  resid_garma_[1] = 0"

  if (max_order >= 2) {
    for (i in 2:max_order) {
      jagscode = paste0(
        jagscode, "\n  resid_garma_[", i, "] = ",
        get_terms(i, min(ar_order, i - 1), min(ma_order, i - 1))
      )
    }
  }

  paste0(
    jagscode,
    "\n  for (i_ in ", max_order + 1, ":length(", x_name, ")) {",
    "\n    resid_garma_[i_] = ", get_terms("i_", ar_order, ma_order),
    "\n  }"
  )
}



# Build the observation-threshold formula used by GARMA terms
#
# A threshold supplied with an AR or MA term applies until the next such
# term. The first supplied threshold also applies to earlier observations so
# they can safely be used as lags.
get_garma_threshold_jagscode = function(segments, predictors, par_x) {
  threshold_table = predictors %>%
    dplyr::filter(.data$dpar %in% c("ar", "ma"), !is.na(.data$threshold)) %>%
    dplyr::distinct(.data$segment, .data$threshold) %>%
    dplyr::arrange(.data$segment)

  if (nrow(threshold_table) == 0)
    return("")
  if (anyDuplicated(threshold_table$segment))
    stop_github("Found multiple GARMA thresholds in one segment.")

  threshold_code = stats::setNames(segments$cp_code_form, segments$segment)
  threshold_parts = character(nrow(threshold_table))
  for (i in seq_len(nrow(threshold_table))) {
    lower = if (i == 1) "" else paste0("(", par_x, "[i_] >= ", threshold_code[[as.character(threshold_table$segment[i])]], ") * ")
    upper = if (i == nrow(threshold_table)) "" else paste0("(", par_x, "[i_] < ", threshold_code[[as.character(threshold_table$segment[i + 1])]], ") * ")
    threshold_value = sprintf("%.15g", threshold_table$threshold[i])
    threshold_parts[i] = paste0("  ", lower, upper, threshold_value)
  }

  paste0(
    "\n\n# GARMA observation threshold\n",
    "garma_threshold_[i_] =\n",
    paste0(threshold_parts, collapse = " +\n")
  )
}


# Evaluate or generate a GARMA response series
simulate_garma = function(base_link_mu, ar_list, ma_list, threshold, family,
                          dpars, data = list(), y = NULL, series_id = NULL) {
  assert_arma_thresholds(family, threshold, data)

  if (is.null(series_id))
    series_id = rep(1, length(base_link_mu))
  if (length(series_id) != length(base_link_mu) || anyNA(series_id))
    stop_github("series_id must have one non-missing value per observation.")
  series_runs = rle(as.character(series_id))
  if (anyDuplicated(series_runs$values))
    stop_github("Rows belonging to each series_id must be contiguous.")

  generate_series = is.null(y)

  ar_order = length(ar_list)
  ma_order = length(ma_list)
  n = length(base_link_mu)

  # One recurrence for supplied and generated histories. At each time position,
  # independent series/draws can be evaluated together, including ragged series.
  if (generate_series)
    y = rep(NA_real_, n)
  link_y = family$linkfun(family$garma$observed_r(y, data, threshold))
  resid_abs = numeric(n)
  resid_ma = numeric(n)
  resid_garma = numeric(n)
  link_mu = numeric(n)
  positions = sequence(series_runs$lengths)
  rows_by_position = split(seq_len(n), factor(positions, levels = seq_len(max(positions))))

  for (rows in rows_by_position) {
    position = positions[rows[1]]
    for (lag in seq_len(min(ar_order, position - 1))) {
      resid_garma[rows] = resid_garma[rows] +
        ar_list[[paste0("ar", lag, "_")]][rows] * resid_abs[rows - lag]
    }
    for (lag in seq_len(min(ma_order, position - 1))) {
      resid_garma[rows] = resid_garma[rows] +
        ma_list[[paste0("ma", lag, "_")]][rows] * resid_ma[rows - lag]
    }
    link_mu[rows] = base_link_mu[rows] + resid_garma[rows]

    # Only absent history values are generated here. Conditional predictions
    # are separate response draws and never replace the supplied history.
    missing_rows = rows[is.na(y[rows])]
    if (length(missing_rows) > 0) {
      row_dpars = lapply(dpars, function(x) if (length(x) == 1) x else x[missing_rows])
      row_dpars$mu = family$linkinv(link_mu[missing_rows])
      row_data = lapply(data, function(x) if (length(x) == 1) x else x[missing_rows])
      y[missing_rows] = family$r$rng(length(missing_rows), row_dpars, row_data)
      link_y[missing_rows] = family$linkfun(family$garma$observed_r(
        y[missing_rows], row_data, threshold[missing_rows]
      ))
    }
    resid_abs[rows] = link_y[rows] - base_link_mu[rows]
    resid_ma[rows] = link_y[rows] - link_mu[rows]
  }

  list(
    y = y,
    mu = family$linkinv(link_mu),
    link_mu = link_mu,
    resid_garma = resid_garma,
    resid_abs = resid_abs,
    resid_ma = resid_ma
  )
}
