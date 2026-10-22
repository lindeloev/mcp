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


# Warn when a segment does not include the ar()/ma() of the preceding segment, which then ends
# there. ar() continued into later segments in mcp v0.3.4. Including ar(0) silences the warning.
warn_arma_not_included = function(components) {
  included = get_included_arma(components)
  next_rows = included[included$order > 0 & included$segment < length(components), ]
  next_rows$segment = next_rows$segment + 1
  next_rows = dplyr::anti_join(next_rows, included, by = c("dpar", "segment"))
  for (i in seq_len(nrow(next_rows))) {
    dpar = next_rows$dpar[i]
    segment = next_rows$segment[i]
    warning(
      "Segment ", segment, " does not include ", dpar, "(), so it has no ", dpar, "(", next_rows$order[i], ") component. If this is intentional, you can ignore this warning or add `", dpar, "(0)` to the segment to silence it. ",
      "To join the ", dpar, "()-component from segment ", segment - 1, ", include `", dpar, "(", next_rows$order[i], ", 0)` in segment ", segment, ".",
      if (dpar == "ar") " Before mcp v0.4.0, ar() continued into later segments.",
      call. = FALSE
    )
  }
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


# Validate the AR/MA configuration of a model before sampling
check_arma_model = function(family, predictors, segments, data, par_x, series = NULL) {
  if (is.null(family$garma))
    stop(
      "family = ", family$family, "(link = \"", family$link,
      "\") does not define the GARMA behavior required by ar() or ma()."
    )

  response_data = get_family_response_data(family, segments, data)
  assert_arma_thresholds(family, predictors$threshold[predictors$dpar %in% c("ar", "ma")], response_data)

  x_by_series = split(data[[par_x]], if (is.null(series)) 1 else data[[series]])
  x_unordered = any(vapply(
    x_by_series,
    function(x) is.unsorted(x) && is.unsorted(rev(x)),
    logical(1)
  ))
  if (x_unordered)
    message("'", par_x, "' is unordered. Please note that ar() and ma() apply in data-frame row order, not the values of '", par_x, "'.")

  invisible(TRUE)
}


# Validate GARMA observation thresholds against the support of the response
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


# Whether an evaluation runs the GARMA recurrence. AR/MA terms change the expected
# response and mu, but not other dpars such as sigma or ar1.
uses_garma = function(fit, arma, type, dpar = NULL) {
  arma && is_arma(fit) && (type != "fitted" || is.null(dpar) || dpar %in% c("epred", "mu"))
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


# The response with missing values replaced by their retained JAGS imputations: a rows x draws matrix.
# - y: The response of the fitted data, one value per row.
# - draw_ids: `.draw` indices of the draws, i.e., rows of the stacked posterior chains.
get_imputed_response = function(fit, y, draw_ids) {
  rows = fit$.internal$imputed_response_rows
  nodes = paste0(mcp_columns(fit)$response, "[", rows, "]")
  y = matrix(y, length(y), length(draw_ids))
  y[rows, ] = t(as.matrix(fit$.internal$imputed_response)[draw_ids, nodes, drop = FALSE])
  y
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


# The AR or MA lags of evaluated linear predictors (see `evaluate_formula()`), ordered by lag,
# e.g., `list(ar1 = , ar2 = )`. An empty list without such terms.
get_lags = function(links, component) {
  keys = grep(paste0("^", component, "[0-9]+$"), names(links), value = TRUE)
  links[keys[order(as.integer(substring(keys, 3)))]]
}


# Check AR stationarity or MA invertibility for each row and draw. Returns a rows x draws matrix.
# - links: Evaluated linear predictors from `evaluate_formula()`, including "ar1", "ma1", etc.
# - component: Either "ar" or "ma".
arma_root_violations = function(links, component) {
  component = rlang::arg_match0(component, c("ar", "ma"))
  lags = get_lags(links, component)
  if (length(lags) == 0)
    return(FALSE)
  coefficients = do.call(cbind, lapply(lags, as.vector))  # One column per lag
  violations = if (ncol(coefficients) == 1) {
    !is.finite(coefficients[, 1]) | abs(coefficients[, 1]) >= 1
  } else {
    apply(coefficients, 1, function(x) {
      polynomial = c(1, if (component == "ar") -x else x)
      any(!is.finite(x)) || any(Mod(polyroot(polynomial)) <= 1)
    })
  }
  matrix(violations, nrow(lags[[1]]))
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
  components = intersect(c("ar", "ma"), unique(get_fit_model_tables(fit)$predictors$dpar))
  enabled = vapply(components, function(component) {
    !is.null(diagnostics[[component]])
  }, logical(1))
  if (!any(enabled))
    return(invisible(NULL))

  rows = unique(round(seq(1, nrow(fit$data), length.out = min(nrows, nrow(fit$data)))))
  newdata = fit$data[rows, , drop = FALSE]
  draws = as.matrix(.subset2(fit, "mcmc_post"))
  # Spread the check over all retained post-warmup draws and chains.
  keep = unique(round(seq(1, nrow(draws), length.out = min(ndraws, nrow(draws)))))
  smoke_fit = fit
  smoke_fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws[keep, , drop = FALSE]))
  draws = mcp_draws(smoke_fit, population = TRUE, group = TRUE, pivot = FALSE)
  effects = unpack_group_effects(fit, pars = TRUE)$effects
  pars = draws_to_pars(draws, newdata, split(effects$name, effects$group_col))
  links = evaluate_formula(fit, get_model_data(fit, newdata), pars, nrow(draws))

  # Maximum violation rate across rows
  probabilities = vapply(components, function(component) {
    max(rowMeans(arma_root_violations(links, component)))
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
# - links: Evaluated linear predictors from `evaluate_formula()`, including "ar1", "ma1", etc.
warn_arma_simulation = function(links) {
  bad = vapply(c("ar", "ma"), function(component) {
    any(arma_root_violations(links, component))
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



# Build the observation-threshold node used by GARMA terms. One threshold applies to all segments.
get_garma_threshold_jagscode = function(predictors) {
  threshold = unique(stats::na.omit(predictors$threshold))
  if (length(threshold) == 0)
    return("")

  paste0("\n\n# GARMA observation threshold\ngarma_threshold_[i_] = ", sprintf("%.15g", threshold))
}


#' Evaluate or generate GARMA response series
#'
#' Rows are time points, contiguous within each series, and columns are independent draws.
#' At each time position, all series and draws are evaluated together.
#'
#' @param link_mu The linear predictor of mu without AR/MA terms (rows x draws).
#' @param ar,ma Lists of AR and MA coefficients (rows x draws), one per lag.
#' @param threshold Observation threshold that keeps the link of the response finite.
#' @param dpars Distributional parameters on the response scale (rows x draws) to generate responses.
#' @param data Response auxiliaries, e.g., `trials`, one value per row.
#' @param y The observed response, one value per row or rows x draws. `NULL` generates fresh series.
#' @param series The series of each row, or `NULL` for a single series.
#' @return A list with the response `y` (rows x draws), `link_mu` including the AR/MA terms, and these AR/MA terms (`resid_garma`).
#' @noRd
simulate_garma = function(link_mu, ar, ma, threshold, family, dpars, data = list(), y = NULL, series = NULL) {
  assert_arma_thresholds(family, threshold, data)
  if (is.null(series))
    series = rep(1, nrow(link_mu))
  runs = rle(match(series, unique(series)))
  if (anyDuplicated(runs$values))
    stop_github("Rows belonging to each series must be contiguous.")

  generate = is.null(y)
  y = matrix(if (generate) NA_real_ else y, nrow(link_mu), ncol(link_mu))
  link_y = matrix(family$linkfun(family$garma$observed_r(y, data, threshold)), nrow(y))
  resid_abs = resid_ma = resid_garma = matrix(0, nrow(y), ncol(y))

  position = sequence(runs$lengths)
  for (rows in split(seq_along(position), position)) {
    p = position[rows[1]]
    for (lag in seq_len(min(length(ar), p - 1)))
      resid_garma[rows, ] = resid_garma[rows, ] + ar[[lag]][rows, ] * resid_abs[rows - lag, ]
    for (lag in seq_len(min(length(ma), p - 1)))
      resid_garma[rows, ] = resid_garma[rows, ] + ma[[lag]][rows, ] * resid_ma[rows - lag, ]

    # A fresh series draws each response given the generated history
    if (generate) {
      row_dpars = lapply(dpars, function(value) value[rows, , drop = FALSE])
      row_dpars$mu = family$linkinv(link_mu[rows, , drop = FALSE] + resid_garma[rows, , drop = FALSE])
      row_data = lapply(data, function(value) value[rows])
      y[rows, ] = family$r$rng(length(rows) * ncol(y), row_dpars, row_data)
      link_y[rows, ] = family$linkfun(family$garma$observed_r(y[rows, ], row_data, threshold))
    }
    resid_abs[rows, ] = link_y[rows, ] - link_mu[rows, ]
    resid_ma[rows, ] = resid_abs[rows, ] - resid_garma[rows, ]
  }

  list(y = y, link_mu = link_mu + resid_garma, resid_garma = resid_garma)
}
