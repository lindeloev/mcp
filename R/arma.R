# ABOUT: AR/MA (autoregressive/moving-average, "GARMA") specific code:
# formula parsing, JAGS code generation, R-side simulation/evaluation, and
# root-condition (stationarity/invertibility) checks.
# ------------------------------------------------------------

#' Unpack arma order and formula
#'
#' @aliases unpack_arma
#' @keywords internal
#' @noRd
#' @param form_str_in A character such as `"ar(number)"`, `"ma(number)"`, or
#'   either component with a second formula argument.
#' @return A list with `$order`, `$form_str` (e.g., `"ar(formula)"`), and
#'   `$boundary`. The component formula is 1 if no formula is given.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_arma = function(form_str_in) {
  if (length(form_str_in) == 0) {
    return(list(
      order = NA,
      form_str = NA
    ))
  } else if (length(form_str_in) > 1) {
    stop("Only one of these allowed per segment: ", form_str_in)
  }

  component_call = str2lang(form_str_in)
  component = as.character(component_call[[1]])
  component = rlang::arg_match0(component, c("ar", "ma"))

  component_args = as.list(component_call)[-1]
  component_arg_names = names(component_args)
  if (is.null(component_arg_names))
    component_arg_names = rep("", length(component_args))

  if (length(component_args) == 0 || component_arg_names[1] %notin% c("", "order"))
    stop("The first argument to ", component, "() must be its order.")

  boundary_index = which(component_arg_names == "boundary")
  if (length(boundary_index) > 1)
    stop("Only one `boundary` value is allowed in ", component, "().")

  formula_index = setdiff(seq_along(component_args), c(1, boundary_index))
  if (length(formula_index) > 1 || any(component_arg_names[formula_index] %notin% c("", "formula")))
    stop(component, "() accepts only `order`, an optional formula, and `boundary`.")

  # GET ORDER
  order_str = paste(deparse(component_args[[1]], width.cutoff = 500), collapse = "")
  order = suppressWarnings(as.numeric(order_str))

  # Check the order
  if (is.na(order))
    stop("Wrong specification of order in '", form_str_in, "'. Must be ", component, "(order) or ", component, "(order, formula) where order is a positive integer.")
  checkmate::assert_int(order, lower = 1, .var.name = form_str_in)

  # GET FORMULA AND BOUNDARY
  if (length(formula_index) == 1) {
    formula_str = paste(deparse(component_args[[formula_index]], width.cutoff = 500), collapse = "")
    form_str = paste0(component, "(", formula_str, ")")
  } else {
    # If there is no formula, use an intercept-only component formula.
    form_str = paste0(component, "(1)")
  }

  if (length(boundary_index) == 1) {
    boundary_str = paste(deparse(component_args[[boundary_index]], width.cutoff = 500), collapse = "")
    boundary = suppressWarnings(as.numeric(boundary_str))
    if (length(boundary) != 1 || is.na(boundary) || !is.finite(boundary) || boundary <= 0 || boundary >= 1)
      stop("`boundary` in ", component, "() must be one number between 0 and 1.")
  } else {
    # Defer the default until AR and MA have been parsed together for a segment.
    boundary = NA_real_
  }

  # Return
  list(
    order = order,
    form_str = form_str,
    boundary = boundary
  )
}


# Returns the requested AR/MA order, or NA if the term is absent
get_arma_order = function(predictors, term) {
  term = rlang::arg_match0(term, c("ar", "ma"))
  orders = predictors$order[predictors$dpar == term]
  if (length(orders) == 0) NA else max(orders, na.rm = TRUE)
}


#' Check if this is an AR/MA model
#'
#' @aliases is_arma
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return TRUE or FALSE
is_arma = function(fit) {
  length(fit$pars$arma) > 0
}


#' Warn about model-checking limitations for AR/MA models
#'
#' @keywords internal
#' @noRd
#' @param fit An mcpfit object.
#' @param arma Whether AR and MA effects are included in the evaluation.
#' @param check One of `"ppc"` or `"information_criterion"`.
#' @return `NULL`, invisibly. Called for the warning side-effect.
warn_arma_check = function(fit, arma, check) {
  if (!arma || !is_arma(fit))
    return(invisible(NULL))

  warning(
    switch(
      check,
      ppc = paste(
        "For AR/MA models, mcp conditions predictions on the observed response history.",
        "These are one-step-ahead conditional predictions, not jointly replicated time series.",
        "Serial summaries such as ACF and run lengths may therefore be misleading.",
        "Joint-series posterior predictive checks are not yet supported."
      ),
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


#' Check AR stationarity or MA invertibility row by row
#'
#' @keywords internal
#' @noRd
#' @param values Evaluated model parameters, including `ar1_`, `ma1_`, etc.
#' @param component Either `"ar"` or `"ma"`.
#' @return A logical vector.
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
#' @keywords internal
#' @noRd
#' @param fit An `mcpfit` object with posterior draws.
#' @param ndraws,nrows Maximum numbers of draws and observed rows to check.
#' @param diagnostics A resolved diagnostics configuration.
#' @return `NULL`, invisibly.
warn_arma_fit = function(fit, ndraws = 500, nrows = 100, diagnostics = list()) {
  diagnostics = resolve_diagnostics(diagnostics)
  model_predictors = get_fit_model_tables(fit)$predictors
  components = intersect(c("ar", "ma"), unique(model_predictors$dpar))
  enabled = vapply(components, function(component) {
    !is.null(diagnostics[[component]])
  }, logical(1))
  if (!any(enabled))
    return(invisible(NULL))

  rows = unique(round(seq(1, nrow(fit$data), length.out = min(nrows, nrow(fit$data)))))
  newdata = fit$data[rows, , drop = FALSE]
  newdata$data_row = seq_len(nrow(newdata))
  group_info = unpack_varying(fit, pars = TRUE)
  draws = as.matrix(.subset2(fit, "mcmc_post"))
  # Spread the check over all retained post-warmup draws and chains.
  keep = unique(round(seq(1, nrow(draws), length.out = min(ndraws, nrow(draws)))))
  smoke_fit = fit
  smoke_fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws[keep, , drop = FALSE]))
  draws = mcp_draws(
    smoke_fit, population = TRUE, varying = length(group_info$cols) > 0
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
    paste0(".pred_", model_predictors$code_name)
  )
  components = intersect(c("ar", "ma"), unique(model_predictors$dpar))
  probabilities = vapply(components, function(component) {
    violations = arma_root_violations(values, component)
    max(tapply(violations, draws_predictors$data_row, mean))
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
    "stationarity or invertibility. ",
    "See `vignette(\"arma\")`.",
    call. = FALSE
  )
  invisible(NULL)
}


#' Warn when fresh-series AR/MA coefficients violate root conditions
#'
#' @keywords internal
#' @noRd
#' @inheritParams arma_root_violations
#' @return `NULL`, invisibly.
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
    " coefficients. See `vignette(\"arma\")`.",
    call. = FALSE
  )
  invisible(NULL)
}


#' Get JAGS code for GARMA residual recursion
#'
#' @aliases get_arma_jagscode get_ar_jagscode
#' @keywords internal
#' @noRd
#' @param ar_order,ma_order Positive integer or `NA` when absent.
#' @param x_name Character. Name of some vector that has the length of the dataset.
#' @return Character JAGS code
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_arma_jagscode = function(ar_order, ma_order, x_name) {
  ar_order = ifelse(is.na(ar_order), 0, ar_order)
  ma_order = ifelse(is.na(ma_order), 0, ma_order)
  checkmate::assert_int(ar_order, lower = 0)
  checkmate::assert_int(ma_order, lower = 0)
  checkmate::assert_string(x_name)
  max_order = max(ar_order, ma_order)
  if (max_order == 0)
    stop_github("get_arma_jagscode() requires a positive AR or MA order.")

  get_terms = function(row, available_ar, available_ma) {
    terms = character()
    if (available_ar > 0) {
      lags = seq_len(available_ar)
      terms = c(terms, paste0("ar", lags, "_[", row, "] * resid_abs_[", row, " - ", lags, "]"))
    }
    if (available_ma > 0) {
      lags = seq_len(available_ma)
      terms = c(terms, paste0("ma", lags, "_[", row, "] * resid_ma_[", row, " - ", lags, "]"))
    }
    paste0(terms, collapse = " +\n              ")
  }

  jagscode = "
  # Apply GARMA recursion to link-scale residuals
  resid_arma_[1] = 0"

  if (max_order >= 2) {
    for (i in 2:max_order) {
      jagscode = paste0(
        jagscode, "\n  resid_arma_[", i, "] = ",
        get_terms(i, min(ar_order, i - 1), min(ma_order, i - 1))
      )
    }
  }

  paste0(
    jagscode,
    "\n  for (i_ in ", max_order + 1, ":length(", x_name, ")) {",
    "\n    resid_arma_[i_] = ", get_terms("i_", ar_order, ma_order),
    "\n  }"
  )
}


# Backwards-compatible internal wrapper for pure AR code generation.
get_ar_jagscode = function(ar_order, x_name) {
  get_arma_jagscode(ar_order, NA, x_name)
}


#' Build the observation-boundary formula used by GARMA terms
#'
#' A boundary supplied with an AR or MA term remains active until the next such
#' term. The first supplied boundary also applies to earlier observations so
#' they can safely be used as lags.
#'
#' @keywords internal
#' @noRd
get_garma_boundary_jagscode = function(segments, predictors, par_x) {
  boundary_table = predictors %>%
    dplyr::filter(.data$dpar %in% c("ar", "ma"), !is.na(.data$boundary)) %>%
    dplyr::distinct(.data$segment, .data$boundary) %>%
    dplyr::arrange(.data$segment)

  if (nrow(boundary_table) == 0)
    return("")
  if (anyDuplicated(boundary_table$segment))
    stop_github("Found multiple GARMA boundaries in one segment.")

  boundary_code = stats::setNames(segments$cp_code_form, segments$segment)
  boundary_parts = character(nrow(boundary_table))
  for (i in seq_len(nrow(boundary_table))) {
    lower = if (i == 1) "" else paste0("(", par_x, "[i_] >= ", boundary_code[[as.character(boundary_table$segment[i])]], ") * ")
    upper = if (i == nrow(boundary_table)) "" else paste0("(", par_x, "[i_] < ", boundary_code[[as.character(boundary_table$segment[i + 1])]], ") * ")
    boundary_value = sprintf("%.15g", boundary_table$boundary[i])
    boundary_parts[i] = paste0("  ", lower, upper, boundary_value)
  }

  paste0(
    "\n\n# GARMA observation boundary\n",
    "garma_boundary_[i_] =\n",
    paste0(boundary_parts, collapse = " +\n")
  )
}


#' Evaluate or generate a GARMA response series
#'
#' @keywords internal
#' @noRd
simulate_garma = function(base_link_mu, ar_list, ma_list, boundary, family,
                          dpars, data = list(), y = NULL, series_id = NULL) {
  if (is.null(series_id))
    series_id = rep(1, length(base_link_mu))
  if (length(series_id) != length(base_link_mu) || anyNA(series_id))
    stop_github("series_id must have one non-missing value per observation.")

  generate_series = is.null(y)
  if (generate_series && !is.null(family$garma$generate_message))
    message(family$garma$generate_message)

  ar_order = length(ar_list)
  ma_order = length(ma_list)
  resid_abs = numeric(length(base_link_mu))
  resid_ma = numeric(length(base_link_mu))
  resid_arma = numeric(length(base_link_mu))
  link_mu = numeric(length(base_link_mu))
  mu = numeric(length(base_link_mu))
  if (generate_series)
    y = numeric(length(base_link_mu))

  for (rows in split(seq_along(base_link_mu), series_id)) {
    for (position in seq_along(rows)) {
      row = rows[position]
      for (lag in seq_len(min(ar_order, position - 1))) {
        resid_arma[row] = resid_arma[row] +
          ar_list[[paste0("ar", lag, "_")]][row] * resid_abs[rows[position - lag]]
      }
      for (lag in seq_len(min(ma_order, position - 1))) {
        resid_arma[row] = resid_arma[row] +
          ma_list[[paste0("ma", lag, "_")]][row] * resid_ma[rows[position - lag]]
      }

      link_mu[row] = base_link_mu[row] + resid_arma[row]
      mu[row] = family$linkinv(link_mu[row])
      generate_observation = generate_series || is.na(y[row])
      if (generate_observation) {
        row_dpars = lapply(dpars, function(x) if (length(x) == 1) x else x[row])
        row_dpars$mu = mu[row]
        row_data = lapply(data, function(x) if (length(x) == 1) x else x[row])
        y[row] = family$r$rng(1, row_dpars, row_data)
      }

      row_data = lapply(data, function(x) if (length(x) == 1) x else x[row])
      garma_y = get_garma_observed(y[row], family, boundary[row], row_data)
      garma_link_y = family$linkfun(garma_y)
      resid_abs[row] = garma_link_y - base_link_mu[row]
      resid_ma[row] = garma_link_y - link_mu[row]
    }
  }

  list(
    y = y,
    mu = mu,
    link_mu = link_mu,
    resid_arma = resid_arma,
    resid_abs = resid_abs,
    resid_ma = resid_ma
  )
}


#' Simulate/evaluate autoregressive residuals
#'
#' Developer note: some of the eval(parse(text = ...)) here could probably be replaced with inner products (%*%).
#' @aliases simulate_ar
#' @keywords internal
#' @noRd
#' @param sigma_ Numeric vector of innovations
#' @param ar_list List with numerical vectors, list(ar1_ = c(...), ar2_ = c(...))
#' @param resid_abs NULL or numerical vector of observed minus fitted residuals
#'   on the GARMA link scale.
#' @param series_id Optional vector identifying independent series. Lagged
#'   residuals never cross from one series to another.
#' @return List with
#'   * `resid_ar`: the ARMA part of the residuals
#'   * `resid_sigma`: the innovations.
#'
#'   Note that `resid_abs = resid_ar + resid_sigma`.
#' @seealso get_ar_jagscode
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
simulate_ar = function(sigma_, ar_list, resid_abs = NULL, series_id = NULL) {
  # Check inputs
  checkmate::assert_numeric(sigma_, any.missing = FALSE)
  checkmate::assert_list(ar_list)
  checkmate::assert_numeric(resid_abs, any.missing = FALSE, null.ok = TRUE)
  if (length(grep("^ar[0-9]+_$", names(ar_list))) != length(ar_list))
    stop_github("Not all names(ar_list) are arx_.")
  if (is.null(series_id))
    series_id = rep(1, length(sigma_))
  if (length(series_id) != length(sigma_) || anyNA(series_id))
    stop_github("series_id must have one non-missing value per residual.")
  if (length(resid_abs) > 0 && length(resid_abs) != length(sigma_))
    stop_github("resid_abs and sigma_ must have the same length.")

  ar_order = length(ar_list)

  # resid_ is the observed residual from y_
  # resid_ is split into the innovation and AR() part. So resid_ = resid_ar + resid_sigma
  resid_sigma = stats::rnorm(length(sigma_), 0, sigma_)
  generate_residuals = length(resid_abs) == 0
  if (generate_residuals) {
    message("Generating residuals for AR(N) model since the response column/argument was not provided.")
    resid_abs = numeric(length(sigma_))
  }

  resid_ar = numeric(length(resid_abs))
  for (rows in split(seq_along(resid_abs), series_id)) {
    for (position in seq_along(rows)) {
      row = rows[position]
      available_lags = seq_len(min(ar_order, position - 1))
      for (lag in available_lags) {
        resid_ar[row] = resid_ar[row] +
          ar_list[[paste0("ar", lag, "_")]][row] * resid_abs[rows[position - lag]]
      }
      if (generate_residuals)
        resid_abs[row] = resid_sigma[row] + resid_ar[row]
    }
  }

  # Return
  list(
    resid_ar = resid_ar,
    resid_sigma = resid_sigma
  )
}
