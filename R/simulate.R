# ABOUT: These functions take data/samples and run them through
# formula code to make model predictions
# ------------------------------------------------------------

#' Make the levels of categorical predictors match the original data
#' @aliases relevel_newdata
#' @keywords internal
#' @noRd
#' @inheritParams add_rhs_predictors
#' @return `newdata` with all categorical columns as factors that have identical levels to the original data.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
relevel_newdata = function(newdata, fit) {
  # Check that the necessary data is available
  rhs_vars = get_rhs_vars(fit$model)
  assert_data_cols(newdata, cols = rhs_vars, fail_funcs = c(is.na, is.nan, is.infinite))

  # Make sure to carry over the exact level-structure of the original model
  for (col_name in rhs_vars) {
    org_col = fit$data[, col_name]
    new_col = newdata[, col_name]
    if (is.character(org_col) | is.factor(org_col)) {
      new_col = factor(new_col, levels = levels(factor(org_col)))

      # Helpful error
      if (any(is.na(new_col))) {
        new_levels = newdata[is.na(new_col), col_name] %>% as.character() %>% unique()
        stop("Got novel values (", and_collapse(new_levels), ") for column ", col_name, ". Only values used during fitting are allowed.")
      }
    }

    newdata[, col_name] = new_col
  }

  newdata
}


#' Add predictors to `newdata`
#' @aliases add_rhs_predictors
#' @keywords internal
#' @noRd
#' @param newdata A data.frame that contains the predictor columns from
#'   `fit$data` and `fit$model`.
#' @param fit An `mcpfit` object.
#' @return `newdata` with additional (dummy-coded) columns that make up the design matrix required for the model.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
add_rhs_predictors = function(newdata, fit) {
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_data_frame(newdata)

  # Make categorical predictors match the original data
  newdata = newdata %>%
    as.data.frame() %>%
    relevel_newdata(fit)

  # Get predictor matrix
  predictors = get_fit_model_tables(fit)$predictors
  new_predictors = get_predictors(
    fit$model, newdata, fit$family, fit$pars$x, check_rank = FALSE
  )
  predictor_matrix = get_predictor_matrix(new_predictors)

  # Check that the variable structure matches
  if (nrow(predictors) != nrow(new_predictors) ||
      any(predictors$code_name != new_predictors$code_name))
    stop_github("The new predictors table does not match the fitted model.")

  # All permutations of rows in newdata and parameters
  as.data.frame(predictor_matrix) %>%
    magrittr::set_colnames(paste0(".pred_", colnames(predictor_matrix))) %>%
    dplyr::bind_cols(newdata)
}


#' Returns parameters needed for simulation
#' @aliases get_sim_pars
#' @keywords internal
#' @return Character vector
get_sim_pars = function(predictors, pars) {
  c(
    pars$cp,  # cp_1, cp_2, etc.
    predictors$code_name,  # mu, sigma, ar, etc.
    pars$varying
  )
}



#' Evaluate model and return a list of dpars
#'
#' This is currently hard-coded to be run from `simulate_vectorized`.
#' It serves to scope the evaluation of the model to prevent name conflicts.
#'
#' @keywords internal
#' @noRd
#' @param fit An `mcpfit` object.
#' @param args args from `simulate_vectorized`
#' @return `data.frame` with one column per dpar
#' @noRd
evaluate_model_dpars = function(fit, args, pred_pars) {
  # Generate more predictors
  pred_args = args[names(args) %in% pred_pars]
  rhs_matrix_ = do.call(cbind, pred_args)
  rhs_matrix_ = rhs_matrix_[, match(pred_pars, colnames(rhs_matrix_)), drop = FALSE]  # Same order as predictors$code_name no matter order of args

  cp_0 = -Inf
  assign(paste0("cp_", length(fit$model)), Inf)  # e.g., cp_3 = Inf

  # Evaluate model in environment
  out = new.env()
  with(out, eval(str2expression(fit$.internal$formula_r)))

  out %>%
    as.list() %>%
    as.data.frame() %>%
    dplyr::select(-dplyr::starts_with("x_local_"), -dplyr::any_of("ar_"))
}


#' Transform link-scale predictors to distribution-scale parameters
#'
#' Mirrors the deterministic dpar transformations in generated JAGS code. Both
#' representations are kept only during R-side evaluation; neither
#' observation-level vector is monitored in JAGS.
#'
#' @keywords internal
#' @noRd
add_response_dpars = function(dpar_values, family) {
  for (dpar in family$dpar_specs$dpar) {
    spec = get_dpar_spec(family, dpar)
    link_name = paste0("link_", dpar, "_")
    response_name = paste0(dpar, "_")

    if (link_name %notin% names(dpar_values))
      stop_github("Missing link-scale values for dpar '", dpar, "'.")

    response = get_link_function(spec$link, inverse = TRUE)(dpar_values[[link_name]])
    if (!is.na(spec$lower))
      response = pmax(spec$lower, response)
    dpar_values[[response_name]] = response
  }

  dpar_values
}


#' Transform observations to a finite GARMA link domain
#'
#' @keywords internal
#' @noRd
get_garma_observed = function(y, family, boundary, data = list()) {
  if (is.null(family$garma))
    stop_github("GARMA observation transformation is unavailable for family = ", family$family, "().")
  family$garma$observed_r(y, data, boundary)
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


#' Vectorized R-side run of the full model.
#'
#' Used internally in mcp.
#' Parameters not documented here are documented in `pp_eval` (without dot prefix).
#'
#' @aliases simulate_vectorized
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param ... Parameter names (e.g., `cp_1 = c(4.3, 4.5, 6.2), Intercept_1 = c(11.2, 12.1, 10.9)`, etc.)
#'   and columns from `predictors$matrix_data` prefixed with ".pred_" (e.g.,
#'   `.pred_Intercept_1 = c(1, 1, 1)`).
#' @return Vector with same length as inputs.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
simulate_vectorized = function(fit, ..., .type = "predict", .rate = FALSE, .dpar = "epred", .arma = TRUE, .scale = "response") {
  ###########
  # ASSERTS #
  ###########
  checkmate::assert_class(fit, "mcpfit")
  if (!is.mcpfamily(fit$family))
    fit$family = mcpfamily(fit$family)
  model_tables = get_fit_model_tables(fit)
  predictors = model_tables$predictors

  # Assert that the ellipsis contains the expected argument names
  param_pars = get_sim_pars(predictors, fit$pars)
  pred_pars = paste0(".pred_", predictors$code_name)
  operation = switch(.type, fitted = "epred", loglik = "log_lik", predict = "rng")
  is_arma = any(predictors$dpar %in% c("ar", "ma"))
  aux_operations = c(operation, if (is_arma && .arma) "garma")
  aux_columns = get_family_aux_columns(fit$family, model_tables$segments, aux_operations)
  data_pars = c(fit$pars$x, stats::na.omit(unname(aux_columns)))
  expected_arg_names = c(param_pars, pred_pars, data_pars)

  args = list(...)
  missing_args = dplyr::setdiff(expected_arg_names, names(args))
  if (length(missing_args) > 0)
    stop_github("Missing the following arguments: ", and_collapse(missing_args))

  # Other args
  assert_typescale(.type, .scale)
  checkmate::assert_flag(.rate)

  .dpar = assert_dpar(.dpar, fit = fit, type = .type)
  checkmate::assert_flag(.arma)


  ##################################################
  # EVALUATE MODEL AND RETURN VIA FAMILY AND .TYPE #
  ##################################################
  dpar_values = evaluate_model_dpars(fit, args, pred_pars)
  uses_link_dpars = "link_mu_" %in% names(dpar_values)
  if (uses_link_dpars)
    dpar_values = add_response_dpars(dpar_values, fit$family)

  # Prepare for stuff needing .ydata
  has_ydata = is.null(args[[fit$pars$y]]) == FALSE
  if (has_ydata)
    dpar_values$.ydata = args[[fit$pars$y]]
  if (.type == "loglik" & has_ydata == FALSE)
    stop(".ydata must be non-NULL for .type = 'loglik'.")
  .dpar = paste0(.dpar, "_")

  if (!uses_link_dpars && .scale == "response" && .dpar %in% c("epred_", "mu_"))
    dpar_values$mu_ = fit$family$linkinv(dpar_values$mu_)

  response_data = get_family_response_data(fit$family, model_tables$segments, args)
  dpars = stats::setNames(lapply(fit$family$dpars, function(dpar) {
    dpar_values[[paste0(dpar, "_")]]
  }), fit$family$dpars)

  # Simply return for fitted dpars
  if (.dpar %notin% c("epred_", "mu_") & .type == "fitted") {
    link_dpar = paste0("link_", .dpar)
    if (uses_link_dpars && .scale == "linear" && link_dpar %in% names(dpar_values))
      return(dpar_values[[link_dpar]])
    return(dpar_values[[.dpar]])
  }

  # GARMA is defined on the link scale. The observed response is clipped only
  # where needed to keep log and logit transformations finite.
  if (is_arma && .arma == TRUE) {
    base_link_mu = if (uses_link_dpars) dpar_values$link_mu_ else fit$family$linkfun(dpar_values$mu_)
    ar_list = dplyr::select(dpar_values, dplyr::matches("^ar[0-9]+_$"))
    ma_list = dplyr::select(dpar_values, dplyr::matches("^ma[0-9]+_$"))
    boundary = dpar_values$garma_boundary_
    if (is.null(boundary))
      boundary = rep(0.1, length(base_link_mu))
    if (!has_ydata && .type != "predict")
      stop("The response is required to evaluate GARMA terms.")

    arma_result = simulate_garma(
      base_link_mu, ar_list, ma_list, boundary, fit$family,
      dpars = dpars, data = response_data,
      y = if (has_ydata) dpar_values$.ydata else NULL,
      series_id = args[[".draw"]]
    )
    if (!has_ydata)
      return(fit$family$response$observed(arma_result$y, response_data, .rate))

    dpar_values$link_mu_ = arma_result$link_mu
    dpar_values$mu_ = arma_result$mu
    dpars$mu = arma_result$mu
  }

  if (.type == "fitted" && .scale == "linear") {
    if (uses_link_dpars)
      return(dpar_values$link_mu_)
    return(dpar_values$mu_)
  }

  if (.type == "fitted")
    return(fit$family$r$epred(dpars, response_data, rate = .rate))
  if (.type == "loglik")
    return(fit$family$r$log_lik(dpar_values$.ydata, dpars, response_data))
  if (.type == "predict")
    return(fit$family$r$rng(length(dpars$mu), dpars, response_data, rate = .rate))
}


#' User-friendy interface to simulate_vectorized
#'
#' Does the following:
#'  * Converts factors in `newdata` to dummies with correct levels.
#'  * Checks and binds parameter values (`...`)
#'  * Call `simulate_vectorized`
#'  * Add "simulated" attribute and returns
#'
#' Parameters not documented here are documented in `pp_eval` (without dot prefix).
#'
#' @aliases simulate_atomic
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param newdata A data.frame or tibble that should contain the same columns as `fit$data`
#' @param ... Parameter values of length 1, e.g., `cp_1 = 80, Intercept_1 = -22.5` etc.
#' @return Numeric vector with attribute "simulated".
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
simulate_atomic = function(fit,
                        newdata,
                        ...,
                        .type = "predict",
                        .rate = FALSE,
                        .dpar = NULL,
                        .arma = TRUE,
                        .scale = "response") {

  # Check some inputs.
  # Remaining values are asserted in simulate_vectorized()
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_data_frame(newdata)
  args = list(...)
  model_predictors = get_fit_model_tables(fit)$predictors
  expected_args = get_sim_pars(model_predictors, fit$pars)
  if (is.null(names(args)) | any(names(args) == ""))
    stop("All arguments must be named.")
  checkmate::assert_subset(
    names(args),
    expected_args,
    .var.name = "names of arguments in `...`"
  )
  lapply(args, checkmate::assert_numeric, any.missing = FALSE)
  lapply(args, function(x) stopifnot(length(x) == 1 | length(x) == nrow(newdata)))

  # Remove response column if present - it is to be simulated
  if (fit$pars$y %in% colnames(newdata))
    newdata = dplyr::select(newdata, -dplyr::all_of(fit$pars$y))

  # Get permutations
  predictor_data = add_rhs_predictors(newdata, fit)
  pred_param_grid = cbind(predictor_data, args)  # Use tidyr::expand_grid() if any args have length > 1.
  if (.type == "predict" && .arma && is_arma(fit)) {
    values = evaluate_model_dpars(
      fit, as.list(pred_param_grid),
      paste0(".pred_", model_predictors$code_name)
    )
    warn_arma_simulation(values)
  }

  # Get y
  simulated_y = pred_param_grid %>%
    dplyr::mutate(
      .simulated_y = rlang::exec(simulate_vectorized,
                               fit,
                               !!!pred_param_grid,
                               .type = .type,
                               .rate = .rate,
                               .dpar = .dpar,
                               .arma = .arma,
                               .scale = .scale)
    ) %>%
    dplyr::pull(.data$.simulated_y)

  # add_simulated etc. and return
  attr(simulated_y, "simulated") = args  # Set as attribute
  class(attr(simulated_y, "simulated")) = c("mcplist", "list")  # for nicer printing
  simulated_y
}


#' Wrapper for `simulate_atomic()` with named arguments
#'
#' Auto-completion of named arguments  makes it easier to call from, e.g., RStudio.
#' This is typically stored as `fit$simulate()`.
#'
#' @aliases get_fitsimulate
#' @keywords internal
#' @noRd
#' @param pars A list of model parameters, typically from `fit$pars`
#' @return An R function.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_fitsimulate = function(pars) {
  # List of argument names
  sim_pars = c(pars$cp, pars$fixed[!stringr::str_ends(pars$fixed, "_sd")])  # Remove hyperparameter on varying effects from pars$reg since it is not used for simulation

  args_required = c(sim_pars, pars$sigma, pars$arma)
  args_default = c(pars$varying)
  args_all = c(args_required, args_default)

  args_withdefault = paste0(args_default, " = 0")
  args_withdefault = args_withdefault[args_withdefault != " = 0"]  # remove empty strings

  # Build function
  fitsimulate_code = paste0("function(fit, newdata, ", paste0(c(args_required, args_withdefault), collapse = ", "), ",
  .type = 'predict',
  .rate = FALSE,
  .dpar = 'epred',
  .arma = TRUE,
  .scale = 'response') {

  if (is.numeric(fit))
    stop('`fit` must be an `mcpfit` object. fit$simulate() had many breaking changes in mcp v0.4, to accomodate multiple regression models.')

  result = simulate_atomic(fit, newdata, ", paste0(args_all, " = ", args_all, collapse = ", "), ", .type = .type, .rate = .rate, .dpar = .dpar, .arma = .arma, .scale = .scale)
  return(result)
}")

  eval(parse(text = fitsimulate_code))
}


#' Simulate/evaulate autoregressive residuals
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
