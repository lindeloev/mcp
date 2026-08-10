# ABOUT: These functions take data/draws and run them through
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
  group_cols = unique(stats::na.omit(
    get_fit_model_tables(fit)$group_effects$group_col
  ))
  rhs_vars = setdiff(get_rhs_vars(fit$model), group_cols)
  assert_data_cols(newdata, cols = rhs_vars, fail_funcs = c(is.na, is.nan, is.infinite))

  # Make sure to carry over the exact level-structure of the original model
  for (col_name in intersect(c(rhs_vars, group_cols), names(newdata))) {
    org_col = fit$data[, col_name]
    new_col = newdata[, col_name]
    if (is.character(org_col) | is.factor(org_col)) {
      new_col = factor(
        new_col,
        levels = levels(factor(org_col)),
        ordered = is.ordered(org_col)
      )

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


#' Evaluate fitted predictor designs on new data
#'
#' @keywords internal
#' @noRd
#' @param table A fitted population-predictor or group-effects table.
#' @param design_specs Named fitted specifications from `get_predictor_tables()`.
#' @inheritParams add_rhs_predictors
#' @return `table` with `matrix_data` evaluated on `newdata`.
evaluate_fitted_designs = function(table, design_specs, newdata, par_x) {
  if (nrow(table) == 0 || "design_id" %notin% names(table))
    return(table)

  design_ids = unique(stats::na.omit(table$design_id))
  for (design_id in design_ids) {
    rows = which(table$design_id == design_id)
    spec = design_specs[[design_id]]
    if (is.null(spec))
      stop_github("Missing fitted design specification '", design_id, "'.")

    # Recreate the component matrix relevant here
    component_matrix = get_fitted_design(data = newdata, spec = spec)$matrix
    component_matrix = component_matrix[
      , table$design_col[rows], drop = FALSE
    ]

    # Bare par_x terms are represented relative to segment onset elsewhere in
    # mcp. Divide those factors out exactly as during fitting, then replace the
    # stored fitting-data columns with their newdata values.
    component_matrix = remove_x_factors(
      component_matrix, table$x_factor[rows], newdata[, par_x]
    )
    table$matrix_data[rows] = unname(as.list(as.data.frame(component_matrix)))
  }

  table
}


#' Rebuild predictors for fitted objects without stored design specifications
#'
#' @keywords internal
#' @noRd
#' @inheritParams add_rhs_predictors
#' @return A predictor matrix.
#' @section Removal:
#' This compatibility path can be removed for mcp 1.0 together with support
#' for fitted objects that predate stored design specifications.
get_legacy_predictor_matrix = function(newdata, fit, model_tables) {
  design_data = newdata
  group_cols = unique(stats::na.omit(model_tables$group_effects$group_col))
  for (group_col in setdiff(group_cols, names(design_data)))
    design_data[[group_col]] = fit$data[[group_col]][1]

  new_tables = get_predictor_tables(
    fit$model, design_data, fit$family, fit$pars$x, check_rank = FALSE
  )
  new_group_effects = get_group_effects(
    model_tables$cps, new_tables$group_effects
  )
  get_predictor_matrix(new_tables$predictors, new_group_effects)
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

  # Evaluate the fitted design on newdata
  model_tables = get_fit_model_tables(fit)
  predictors = model_tables$predictors
  group_effects = model_tables$group_effects
  design_specs = model_tables$design_specs
  if (is.null(design_specs)) {
    # TODO(mcp 1.0): Remove this branch with legacy fitted-object support.
    predictor_matrix = get_legacy_predictor_matrix(newdata, fit, model_tables)
  } else {
    predictors = evaluate_fitted_designs(
      predictors, design_specs, newdata, fit$pars$x
    )
    group_effects = evaluate_fitted_designs(
      group_effects, design_specs, newdata, fit$pars$x
    )
    predictor_matrix = get_predictor_matrix(predictors, group_effects)
  }

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
    pars$group
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


#' Check ordering of realized group-level change points
#'
#' @keywords internal
#' @noRd
assert_ordered_group_cps = function(cps, args) {
  if (nrow(cps) < 2 || !any(cps$varying))
    return(invisible(NULL))

  boundaries = lapply(seq_len(nrow(cps)), function(i) {
    value = args[[cps$name[i]]]
    if (cps$varying[i])
      value = value + args[[cps$group_name[i]]]
    value
  })
  n = max(lengths(boundaries))
  boundaries = vapply(
    boundaries, rep, numeric(n), length.out = n
  )
  if (any(boundaries[, -1, drop = FALSE] <= boundaries[, -ncol(boundaries), drop = FALSE]))
    stop("Realized group-level change points must remain ordered, i.e., between population-level change points. Consider using smaller `sd`.")

  invisible(NULL)
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
simulate_vectorized = function(fit, ..., .type = "predict", .rate = FALSE, .dpar = "epred", .arma = TRUE, .scale = "response", .include_fitted = FALSE) {
  ###########
  # ASSERTS #
  ###########
  checkmate::assert_class(fit, "mcpfit")
  if (!is.mcpfamily(fit$family))
    fit$family = mcpfamily(fit$family)
  model_tables = get_fit_model_tables(fit)
  predictors = model_tables$predictors
  group_effects = model_tables$group_effects

  # Assert that the ellipsis contains the expected argument names
  param_pars = get_sim_pars(predictors, fit$pars)
  pred_pars = paste0(
    ".pred_", get_predictor_design_names(predictors, group_effects)
  )
  operation = switch(.type, fitted = "epred", loglik = "log_lik", predict = "rng")
  has_arma_terms = any(predictors$dpar %in% c("ar", "ma"))
  aux_operations = c(operation, if (has_arma_terms && .arma) "garma")
  aux_columns = get_family_aux_columns(fit$family, model_tables$segments, aux_operations)
  data_pars = c(fit$pars$x, fit$pars$series, stats::na.omit(unname(aux_columns)))
  expected_arg_names = c(param_pars, pred_pars, data_pars)

  args = list(...)
  missing_args = dplyr::setdiff(expected_arg_names, names(args))
  if (length(missing_args) > 0)
    stop_github("Missing the following arguments: ", and_collapse(missing_args))
  assert_ordered_group_cps(model_tables$cps, args)

  # Other args
  assert_typescale(.type, .scale)
  checkmate::assert_flag(.rate)

  .dpar = assert_dpar(.dpar, fit = fit, type = .type)
  checkmate::assert_flag(.arma)
  checkmate::assert_flag(.include_fitted)
  if (.include_fitted && .type != "predict")
    stop_github("`.include_fitted` requires `.type = 'predict'`.")


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
  if (has_arma_terms && .arma == TRUE) {
    base_link_mu = if (uses_link_dpars) dpar_values$link_mu_ else fit$family$linkfun(dpar_values$mu_)
    ar_list = dplyr::select(dpar_values, dplyr::matches("^ar[0-9]+_$"))
    ma_list = dplyr::select(dpar_values, dplyr::matches("^ma[0-9]+_$"))
    boundary = dpar_values$garma_boundary_
    if (is.null(boundary))
      boundary = rep(0.1, length(base_link_mu))
    if (!has_ydata && .type != "predict")
      stop("The response is required to evaluate GARMA terms.")

    series_id = args[[".draw"]]
    if (!is.null(fit$pars$series)) {
      if (is.null(series_id))
        series_id = rep(1, length(args[[fit$pars$series]]))
      series_id = interaction(series_id, args[[fit$pars$series]], drop = TRUE)
    }

    garma_result = simulate_garma(
      base_link_mu, ar_list, ma_list, boundary, fit$family,
      dpars = dpars, data = response_data,
      y = if (has_ydata) dpar_values$.ydata else NULL,
      series_id = series_id
    )
    if (!has_ydata)
      return(fit$family$response$observed(garma_result$y, response_data, .rate))

    dpar_values$link_mu_ = garma_result$link_mu
    dpar_values$mu_ = garma_result$mu
    dpars$mu = garma_result$mu
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
  if (.type == "predict") {
    predicted = fit$family$r$rng(length(dpars$mu), dpars, response_data, rate = .rate)
    if (.include_fitted)
      attr(predicted, "fitted") = fit$family$r$epred(dpars, response_data, rate = .rate)
    return(predicted)
  }
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
  assert_arma_series(newdata, fit$pars$series)
  args = list(...)
  model_tables = get_fit_model_tables(fit)
  model_predictors = model_tables$predictors
  model_group_effects = model_tables$group_effects
  expected_args = c(
    setdiff(get_sim_pars(model_predictors, fit$pars), model_group_effects$name),
    model_group_effects$sd_name
  )
  if (is.null(names(args)) | any(names(args) == ""))
    stop("All arguments must be named.")
  checkmate::assert_subset(
    names(args),
    expected_args,
    .var.name = "names of arguments in `...`"
  )
  lapply(args, checkmate::assert_numeric, any.missing = FALSE)
  lapply(args, function(x) stopifnot(length(x) == 1 | length(x) == nrow(newdata)))
  for (sd_name in model_group_effects$sd_name)
    checkmate::assert_number(args[[sd_name]], lower = 0, .var.name = sd_name)

  # Remove response column if present - it is to be simulated
  if (fit$pars$y %in% colnames(newdata))
    newdata = dplyr::select(newdata, -dplyr::all_of(fit$pars$y))

  # Simulate one deviation per grouping level, then map to rows. Change-point
  # deviations are exactly centered, matching the fitted parameterization.
  simulated = args
  for (i in seq_len(nrow(model_group_effects))) {
    effect = model_group_effects[i, ]
    assert_data_cols(newdata, effect$group_col)
    group = newdata[[effect$group_col]]
    group_levels = unique(group)
    group_deviations = stats::rnorm(length(group_levels), 0, args[[effect$sd_name]])
    if (effect$part == "cp")
      group_deviations = group_deviations - mean(group_deviations)
    args[[effect$name]] = group_deviations[match(group, group_levels)]
    simulated[[effect$name]] = args[[effect$name]]
    args[[effect$sd_name]] = NULL
  }
  assert_ordered_group_cps(model_tables$cps, args)

  # Get permutations
  predictor_data = add_rhs_predictors(newdata, fit)
  pred_param_grid = cbind(predictor_data, args)  # Use tidyr::expand_grid() if any args have length > 1.
  if (.type == "predict" && .arma && is_arma(fit)) {
    values = evaluate_model_dpars(
      fit, as.list(pred_param_grid),
      paste0(
        ".pred_",
        get_predictor_design_names(model_predictors, model_group_effects)
      )
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
  attr(simulated_y, "simulated") = simulated  # Set as attribute
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
#' @param pars A list of model parameters, typically from `fit$pars`.
#' @param group_effects The output of `get_group_effects()`.
#' @return An R function.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_fitsimulate = function(pars, group_effects) {
  # List of argument names
  sim_pars = c(pars$cp, pars$mu)
  predictor_group_effects = group_effects[group_effects$part == "predictor", , drop = FALSE]
  cp_group_effects = group_effects[group_effects$part == "cp", , drop = FALSE]

  args_required = c(sim_pars, pars$sigma, pars$arma)
  args_default = c(cp_group_effects$sd_name, predictor_group_effects$sd_name)
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

  if (!inherits(fit, 'mcpfit'))
    stop('`fit` must be an `mcpfit` object. fit$simulate() had breaking changes in mcp v0.4.0. Signature is fit$simulate(fit, newdata, ...).', call. = FALSE)
  if (missing(newdata) || !is.data.frame(newdata))
    stop('`newdata` must be a data.frame or tibble. fit$simulate() had breaking changes in mcp v0.4.0. Signature is fit$simulate(fit, newdata, ...).', call. = FALSE)

  result = simulate_atomic(fit, newdata, ", paste0(args_all, " = ", args_all, collapse = ", "), ", .type = .type, .rate = .rate, .dpar = .dpar, .arma = .arma, .scale = .scale)
  return(result)
}")

  eval(parse(text = fitsimulate_code))
}
