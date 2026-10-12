# ABOUT: These functions take data/draws and run them through
# formula code to make model predictions
# ------------------------------------------------------------

#' Make the levels of categorical predictors match the original data
#'
#' @inheritParams add_rhs_predictors
#' @return `newdata` with all categorical columns as factors that have identical levels to the original data.
#' @noRd
relevel_newdata = function(newdata, fit) {
  # Check that the necessary data is available
  rhs_vars = get_predictor_cols(fit)
  assert_data_cols(newdata, cols = rhs_vars, fail_funcs = c(is.na, is.nan, is.infinite))

  # Make sure to carry over the exact level-structure of the original model
  for (col_name in intersect(rhs_vars, names(newdata))) {
    org_col = fit$data[, col_name]
    new_col = newdata[, col_name]

    if (is.factor(org_col)) {
      new_col = factor(
        new_col,
        levels = levels(org_col),
        ordered = is.ordered(org_col)
      )

      # Helpful error
      if (any(is.na(new_col))) {
        new_levels = newdata[is.na(new_col), col_name] %>% as.character() %>% unique()
        stop("Got novel values (", and_collapse(new_levels), ") for column ", col_name, ". Only values used during fitting are allowed.")
      }
    } else if (is.character(org_col)) {
      # Keep character type for grouping columns; factorize predictor covariates for model.matrix
      novel = setdiff(as.character(new_col), unique(org_col))
      if (length(novel) > 0) {
        stop("Got novel values (", and_collapse(novel), ") for column ", col_name, ". Only values used during fitting are allowed.")
      }
      if (col_name %in% rhs_vars) {
        new_col = factor(new_col, levels = levels(factor(org_col)))
      } else {
        new_col = as.character(new_col)
      }
    }

    newdata[, col_name] = new_col
  }

  newdata
}


# Evaluate fitted design specifications on new data.
evaluate_design_specs = function(design_specs, newdata) {
  lapply(design_specs, function(spec) {
    design_data = newdata
    if (!is.null(spec$local_x_name))
      design_data[[spec$local_x_name]] = 1
    get_fitted_design(data = design_data, spec = spec)
  })
}


# Populate matrix data columns in a predictor or group-effects table.
populate_design_columns = function(table, evaluated_designs) {
  table$matrix_data = vector("list", nrow(table))
  for (design_id in unique(stats::na.omit(table$design_id))) {
    rows = which(table$design_id == design_id)
    design = evaluated_designs[[design_id]]
    if (is.null(design))
      stop_github("Missing evaluated design '", design_id, "'.")

    component_matrix = design$matrix[, table$design_col[rows], drop = FALSE]
    table$matrix_data[rows] = unname(as.list(as.data.frame(component_matrix)))
  }

  table
}


#' Add predictors to `newdata`
#'
#' @param newdata A data.frame that contains the predictor columns from
#'   `fit$data` and `fit$model`.
#' @param fit An `mcpfit` object.
#' @return `newdata` with additional (dummy-coded) columns that make up the design matrix required for the model.
#' @noRd
add_rhs_predictors = function(newdata, fit) {
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_data_frame(newdata)

  # Make categorical predictors match the original data
  newdata = newdata %>%
    as.data.frame() %>%
    relevel_newdata(fit)

  # Evaluate each fitted design once
  model_tables = get_fit_model_tables(fit)
  predictors = model_tables$predictors
  group_effects = model_tables$group_effects
  design_specs = model_tables$design_specs
  evaluated_designs = evaluate_design_specs(design_specs, newdata)

  # Populate coefficient columns and build predictor matrix
  predictors = populate_design_columns(predictors, evaluated_designs)
  group_effects = populate_design_columns(group_effects, evaluated_designs)
  predictor_matrix = get_predictor_matrix(predictors, group_effects)

  # Collect offset columns from the evaluated designs
  offset_cols = list()
  for (id in names(design_specs)) {
    spec = design_specs[[id]]
    if (isTRUE(spec$has_offset))
      offset_cols[[spec$offset_name]] = as.numeric(evaluated_designs[[id]]$offset)
  }

  # All permutations of rows in newdata and parameters
  pred_df = as.data.frame(predictor_matrix) %>%
    magrittr::set_colnames(paste0(".pred_", colnames(predictor_matrix)))
  if (length(offset_cols) > 0)
    pred_df = dplyr::bind_cols(pred_df, as.data.frame(offset_cols))

  dplyr::bind_cols(pred_df, newdata)
}


#' Returns parameters needed for simulation
#'
#' @return Character vector
#' @noRd
get_sim_pars = function(cps, predictors, group_effects) {
  c(
    cps$name,
    predictors$name,
    group_effects$name
  )
}



#' Evaluate model and return a list of dpars
#'
#' Evaluates the R translation of `fit$.internal$formula_jags`. Parameters and
#' data are looked up in `args`, while model nodes are assigned in a separate
#' environment, so neither can overwrite the other.
#'
#' @param fit An `mcpfit` object.
#' @param args args from `simulate_vectorized`
#' @param pred_pars Design-column argument names in model-matrix order.
#' @return `data.frame` with one column per dpar
#' @noRd
evaluate_model_dpars = function(fit, args, pred_pars) {
  inputs = list2env(args, parent = topenv())
  inputs$rhs_matrix_ = do.call(cbind, args[pred_pars])  # Design columns in model-matrix order
  inputs$cp_0 = -Inf  # Outer segments extend beyond the data in R
  assign(paste0("cp_", length(fit$model)), Inf, envir = inputs)  # e.g., cp_3 = Inf

  # Evaluate model in environment
  out = new.env(parent = inputs)
  eval(str2expression(get_formula_r(fit$.internal$formula_jags)), envir = out)

  out %>%
    as.list() %>%
    as.data.frame() %>%
    dplyr::select(-dplyr::starts_with("x_local_"), -dplyr::any_of("ar_"))
}


# Draw independent normal values inside possibly different finite intervals.
# Used for the default group-level change-point hierarchy in fit$simulate().
simulate_truncated_normal = function(mean, sd, lower, upper) {
  if (sd == 0) {
    if (any(mean <= lower | mean >= upper))
      stop("A zero change-point SD places a group change point outside its ordered interval.")
    return(mean)
  }
  lower_prob = stats::pnorm(lower, mean, sd)
  upper_prob = stats::pnorm(upper, mean, sd)
  if (any(lower_prob >= upper_prob))
    stop("Could not simulate an ordered group-level change point; its truncated normal interval has negligible probability.")
  stats::qnorm(stats::runif(length(mean), lower_prob, upper_prob), mean, sd)
}


# Transform link-scale predictors to distribution-scale parameters
#
# Mirrors the deterministic dpar transformations in generated JAGS code. Both
# representations are kept only during R-side evaluation; neither
# observation-level vector is monitored in JAGS.
add_response_dpars = function(dpar_values, family) {
  for (dpar in family$dpar_specs$dpar) {
    spec = get_dpar_spec(family, dpar)
    link_name = paste0("link_", dpar, "_")
    response_name = paste0(dpar, "_")

    if (link_name %notin% names(dpar_values))
      stop_github("Missing link-scale values for dpar '", dpar, "'.")

    response = get_link_function(spec$link, inverse = TRUE)(dpar_values[[link_name]])
    if (!is.na(spec$floor))
      response = pmax(spec$floor, response)
    dpar_values[[response_name]] = response
  }

  dpar_values
}


#' Vectorized R-side run of the full model.
#'
#' Used internally in mcp.
#' Parameters not documented here are documented in `pp_eval` (without dot prefix).
#'
#' @inheritParams mcp
#' @param ... Parameter names (e.g., `cp_1 = c(4.3, 4.5, 6.2), Intercept_1 = c(11.2, 12.1, 10.9)`, etc.)
#'   and columns from `predictors$matrix_data` prefixed with ".pred_" (e.g.,
#'   `.pred_Intercept_1 = c(1, 1, 1)`).
#' @return Vector with same length as inputs.
#' @noRd
simulate_vectorized = function(fit, ..., .type = "predict", .rate = FALSE, .dpar = "epred", .arma = TRUE, .scale = "response", .return_fitted = FALSE) {
  ###########
  # ASSERTS #
  ###########
  checkmate::assert_class(fit, "mcpfit")
  model_tables = get_fit_model_tables(fit)
  predictor_occurrences = model_tables$predictors
  group_occurrences = model_tables$group_effects
  data_columns = mcp_columns(fit)

  # Assert that the ellipsis contains the expected argument names
  param_pars = get_sim_pars(model_tables$cps, definition_rows(predictor_occurrences), definition_rows(group_occurrences))
  pred_pars = paste0(
    ".pred_", get_predictor_design_names(predictor_occurrences, group_occurrences)
  )
  operation = switch(.type, fitted = "epred", loglik = "log_lik", predict = "rng")
  has_arma_terms = any(predictor_occurrences$dpar %in% c("ar", "ma"))
  aux_operations = c(operation, if (has_arma_terms && .arma) "garma")
  if (.type == "fitted" && (.rate || (!is.null(.dpar) && .dpar != "epred")))
    aux_operations = setdiff(aux_operations, "epred")
  aux_columns = get_family_aux_columns(fit$family, model_tables$segments, aux_operations)
  evaluating_arma = uses_garma(fit, .arma, .type, .dpar)
  data_pars = c(data_columns$par_x, if (evaluating_arma) data_columns$series, stats::na.omit(unname(aux_columns)))
  offset_pars = unname(unlist(lapply(
    model_tables$design_specs,
    function(s) if (isTRUE(s$has_offset)) s$offset_name else NULL
  )))
  expected_arg_names = c(param_pars, pred_pars, data_pars, offset_pars)

  args = list(...)
  missing_args = dplyr::setdiff(expected_arg_names, names(args))
  if (length(missing_args) > 0)
    stop_github("Missing the following arguments: ", and_collapse(missing_args))

  # Other args
  assert_typescale(.type, .scale)
  checkmate::assert_flag(.rate)

  .dpar = assert_dpar(.dpar, fit = fit, type = .type)
  checkmate::assert_flag(.arma)
  checkmate::assert_flag(.return_fitted)
  if (.return_fitted && .type != "predict")
    stop_github("`.return_fitted` requires `.type = 'predict'`.")


  ##################################################
  # EVALUATE MODEL AND RETURN VIA FAMILY AND .TYPE #
  ##################################################
  dpar_values = evaluate_model_dpars(fit, args, pred_pars)
  dpar_values = add_response_dpars(dpar_values, fit$family)

  # Prepare for stuff needing .ydata
  has_ydata = is.null(args[[data_columns$response]]) == FALSE
  if (has_ydata)
    dpar_values$.ydata = args[[data_columns$response]]
  if (.type == "loglik" && has_ydata == FALSE)
    stop(".ydata must be non-NULL for .type = 'loglik'.")
  .dpar = paste0(.dpar, "_")

  response_data = get_family_response_data(fit$family, model_tables$segments, args)
  dpars = stats::setNames(lapply(fit$family$dpars, function(dpar) {
    dpar_values[[paste0(dpar, "_")]]
  }), fit$family$dpars)

  # Simply return for fitted dpars
  if (.dpar %notin% c("epred_", "mu_") && .type == "fitted") {
    link_dpar = paste0("link_", .dpar)
    if (.scale == "linear" && link_dpar %in% names(dpar_values))
      return(dpar_values[[link_dpar]])
    return(dpar_values[[.dpar]])
  }

  # GARMA is defined on the link scale. The observed response is clipped only
  # where needed to keep log and logit transformations finite.
  if (has_arma_terms && .arma == TRUE) {
    base_link_mu = dpar_values$link_mu_
    ar_list = dplyr::select(dpar_values, dplyr::matches("^ar[0-9]+_$"))
    ma_list = dplyr::select(dpar_values, dplyr::matches("^ma[0-9]+_$"))
    threshold = dpar_values$garma_threshold_
    if (is.null(threshold))
      threshold = rep(0.1, length(base_link_mu))
    if (!has_ydata && .type != "predict")
      stop("The response is required to evaluate GARMA terms.")

    series_id = args[[".draw"]]
    if (!is.null(data_columns$series)) {
      if (is.null(series_id))
        series_id = rep(1, length(args[[data_columns$series]]))
      series_id = interaction(series_id, args[[data_columns$series]], drop = TRUE)
    }

    garma_result = simulate_garma(
      base_link_mu, ar_list, ma_list, threshold, fit$family,
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

  if (.type == "fitted" && .scale == "linear")
    return(dpar_values$link_mu_)

  if (.type == "fitted") {
    if (.dpar == "mu_")
      return(dpar_values$mu_)
    return(fit$family$r$epred(dpars, response_data, rate = .rate))
  }
  if (.type == "loglik")
    return(fit$family$r$log_lik(dpar_values$.ydata, dpars, response_data))
  if (.type == "predict") {
    predicted = fit$family$r$rng(length(dpars$mu), dpars, response_data, rate = .rate)
    if (.return_fitted)
      attr(predicted, "fitted") = fit$family$r$epred(dpars, response_data, rate = .rate)
    attr(predicted, "dpars") = dpars
    attr(predicted, "response_data") = response_data
    return(predicted)
  }
}


#' User-friendly interface to simulate_vectorized
#'
#' Does the following:
#'  * Converts factors in `newdata` to dummies with correct levels.
#'  * Checks and binds parameter values (`...`)
#'  * Call `simulate_vectorized`
#'  * Add "simulated" attribute and returns
#'
#' Parameters not documented here are documented in `pp_eval` (without dot prefix).
#'
#' @inheritParams mcp
#' @param newdata A data.frame or tibble that should contain the same columns as `fit$data`
#' @param ... Parameter values of length 1, e.g., `cp_1 = 80, Intercept_1 = -22.5` etc.
#' @return Numeric vector with attribute "simulated".
#' @noRd
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
  data_columns = mcp_columns(fit)
  if (uses_garma(fit, .arma, .type, .dpar))
    assert_arma_series(newdata, data_columns$series)
  args = list(...)
  model_tables = get_fit_model_tables(fit)
  model_predictor_definitions = definition_rows(model_tables$predictors)
  model_group_definitions = definition_rows(model_tables$group_effects)
  model_predictors = model_tables$predictors
  model_group_effects = model_tables$group_effects
  assert_model_data(
    newdata, data_columns$par_x,
    group_cols = stats::na.omit(model_group_effects$group_col)
  )
  assert_response_data(fit$family, model_tables$segments, newdata)
  expected_args = c(
    setdiff(get_sim_pars(model_tables$cps, model_predictor_definitions, model_group_definitions), model_group_definitions$name),
    model_group_definitions$sd_name
  )
  if (is.null(names(args)) || any(names(args) == ""))
    stop("All arguments must be named.")
  checkmate::assert_subset(
    names(args),
    expected_args,
    .var.name = "names of arguments in `...`"
  )
  lapply(args, checkmate::assert_numeric, any.missing = FALSE)
  lapply(args, function(x) stopifnot(length(x) == 1 | length(x) == nrow(newdata)))
  assert_ordered_population_cps(model_tables$cps, args)
  for (sd_name in model_group_definitions$sd_name)
    checkmate::assert_number(args[[sd_name]], lower = 0, .var.name = sd_name)

  # Remove response column if present - it is to be simulated
  if (data_columns$response %in% colnames(newdata))
    newdata = dplyr::select(newdata, -dplyr::all_of(data_columns$response))

  simulated = args
  cp_effects = model_group_definitions[model_group_definitions$part == "cp", , drop = FALSE]

  # All group-level CPs share one grouping factor. Set up its level mapping once.
  if (nrow(cp_effects) > 0) {
    cps = model_tables$cps
    group = newdata[[cp_effects$group_col[1]]]
    group_levels = unique(group)
    group_rows = match(group_levels, group)
    group_index = match(group, group_levels)
    cp_values = lapply(args[cps$name], function(value) {
      rep(value, length.out = nrow(newdata))[group_rows]
    })

    # Mirror JAGS: sample absolute locations in CP order, using a group-level
    # predecessor as the next lower bound and the next fixed CP as upper anchor.
    x_range = range(fit$data[[data_columns$par_x]], na.rm = TRUE)
    cp_locations = list()
    for (i in seq_len(nrow(cp_effects))) {
      effect = cp_effects[i, ]
      assert_generatable_group_prior(fit, effect)
      j = match(effect$population_name, cps$name)
      population = cp_values[[effect$population_name]]

      lower = if (j == 1) {
        x_range[1]
      } else {
        previous = cps[j - 1, ]
        if (!is.na(previous$group_col)) cp_locations[[previous$group_name]] else cp_values[[previous$name]]
      }

      next_fixed = which(seq_len(nrow(cps)) > j & is.na(cps$group_col))[1]
      upper = if (is.na(next_fixed)) {
        x_range[2]
      } else {
        cp_values[[cps$name[next_fixed]]]
      }

      location = simulate_truncated_normal(population, args[[effect$sd_name]], lower, upper)
      cp_locations[[effect$name]] = location
      args[[effect$name]] = (location - population)[group_index]
      simulated[[effect$name]] = args[[effect$name]]
      args[[effect$sd_name]] = NULL
    }
  }

  # Predictor group effects retain the ordinary independent normal hierarchy.
  predictor_effects = model_group_definitions[model_group_definitions$part != "cp", , drop = FALSE]
  for (i in seq_len(nrow(predictor_effects))) {
    effect = predictor_effects[i, ]
    assert_data_cols(newdata, effect$group_col)
    group = newdata[[effect$group_col]]
    group_levels = unique(group)

    # Identify alias group-level effects; otherwise add new.
    target = fit$prior[[effect$name]]
    if (!is.null(target) && target %in% names(args)) {
      args[[effect$name]] = args[[target]]
    } else {
      assert_generatable_group_prior(fit, effect)
      group_deviations = stats::rnorm(length(group_levels), 0, args[[effect$sd_name]])
      args[[effect$name]] = group_deviations[match(group, group_levels)]
    }
    simulated[[effect$name]] = args[[effect$name]]
    args[[effect$sd_name]] = NULL
  }

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

  # Get y. Internal attributes (e.g., dpars) are not returned to the user.
  simulated_y = rlang::exec(simulate_vectorized, fit, !!!pred_param_grid,
    .type = .type, .rate = .rate, .dpar = .dpar, .arma = .arma, .scale = .scale)
  attributes(simulated_y) = NULL

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
#' @param parameters Output of `get_pars_table()`.
#' @return An R function.
#' @noRd
get_fitsimulate = function(parameters) {
  # Population-level parameters in mcp_pars() order. Group-level SDs default to 0 (no group-level deviations).
  parameters = parameters[parameters$scope == "population", ]
  is_sd = parameters$role == "group_sd"
  args_required = parameters$parameter[!is_sd]
  args_pars = stats::setNames(rep(alist(x = )["x"], nrow(parameters)), parameters$parameter)
  args_pars[is_sd] = list(0)

  # Assemble formal arguments
  args_formals = c(
    alist(fit = , newdata = ),
    args_pars,
    alist(.type = "predict", .rate = FALSE, .dpar = "epred", .arma = TRUE, .scale = "response")
  )

  # Initiate with some formals to make R CMD Check happy.
  sim_fn = function(fit, newdata) {
    if (!inherits(fit, "mcpfit"))
      stop(legacy_mcp_message("`fit$simulate()` now requires the fit as its first argument. Use `fit$simulate(fit, newdata, ...)`."), call. = FALSE)
    if (missing(newdata) || !is.data.frame(newdata))
      stop(legacy_mcp_message("`fit$simulate()` now requires a data.frame or tibble as its second argument. Use `fit$simulate(fit, newdata, ...)`."), call. = FALSE)
    check_mcpfit_version(fit)
    warn_custom_jags_code(fit)

    lapply(args_required, get, envir = environment(), inherits = FALSE)
    do.call(simulate_atomic, as.list(environment(), all.names = TRUE))
  }

  formals(sim_fn) = as.pairlist(args_formals)

  # Enclose only the argument names, not the data-sized tables passed to this function
  environment(sim_fn) = list2env(list(args_required = args_required), parent = topenv())
  sim_fn
}
