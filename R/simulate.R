# ABOUT: These functions evaluate the model on data and parameter values in R, the
# counterpart of the likelihood in fit$jags_code, for predictions and fit$simulate()
# ------------------------------------------------------------

# Make categorical columns in `newdata` factors with levels identical to the original data
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
      # Factorize character predictors for model.matrix
      novel = setdiff(as.character(new_col), unique(org_col))
      if (length(novel) > 0) {
        stop("Got novel values (", and_collapse(novel), ") for column ", col_name, ". Only values used during fitting are allowed.")
      }
      new_col = factor(new_col, levels = levels(factor(org_col)))
    }

    newdata[, col_name] = new_col
  }

  newdata
}


# Populate matrix data columns in a `predictors` or `group_predictors` table.
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


#' Data used to evaluate the model on `newdata`
#'
#' @param fit An `mcpfit` object.
#' @param newdata A data frame with the predictors of the model.
#' @return A list with one value (or matrix row) per row of `newdata`:
#'   * `x`: The change-point variable.
#'   * `matrix`: The fitted design matrix, with columns in `matrix_col` order.
#'   * `offsets`: Named list of offsets, e.g., `offset_mu_1_`.
#'   * `series`: The AR/MA series, or `NULL`.
#'   * `response`: Named list of family response auxiliaries, e.g., `trials`.
#' @noRd
get_model_data = function(fit, newdata) {
  tables = get_fit_model_tables(fit)
  columns = mcp_columns(fit)
  newdata = relevel_newdata(as.data.frame(newdata), fit)

  # Evaluate each fitted design once. x-terms are multiplied by x measured from the change point later.
  designs = lapply(tables$design_specs, function(spec) {
    if (!is.null(spec$local_x_name))
      newdata[[spec$local_x_name]] = 1
    get_fitted_design(data = newdata, spec = spec)
  })
  offset_specs = Filter(function(spec) isTRUE(spec$has_offset), tables$design_specs)

  list(
    x = newdata[[columns$par_x]],
    matrix = get_predictor_matrix(
      populate_design_columns(tables$predictors, designs),
      populate_design_columns(tables$group_predictors, designs)
    ),
    offsets = stats::setNames(
      lapply(names(offset_specs), function(id) as.numeric(designs[[id]]$offset)),
      vapply(offset_specs, `[[`, character(1), "offset_name")
    ),
    series = if (!is.null(columns$series)) newdata[[columns$series]],
    response = get_family_response_data(fit$family, tables$segments, newdata)
  )
}


#' Evaluate the model for each row of data and each draw of parameter values
#'
#' Underlies `pp_eval()` and `fit$simulate()`.
#'
#' @inheritParams evaluate_formula
#' @param fit An `mcpfit` object.
#' @param y The response, one value per row or a rows x draws matrix. AR/MA terms are
#'   conditional on it. `NULL` generates fresh AR/MA series when `type = "predict"`.
#' @param type One of `"predict"`, `"fitted"`, or `"loglik"`.
#' @param rate,dpar,arma,scale See `pp_eval()`. `arma` is logical.
#' @return A list with `value`, a rows x draws matrix of the requested `type`, and `dpars`,
#'   the rows x draws distributional parameters on the response scale.
#' @noRd
evaluate_model = function(fit, data, pars, n_draws, y = NULL, type = "predict",
                          rate = FALSE, dpar = "epred", arma = TRUE, scale = "response") {
  family = fit$family
  tables = get_fit_model_tables(fit)
  links = evaluate_formula(tables, family, data, pars, n_draws)

  # Distributional parameters on the response scale, with numerical floors as in the JAGS code
  dpars = lapply(stats::setNames(nm = family$dpar_specs$dpar), function(name) {
    spec = get_dpar_spec(family, name)
    value = get_link_function(spec$link, inverse = TRUE)(links[[name]])
    if (!is.na(spec$floor))
      value[] = pmax(spec$floor, value)
    value
  })

  # Dpars other than the expected response, e.g., sigma or ar1, do not depend on AR/MA terms
  if (type == "fitted" && dpar %notin% c("epred", "mu")) {
    value = if (scale == "response" && dpar %in% names(dpars)) dpars[[dpar]] else links[[dpar]]
    return(list(value = value, dpars = dpars))
  }

  # AR/MA terms act on the link scale, conditional on `y` or generating fresh series
  generated = NULL
  if (arma && is_arma(fit)) {
    if (is.null(y) && type != "predict")
      stop("The response is required to evaluate GARMA terms.")
    threshold = unique(stats::na.omit(tables$predictors$threshold))
    garma = simulate_garma(
      links$mu, get_lags(links, "ar"), get_lags(links, "ma"), threshold, family, dpars,
      data = data$response, y = y, series = data$series
    )
    links$mu = garma$link_mu
    dpars$mu = family$linkinv(garma$link_mu)
    if (is.null(y))
      generated = garma$y
  }

  value = switch(type,
    fitted = if (scale == "linear") links$mu else if (dpar == "mu") dpars$mu else family$r$epred(dpars, data$response, rate = rate),
    loglik = (if (is.null(data$response$weights)) 1 else data$response$weights) * family$r$log_lik(y, dpars, data$response),
    predict = if (is.null(generated)) {
      family$r$rng(length(dpars$mu), dpars, data$response, rate = rate)
    } else {
      family$response$observed(generated, data$response, rate)
    }
  )
  list(value = matrix(value, nrow(links$mu), n_draws), dpars = dpars)
}


# Draw independent normal values inside possibly different finite intervals.
# Used for the default group-level change-point hierarchy in fit$simulate().
simulate_truncated_normal = function(mean, sd, lower, upper) {
  if (sd == 0) {
    if (any(mean <= lower | mean >= upper))
      stop("A zero change-point SD places a group-specific change point outside its bounds.")
    return(mean)
  }
  lower_prob = stats::pnorm(lower, mean, sd)
  upper_prob = stats::pnorm(upper, mean, sd)
  if (any(lower_prob >= upper_prob))
    stop("Could not simulate an ordered group-level change point; its truncated normal interval has negligible probability.")
  stats::qnorm(stats::runif(length(mean), lower_prob, upper_prob), mean, sd)
}


#' Simulate the model for parameter values: the body of `fit$simulate()`
#'
#' Checks the parameter values, simulates group-level deviations from their SDs, and
#' evaluates the model with `evaluate_model()`.
#'
#' Parameters not documented here are documented in `pp_eval` (without dot prefix).
#'
#' @inheritParams mcp
#' @param newdata A data.frame or tibble that should contain the same columns as `fit$data`
#' @param ... Parameter values of length 1 or `nrow(newdata)`, e.g., `cp_1 = 80, Intercept_1 = -22.5` etc.
#' @return Numeric vector with attribute "simulated".
#' @noRd
fit_simulate = function(fit,
                        newdata,
                        ...,
                        .type = "predict",
                        .rate = FALSE,
                        .dpar = NULL,
                        .arma = TRUE,
                        .scale = "response") {

  # Check inputs
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_data_frame(newdata)
  .type = rlang::arg_match0(.type, c("predict", "fitted"))
  assert_typescale(.type, .scale)
  assert_rate(.rate, fit)
  .dpar = assert_dpar(.dpar, fit = fit, type = .type)
  checkmate::assert_flag(.arma)
  data_columns = mcp_columns(fit)
  if (uses_garma(fit, .arma, .type, .dpar))
    assert_arma_series(newdata, data_columns$series)
  args = list(...)
  model_tables = get_fit_model_tables(fit)
  prior_table = fit$.internal$prior_table
  assert_model_data(
    newdata, data_columns$par_x,
    group_cols = stats::na.omit(mcp_pars(fit, scope = "group")$group_col)
  )
  assert_response_data(fit$family, model_tables$segments, newdata)
  if (is.null(names(args)) || any(names(args) == ""))
    stop("All arguments must be named.")
  checkmate::assert_subset(
    names(args),
    mcp_pars(fit, scope = "population")$parameter,
    .var.name = "names of arguments in `...`"
  )
  lapply(args, checkmate::assert_numeric, any.missing = FALSE)
  lapply(args, function(x) stopifnot(length(x) == 1 | length(x) == nrow(newdata)))
  assert_ordered_population_cps(model_tables$cps, args)
  for (sd_name in mcp_pars(fit, role = "group_sd")$parameter)
    checkmate::assert_number(args[[sd_name]], lower = 0, .var.name = sd_name)

  # Remove response column if present - it is to be simulated
  if (data_columns$response %in% colnames(newdata))
    newdata = dplyr::select(newdata, -dplyr::all_of(data_columns$response))

  simulated = args
  cps = model_tables$cps
  group_cps = which(!is.na(cps$group_col))

  # All group-level CPs share one grouping factor. Set up its level mapping once.
  if (length(group_cps) > 0) {
    group = newdata[[cps$group_col[group_cps[1]]]]
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
    for (j in group_cps) {
      cp = cps[j, ]
      assert_generatable_group_prior(fit, cp$group_name, cp$sd_name)
      population = cp_values[[cp$name]]

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

      location = simulate_truncated_normal(population, args[[cp$sd_name]], lower, upper)
      cp_locations[[cp$group_name]] = location
      args[[cp$group_name]] = (location - population)[group_index]
      simulated[[cp$group_name]] = args[[cp$group_name]]
      args[[cp$sd_name]] = NULL
    }
  }

  # Predictor group-level effects keep the ordinary independent normal hierarchy.
  predictor_effects = definition_rows(model_tables$group_predictors)
  for (i in seq_len(nrow(predictor_effects))) {
    effect = predictor_effects[i, ]
    assert_data_cols(newdata, effect$group_col)
    group = newdata[[effect$group_col]]
    group_levels = unique(group)

    # Identify alias group-level effects; otherwise add new.
    target = prior_table$value[match(effect$name, prior_table$parameter)]
    if (target %in% names(args)) {
      args[[effect$name]] = args[[target]]
    } else {
      assert_generatable_group_prior(fit, effect$name, effect$sd_name)
      group_deviations = stats::rnorm(length(group_levels), 0, args[[effect$sd_name]])
      args[[effect$name]] = group_deviations[match(group, group_levels)]
    }
    simulated[[effect$name]] = args[[effect$name]]
    args[[effect$sd_name]] = NULL
  }

  # Evaluate one draw. Parameter values per row are a one-column matrix.
  data = get_model_data(fit, newdata)
  pars = lapply(args, function(value) if (length(value) == 1) value else matrix(value))
  if (.type == "predict" && .arma && is_arma(fit))
    warn_arma_simulation(evaluate_formula(model_tables, fit$family, data, pars, n_draws = 1))
  simulated_y = as.vector(evaluate_model(
    fit, data, pars, n_draws = 1, type = .type, rate = .rate, dpar = .dpar, arma = .arma, scale = .scale
  )$value)

  # add_simulated etc. and return
  attr(simulated_y, "simulated") = simulated  # Set as attribute
  class(attr(simulated_y, "simulated")) = c("mcplist", "list")  # for nicer printing
  simulated_y
}


# Make `fit$simulate()`: a wrapper for `fit_simulate()` with the model parameters as named
# arguments, so they auto-complete in, e.g., RStudio.
# - parameters: Output of `get_pars_table()`.
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
    do.call(fit_simulate, as.list(environment(), all.names = TRUE))
  }

  formals(sim_fn) = as.pairlist(args_formals)

  # Enclose only the argument names, not the data-sized tables passed to this function
  environment(sim_fn) = list2env(list(args_required = args_required), parent = topenv())
  sim_fn
}
