# ABOUT: Posterior predictions, fitted values, residuals, and log-likelihood.
# --------------------------------------------------------------------------

#' Fits and Predictions given Draws and data
#'
#' @aliases pp_eval pp_eval.mcpfit
#' @keywords internal
#' @param varying Group-level effects. One of:
#'   * `TRUE` All group-level deviations.
#'   * `FALSE` No group-level deviations (`c()`).
#'   * `"cp"` or `"predictor"`: All group-level deviations belonging to that part of
#'     the model.
#'   * Character vector: Only include specified group-level parameters.
#' @param object An `mcpfit` object.
#' @param newdata A `tibble` or a `data.frame` containing predictors in the model. 
#' - If `NULL` (default), the original data is used.
#' - For models with `ar()` or `ma()`: `fitted()`, `residuals()`, `log_lik()`,
#'   and `predict()` condition on the response history by default,
#'   so `newdata` must include the response. For `fitted()`, `predict()`, and
#'   `residuals()`, missing response histories are supported only in the original
#'   fitted data, using retained posterior imputations as histories. Predictions
#'   are fresh response draws, including at missing rows. With `conditional = FALSE`,
#'   `predict()` and [`posterior_predict()`][rstantools::posterior_predict] generate
#'   fresh response series recursively, so their `newdata` need only contain
#'   predictors and required response auxiliaries. `log_lik()`
#'   is unavailable when a missing response enters a later observed history.
#' - For models with `y | weights()`: Require the weights column except for `fitted()` and `predict()`.
#' @param summary Summarise at each x-value
#' @param type One of:
#'   - `"fitted"`: return the expected response. When `dpar` names a
#'     distributional parameter (e.g., `"mu"` or `"sigma"`), that parameter is returned instead.
#'     See also `fitted()`.
#'   - `"predict"`: return predicted values (e.g., `y_predict = rnorm(N, y_fitted, sigma_fitted)` for `family = gaussian()`).
#'     See also `predict()`.
#'   - `"residuals"`: observed y-values minus the fitted values. See also `residuals()`.
#'   - `"loglik"`: return the log-likelihood for each draw for each data point. See also `log_lik()`.
#'     Requires `scale = "response"`.
#' @param probs Vector of quantiles (strictly between 0 and 1). Only in effect when `summary == TRUE`.
#' @param rate Logical scalar. For binomial models, return counts (`rate = FALSE`, the default for `fitted()` and `predict()`) or
#'   the observed or expected success proportion (`rate = TRUE`). Predictions and
#'   count-scale fitted values require a trials column in `newdata`.
#'   Distributional parameters such as `dpar = "mu"` evaluate the parameter itself (e.g., success probability)
#'   and are unaffected by `rate`.
#' @param prior Logical. Evaluate prior draws (`TRUE`) instead of posterior draws (`FALSE`, default).
#'   The selected draws must be available; prior-only fits require `prior = TRUE`.
#' @param dpar What distributional parameter to evaluate. This is only relevant when `type == "fitted"`. E.g.,
#'
#'   * `"epred"` (default): Expected response from the full model (or `NULL` for compatibility with brms etc.).
#'   * `"mu"`: The conditional mean (or success probability per trial for binomial/bernoulli models), on the link or response scale.
#'   * `"sigma"`: The standard deviation of the residuals.
#'   * `"ar1"`, `"ar2"`, `"ma1"`, `"ma2"`, etc. depending on which AR or MA
#'     coefficient you want to evaluate.
#' @param arma Whether to include AR and MA effects.
#'   * `TRUE` Compute the GARMA residual recurrence. Requires the response variable in `newdata`.
#'   * `FALSE` Disregard AR and MA effects. For `family = gaussian()`, `predict()` uses only `sigma` for residuals.
#'   For posterior evaluation of the original data, retained JAGS imputations
#'   supply missing GARMA histories. In models with group-level effects, this
#'   currently requires including all such effects (`varying = TRUE`).
#' @param ndraws Integer or `NULL`. Number of posterior draws to return/summarise.
#'   If there are group-level effects, this is the number of draws from each group.
#'   `NULL` means "all". More draws trade speed for accuracy.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @param draws_format One of "tidy" or "matrix". Controls the output format when `summary == FALSE` (for `fitted()`, `predict()`, and `log_lik()`). `residuals()` always returns tidy output.
#' @param samples_format Deprecated. Use `draws_format` instead.
#'   See more under "value"
#' @param scale One of
#'   * `"response"`: return on the observed scale, i.e., after applying the inverse link function.
#'   * `"linear"`: return on the linear-predictor (link) scale, where the linear
#'     trends are modeled.
#'     A linear scale is only applicable when `type == "fitted"` and `dpar` is not `NULL`.
#' @param .include_fitted Internal. Include fitted values with unsummarised predictions.
#' @param .include_dpars Internal. Include distributional parameters and response data as attributes with unsummarised predictions.
#' @param .garma_replicate Internal. For GARMA predictions, generate each
#'   response history recursively instead of conditioning on observed responses.
#' @return
#'   * If `summary = TRUE`: A data frame with the draw mean and SD (`sd`) for
#'     each row in `newdata`. With posterior draws (the default), `sd` is the
#'     posterior predictive SD for `type = "predict"` and the posterior SD of the
#'     evaluated quantity otherwise. With `prior = TRUE`, these are the analogous
#'     prior summaries. If `newdata` is `NULL`, the data in `fit$data` is used.
#'
#'   * If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes` `tibble` with all the posterior
#'     draws (`Nd`) evaluated at each row in `newdata` (`Nn`), i.e., with `Nd x Nn` rows. If there are
#'     group-level effects, the returned data is expanded with the relevant levels for each row.
#'
#'     The return columns are:
#'
#'      - Predictors from `newdata`, plus its response column when supplied.
#'      - Draw descriptors: ".chain", ".iteration", ".draw" (see the `posterior` and `tidybayes` packages), and `data_row`, the row number in the evaluated `newdata`.
#'      - Draw values: one column for each parameter in the model.
#'      - The estimate. Either ".epred", ".prediction", ".residual", or ".loglik" (matching tidybayes/ggdist conventions).
#'
#'   * If `summary = FALSE` and `draws_format = "matrix"`: An `N_draws` X `nrows(newdata)` matrix with fitted/predicted
#'       values (depending on `type`). This format is used by `brms` and it's useful as `yrep` in
#'      `bayesplot::ppc_*` functions.
#' @seealso \code{\link{fitted.mcpfit}} \code{\link{predict.mcpfit}} \code{\link{residuals.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
pp_eval = function(
  object,
  newdata = NULL,
  summary = TRUE,
  type = "fitted",
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  dpar = "epred",
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = 'response',
  .include_fitted = FALSE,
  .include_dpars = FALSE,
  .garma_replicate = FALSE,
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated()
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "pp_eval")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "pp_eval")

  # Recode
  fit = object
  checkmate::assert_class(fit, "mcpfit")
  warn_custom_jags_code(fit)
  if (!is.mcpfamily(fit$family))
    fit$family = mcpfamily(fit$family)
  if (is.null(fit$family$r$cdf))
    fit$family$r$cdf = mcpfamily(fit$family)$r$cdf
  dpar = assert_dpar(dpar, fit = fit, type = type)

  # What data to use
  using_original_data = is.null(newdata) || identical(data.frame(newdata), fit$data)
  if (using_original_data)
    newdata = fit$data

  data_columns = mcp_columns(fit)
  checkmate::assert_flag(.garma_replicate)
  replicate_garma = .garma_replicate && arma && is_arma(fit)
  assert_arma_series(newdata, data_columns$series)
  if (type == "loglik")
    assert_loglik_garma_history(fit, newdata, arma)

  conditional_garma = arma && is_arma(fit) && !replicate_garma &&
    (type %in% c("predict", "residuals") ||
       (type == "fitted" && dpar %in% c("epred", "mu")))
  if (conditional_garma && !using_original_data &&
      data_columns$response %in% names(newdata) && anyNA(newdata[[data_columns$response]]))
    stop(
      "Conditional GARMA evaluation with missing responses is supported only for ",
      "the original fitted data, where retained posterior imputations are available. ",
      "Use `posterior_predict(..., conditional = FALSE)` to generate fresh replicated response series from ",
      "predictor-only `newdata`.",
      call. = FALSE
    )

  response_return = if (data_columns$response %in% colnames(newdata))
    newdata[, data_columns$response, drop = FALSE] else NULL


  ###############
  # FIX NEWDATA #
  ###############
  # Identify grouping columns to exclude based on the `varying` argument
  group_info = unpack_group_effects(fit, pars = varying)
  model_tables = get_fit_model_tables(fit)
  group_cols = unique(stats::na.omit(model_tables$group_effects$group_col))
  exclude_group_cols = setdiff(group_cols, c(group_info$cols, data_columns$series))

  # Determine which auxiliary columns are needed for this operation
  operation = switch(type, predict = "rng", loglik = "log_lik", fitted = "epred", residuals = "epred")
  aux_operations = c(operation, if (arma && is_arma(fit)) "garma")
  if (type == "fitted" && (rate || dpar != "epred"))
    aux_operations = setdiff(aux_operations, "epred")
  aux_columns = get_family_aux_columns(fit$family, model_tables$segments)
  aux_used = names(get_family_aux_columns(fit$family, model_tables$segments, aux_operations))
  unused_aux_columns = unname(aux_columns[names(aux_columns) %notin% aux_used])

  # Build list of required columns and validate presence in newdata
  required_cols = colnames(fit$data)  # Only predictive columns were saved in fit$data
  required_cols = required_cols[required_cols %notin% unused_aux_columns]
  required_cols = required_cols[required_cols %notin% exclude_group_cols]
  response_not_required = type %in% c("fitted", "predict") && !conditional_garma
  if (response_not_required) {
    required_cols = required_cols[required_cols != data_columns$response]
  } else if (data_columns$response %notin% colnames(newdata)) {
    stop("`newdata` must contain a response column named '", data_columns$response, "' for when `arma == TRUE` and/or `type == 'residuals'`")
  }
  assert_data_cols(newdata, required_cols)
  assert_response_data(
    fit$family,
    model_tables$segments,
    newdata,
    response_required = !response_not_required,
    aux_required = aux_used
  )

  # Validate against reserved output namespace
  assert_reserved_output_namespace(colnames(newdata), context = "newdata")

  # Filter newdata columns and attach unique evaluation row index
  kept_cols = colnames(newdata)[colnames(newdata) %notin% exclude_group_cols]
  if (replicate_garma)
    kept_cols = kept_cols[kept_cols != data_columns$response]
  newdata = data.frame(newdata[, kept_cols, drop = FALSE])
  newdata$.mcp_data_row = seq_len(nrow(newdata))  # Evaluation key throughout summaries, matrices, plots, and metrics
  newdata_return = newdata
  if (!is.null(response_return) && data_columns$response %notin% colnames(newdata_return))
    newdata_return[[data_columns$response]] = response_return[[data_columns$response]]

  ########################
  # ASSERTS AND RECODING #
  ########################
  checkmate::assert_flag(summary)
  assert_typescale(type, scale)
  checkmate::assert(
    checkmate::check_flag(probs),
    checkmate::check_numeric(probs, any.missing = FALSE),
    .var.name = "probs"
  )
  if (is.numeric(probs)) {
    checkmate::assert_numeric(probs, min.len = 1, any.missing = FALSE)
    if (any(probs <= 0 | probs >= 1))
      stop("`probs` must be strictly between 0 and 1.")
  }
  if (is.logical(probs) && all(probs == TRUE))
    probs = c(0.025, 0.975)
  checkmate::assert_flag(rate)
  checkmate::assert_flag(prior)
  checkmate::assert_flag(arma)
  checkmate::assert_flag(.include_fitted)
  checkmate::assert_flag(.include_dpars)
  if (.include_fitted && (type != "predict" || summary))
    stop_github("`.include_fitted` requires `type = 'predict'` and `summary = FALSE`.")
  if (.include_dpars && (type != "predict" || summary))
    stop_github("`.include_dpars` requires `type = 'predict'` and `summary = FALSE`.")
  if (.garma_replicate && type != "predict")
    stop_github("`.garma_replicate` requires `type = 'predict'`.")
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)


  ########################
  # GET FITS/PREDICTIONS #
  ########################
  simulate_type = ifelse(type == "residuals", yes = "fitted", no = type)
  if (length(group_info$cols) > 0) {
    # Match group-level draws to each row of data.
    draws = dplyr::left_join(
      add_rhs_predictors(newdata, fit),
      mcp_draws(fit, population = TRUE, varying = varying, prior = prior, ndraws = ndraws),
      by = unique(group_info$cols),
      relationship = "many-to-many"
    )
  } else {
    # Without group-level effects, use all draws for each row of data.
    mcmc_draws = tibble::as_tibble(mcp_draws(fit, population = TRUE, varying = varying, prior = prior, ndraws = ndraws))
    predictors = tibble::as_tibble(add_rhs_predictors(newdata, fit))
    draws = dplyr::cross_join(mcmc_draws, predictors)
  }

  # Complete the conditional history, keeping imputations out of prediction output.
  if (conditional_garma && anyNA(draws[[data_columns$response]])) {
    if (prior)
      stop("Missing GARMA histories require posterior draws. Use `conditional = FALSE` for prior prediction.", call. = FALSE)
    if (!all(model_tables$group_effects$name %in% group_info$pars))
      stop(
        "This model has group-level effects, and its retained missing-response ",
        "histories are conditional on all of them. GARMA evaluation with missing ",
        "responses therefore currently requires `varying = TRUE`.",
        call. = FALSE
      )
    if (is.null(fit$.internal$imputed_response))
      stop(
        "This fit does not retain the missing response draws needed for coherent ",
        "GARMA evaluation. ",
        if (has_custom_jags_code(fit)) {
          "Automatic response imputation is unavailable with custom `jags_code`."
        } else {
          "Refit the model with the current version of mcp."
        },
        call. = FALSE
      )
    missing_rows = is.na(draws[[data_columns$response]])
    draws[[data_columns$response]][missing_rows] = get_imputed_response_draws(fit, draws)[missing_rows]
  }

  # This is the important step!Evaluate the mcp model on newdata and draws.
  # Group-level joins are row-major, while GARMA recurrences require each
  # draw's data rows to be contiguous. Evaluate in draw/data order, then
  # restore the public row order below.
  evaluation_order = if (arma && is_arma(fit)) {
    order(draws$.draw, draws$.mcp_data_row)
  } else {
    seq_len(nrow(draws))
  }
  evaluation_data = draws[evaluation_order, , drop = FALSE]
  evaluate = function() rlang::exec(simulate_vectorized, fit, !!!evaluation_data, .type = simulate_type, .rate = rate, .dpar = dpar, .arma = arma, .scale = scale, .include_fitted = .include_fitted)
  evaluated = if (replicate_garma) suppressMessages(evaluate()) else evaluate()

  # Now more boilerplate stuff...
  fitted_values = attr(evaluated, "fitted")
  dpars_values = attr(evaluated, "dpars")
  response_data_values = attr(evaluated, "response_data")
  attr(evaluated, "fitted") = NULL
  attr(evaluated, "dpars") = NULL
  attr(evaluated, "response_data") = NULL
  restore_order = order(evaluation_order)
  evaluated = evaluated[restore_order]
  if (!is.null(fitted_values)) fitted_values = fitted_values[restore_order]
  if (!is.null(dpars_values)) dpars_values = lapply(dpars_values, function(v) v[restore_order])
  if (!is.null(response_data_values)) response_data_values = lapply(response_data_values, function(v) v[restore_order])
  draws[[type]] = evaluated

  # Plotting can request fitted and predicted values from the same evaluated
  # parameter rows and model evaluation.
  if (.include_fitted)
    draws$fitted = fitted_values

  if (!is.null(response_return))
    draws[[data_columns$response]] = response_return[[data_columns$response]][draws$.mcp_data_row]

  draws = draws %>% dplyr::select(-dplyr::starts_with(".pred_"))

  # Missing outcomes are latent in the fitted JAGS model, but they are not
  # observed-data likelihood contributions. Retain them while evaluating
  # GARMA histories above, then remove them from returned log likelihoods.
  if (type == "loglik") {
    observed_rows = which(!is.na(newdata[, data_columns$response]))
    if (length(observed_rows) == 0)
      stop("Log-likelihood evaluation requires at least one observed response.")
    draws = dplyr::filter(draws, .data$.mcp_data_row %in% observed_rows)
    newdata_return = dplyr::filter(
      newdata_return,
      .data$.mcp_data_row %in% observed_rows
    )
  }


  # Optionally compute residuals
  if (type == "residuals")
    draws = dplyr::mutate(draws, !!type := .data[[data_columns$response]] - .data[[type]])

  # Fail early if group-level joins or another evaluation step duplicated
  # or dropped any joint draw/evaluation-row combinations.
  validate_eval_draws(draws, type)

  # Optionally summarise
  if (summary == TRUE) {
    df_return = draws %>%
      # Summarise for each row in newdata
      dplyr::group_by(.data$.mcp_data_row) %>%
      dplyr::summarise(.groups = "drop",
                       sd = stats::sd(.data[[type]]),
                       !!type := mean(.data[[type]])
      ) %>%

      # Apply original order and put newdata as the first columns
      dplyr::arrange(.data$.mcp_data_row) %>%
      dplyr::left_join(newdata_return, by = ".mcp_data_row", relationship = "one-to-one") %>%
      dplyr::select(dplyr::one_of(colnames(newdata_return)), dplyr::all_of(type), "sd")


    # Quantiles
    if (!isFALSE(probs)) {
      val_col = if (type == "predict" && !is.null(fit$family$r$cdf)) ".predicted" else type
      quantiles = if (type == "predict" && !is.null(fit$family$r$cdf)) {
        get_mixture_quantiles(draws, probs, fit$family, keep = NULL, rate = rate, dpars = dpars_values, response_data = response_data_values)
      } else {
        get_quantiles(draws, probs, type, na.rm = type == "residuals")
      }
      quantiles = quantiles %>%
        dplyr::mutate(quantile = 100 * .data$quantile) %>%
        tidyr::pivot_wider(names_from = "quantile", names_prefix = "Q", values_from = dplyr::all_of(val_col))

      df_return = dplyr::left_join(df_return, quantiles, by = ".mcp_data_row", relationship = "one-to-one")
    }
    return(data.frame(dplyr::select(df_return, -".mcp_data_row")))
  } else if (draws_format == "tidy") {
    value_col = switch(type,
      fitted = ".epred",
      predict = ".prediction",
      residuals = ".residual",
      loglik = ".loglik",
      type
    )
    if (.include_fitted && "fitted" %in% colnames(draws)) {
      draws = dplyr::rename(draws, .epred = "fitted")
    }
    draws = dplyr::rename(draws, !!value_col := dplyr::all_of(type))
    if (.include_dpars) {
      if (!is.null(dpars_values)) attr(draws, "dpars") = dpars_values
      if (!is.null(response_data_values)) attr(draws, "response_data") = response_data_values
    }
    draws$data_row = draws$.mcp_data_row
    draws$.mcp_data_row = NULL
    return(draws)
  } else if (draws_format == "matrix") {
    df_return = tidy_to_matrix(draws, type)
    return(df_return)
  }
}



#' Fitted and predicted values of `mcp` models fits
#'
#' Evaluate the model on data, either summarised (per data-row) or per draw. You
#' can use draws from the prior (`prior = TRUE`), select a distributional
#' parameter with `dpar`, and choose the response or linear-predictor scale with
#' `scale` where applicable.
#'
#' @details
#' `fitted()` and `posterior_epred()` evaluate the same expected responses;
#' `predict()` and `posterior_predict()` evaluate the same response distributions.
#' The `posterior_*()` methods return draws-by-observation matrices, while
#' `fitted()` and `predict()` summarise by default and also offer tidy draws.
#' For binomial models, the default response scale is counts. Use `rate = TRUE`
#' for proportions or `fitted(..., dpar = "mu")` for success probabilities.
#' During migration from v0.3.4, an omitted `rate` warns once per session per
#' function when counts differ from proportions. Explicit `rate` settings do not warn.
#'
#' `residuals(fit)` is equivalent to `fit$data[[mcp_columns(fit)$response]] - fitted(fit, ...)` (or `newdata[[mcp_columns(fit)$response]] - fitted(fit, ...)`),
#' but with fixed arguments for `fitted`: `rate = FALSE, dpar = 'epred', draws_format = 'tidy'`.
#'
#' `log_lik()` defaults to an unsummarised draws-by-observation matrix, as used
#' by `loo` and other posterior workflows. Non-default `varying` and `arma`
#' settings evaluate conditional or counterfactual log-likelihoods (e.g.,
#' omitting random effects or serial correlation); they cannot be used in
#' `loo()` or `waic()` because estimating information criteria for reduced
#' models requires refitting.
#'
#' Missing responses in the original data remain missing in the response column.
#' `fitted()` returns their expected responses, while `predict()` uses retained
#' JAGS imputations for their posterior response draws. In GARMA models these
#' imputations also supply the history used for later fitted and predicted rows.
#'
#' @inheritParams pp_eval
#' @param ... Must be empty. Reserved for future use.
#' @inherit pp_eval return
#' @seealso \code{\link{fitted.mcpfit}} \code{\link{predict.mcpfit}} \code{\link{residuals.mcpfit}} \code{\link{log_lik.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' head(fitted(demo_fit))  # Expected response for each row of demo_fit$data
#' head(residuals(demo_fit))  # Residuals for each row of demo_fit$data
#' log_lik(demo_fit)[1:3, 1:3]  # Log-likelihood at each demo_fit$data
#'
#' # All of the above take a range of arguments. E.g.,:
#' \donttest{
#' head(predict(demo_fit))  # Pointwise posterior predictive
#' head(predict(demo_fit, probs = c(0.1, 0.5, 0.9)))  # Median and 80% posterior predictive interval.
#' head(predict(demo_fit, prior = TRUE))  # Prior predictive
#' head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#' head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#'
#' # Evaluate at novel data
#' novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
#' head(predict(demo_fit, newdata = novel_data, probs = c(0.025, 0.5, 0.975)))
#'
#' # Work with missing responses
#' missing_fit = mcp_example("missing", plot = FALSE)
#' fitted(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Expected responses for missing y
#' fitted(missing_fit, summary = FALSE) |> dplyr::filter(is.na(y)) |> head()  # Same, but draws
#' predict(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Posterior predictive for missing y
#'}
#' @name execute-mcp-model
NULL


#' @aliases predict predict.mcpfit
#' @describeIn execute-mcp-model Predictive Distribution
#' @param conditional Logical. For AR/MA models, condition on observed response
#'   histories (`TRUE`, default) or generate fresh histories recursively (`FALSE`).
#'   Applies equally to prior and posterior draws. Predictive checks
#'   with `pp_check()` generate fresh histories.
#' @export
predict.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = FALSE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  conditional = TRUE,
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "predict.mcpfit")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "predict.mcpfit")
  dots = list(...)
  warn_which_y(dots, "predict")
  dots$which_y = NULL
  if (length(dots) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(names(dots)), call. = FALSE)

  checkmate::assert_flag(conditional)
  if (missing(rate))
    warn_binomial_rate(object, newdata, "predict")

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "predict",
    probs = probs,
    rate = rate,
    prior = prior,
    dpar = NULL,
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    draws_format = draws_format,
    .garma_replicate = !conditional
  )
}


#' @aliases fitted fitted.mcpfit
#' @describeIn execute-mcp-model Expected response
#' @export
fitted.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = FALSE,
  prior = FALSE,
  dpar = "epred",
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = "response",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "fitted.mcpfit")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "fitted.mcpfit")
  dots = list(...)
  warn_which_y(dots, "fitted")
  if ("which_y" %in% names(dots) && missing(dpar))
    dpar = dots$which_y
  dots$which_y = NULL
  if (length(dots) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(names(dots)), call. = FALSE)

  if (missing(rate) && (is.null(dpar) || identical(dpar, "epred")) && scale == "response")
    warn_binomial_rate(object, newdata, "fitted")

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "fitted",
    probs = probs,
    rate = rate,
    prior = prior,
    dpar = dpar,
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    draws_format = draws_format,
    scale = scale
  )
}


#' Posterior prediction draws for `mcpfit` objects
#'
#' Methods for the `{rstantools}` posterior-prediction generics. They return a
#' draws-by-observation matrix and enable `{tidybayes}` workflows such as
#' `add_epred_draws()`, `add_predicted_draws()`, and `add_linpred_draws()`.
#' These methods and workflows require the suggested package `{rstantools}`.
#'
#' @param object An `mcpfit` object.
#' @inheritParams pp_eval
#' @inheritParams predict.mcpfit
#' @param draws,ndraws Number of posterior draws to return. `draws` follows the
#'   `{rstantools}` convention; `ndraws` is the mcp spelling. Supply at most one.
#' @param re.form,re_formula Group-level effects to include. `NULL` includes all
#'   effects and `NA` excludes them.
#' @param dpar Distributional parameter for `posterior_epred()` and
#'   `posterior_linpred()`; `NULL` uses the expected response.
#' @param transform For `posterior_linpred()`, return the inverse-link
#'   transformed expected response instead of the linear predictor.
#' @param seed Optional integer seed for draw selection and posterior prediction.
#' @param ... Must be empty. Reserved for future use.
#' @return A numeric `N_draws` by `nrow(newdata)` matrix.
#' @details For GARMA models, `posterior_predict()` conditions on the observed
#'   response history, just like `predict()`. Use `conditional = FALSE` in either
#'   method to generate fresh response histories recursively. Missing responses
#'   in the original data are filled with retained imputations only to supply
#'   histories. Conditional `posterior_predict()` draws from the response
#'   distributions whose means `posterior_epred()` returns, including at missing
#'   rows. These methods require posterior draws. For prior prediction, use
#'   `predict()` with `prior = TRUE`; `conditional` selects the same behavior.
#'
#'   For binomial models, `posterior_epred()` and `posterior_predict()` (and
#'   corresponding `{tidybayes}` workflows such as `add_epred_draws()`) follow
#'   `{brms}` and `{rstantools}` conventions by returning values on the outcome
#'   count scale (`rate = FALSE`), i.e., expected counts \eqn{E[Y] = n\mu} and
#'   simulated counts in \eqn{\{0, \dots, n\}}, matching `fitted()` and `predict()`.
#'   Use `rate = TRUE` for proportions. To obtain the success probability
#'   parameter \eqn{\mu} on the \eqn{[0, 1]} scale regardless of trial counts, pass
#'   `dpar = "mu"`.
#' @seealso [fitted.mcpfit()], [predict.mcpfit()]
posterior_epred.mcpfit = function(
  object,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  dpar = NULL,
  seed = NULL,
  rate = FALSE,
  ...
) {
  posterior_prediction_matrix(
    object = object,
    newdata = newdata,
    type = "fitted",
    draws = draws,
    ndraws = ndraws,
    re.form = re.form,
    re_formula = re_formula,
    dpar = dpar,
    scale = "response",
    rate = rate,
    seed = seed,
    ...
  )
}


#' @rdname posterior_epred.mcpfit
posterior_predict.mcpfit = function(
  object,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  seed = NULL,
  rate = FALSE,
  conditional = TRUE,
  ...
) {
  posterior_prediction_matrix(
    object = object,
    newdata = newdata,
    type = "predict",
    draws = draws,
    ndraws = ndraws,
    re.form = re.form,
    re_formula = re_formula,
    seed = seed,
    rate = rate,
    conditional = conditional,
    ...
  )
}


#' @rdname posterior_epred.mcpfit
posterior_linpred.mcpfit = function(
  object,
  transform = FALSE,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  dpar = NULL,
  seed = NULL,
  ...
) {
  checkmate::assert_flag(transform)
  posterior_prediction_matrix(
    object = object,
    newdata = newdata,
    type = "fitted",
    draws = draws,
    ndraws = ndraws,
    re.form = re.form,
    re_formula = re_formula,
    dpar = dpar,
    scale = if (transform) "response" else "linear",
    rate = TRUE,
    seed = seed,
    ...
  )
}


# Evaluate posterior prediction draws for rstantools-compatible methods
posterior_prediction_matrix = function(
  object,
  newdata,
  type,
  draws,
  ndraws,
  re.form,
  re_formula,
  dpar = NULL,
  scale = "response",
  rate = FALSE,
  seed = NULL,
  conditional = TRUE,
  ...
) {
  checkmate::assert_flag(conditional)
  checkmate::assert_class(object, "mcpfit")
  mcmclist_draws(object)
  dots = list(...)
  if (length(dots) > 0)
    stop("Unrecognized argument(s): ", and_collapse(names(dots)), call. = FALSE)
  if (!is.null(seed))
    checkmate::assert_int(seed, lower = 1)
  if (!is.null(draws) && !is.null(ndraws))
    stop("Use only one of `draws` and `ndraws`.", call. = FALSE)
  if (!is.null(draws))
    ndraws = draws
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)

  varying = resolve_re_formula(re.form, re_formula)
  if (is.null(dpar))
    dpar = "epred"
  if (!is.null(seed))
    set.seed(seed)
  pp_eval(
    object,
    newdata = newdata,
    summary = FALSE,
    type = type,
    probs = FALSE,
    rate = rate,
    prior = FALSE,
    dpar = if (type == "fitted") dpar else NULL,
    varying = varying,
    arma = TRUE,
    ndraws = ndraws,
    draws_format = "matrix",
    scale = scale,
    .garma_replicate = !conditional
  )
}


# Convert rstantools group-effect syntax to mcp's group-effect selector
resolve_re_formula = function(re.form, re_formula) {
  if (!is.null(re_formula) && !is.null(re.form))
    stop("Use only one of `re.form` and `re_formula`.", call. = FALSE)
  formula = if (!is.null(re_formula)) re_formula else re.form
  if (is.null(formula))
    return(TRUE)
  if (length(formula) == 1 && is.na(formula))
    return(FALSE)
  stop("`re.form`/`re_formula` must be NULL or NA for mcpfit objects.", call. = FALSE)
}


#' @export
log_lik = function(object, ...) UseMethod("log_lik")

#' @aliases log_lik log_lik.mcpfit
#' @describeIn execute-mcp-model Pointwise log-likelihood
#' @export
log_lik.mcpfit = function(
  object,
  newdata = NULL,
  summary = FALSE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "matrix",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "log_lik.mcpfit")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "log_lik.mcpfit")
  dots = list(...)
  warn_which_y(dots, "log_lik")
  dots$which_y = NULL
  if (length(dots) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(names(dots)), call. = FALSE)

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "loglik",
    probs = probs,
    rate = rate,
    prior = prior,
    dpar = NULL,
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    draws_format = draws_format,
    scale = "response"
  )
}


#' @aliases residuals residuals.mcpfit
#' @describeIn execute-mcp-model Residual distribution
#' @export
residuals.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "residuals.mcpfit")
  dots = list(...)
  warn_which_y(dots, "residuals")
  dots$which_y = NULL
  if (length(dots) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(names(dots)), call. = FALSE)

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "residuals",
    probs = probs,
    rate = FALSE,
    prior = prior,
    dpar = NULL,
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    draws_format = "tidy"
  )
}
