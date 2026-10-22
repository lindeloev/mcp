# ABOUT: Posterior predictions, fitted values, residuals, and log-likelihood.
# --------------------------------------------------------------------------

# Evaluate the model on data and draws. Underlies fitted(), predict(), residuals(), log_lik(),
# posterior_*(), plot(), and pp_check(). Public arguments are documented in ?`execute-mcp-model`.
# - type: One of "fitted", "predict", "residuals", or "loglik".
# - .max_cells: Approximate number of draw-row combinations to evaluate at a time. Limits memory use.
pp_eval = function(
  object,
  newdata = NULL,
  summary = TRUE,
  type = "fitted",
  width = 0.95,
  rate = FALSE,
  prior = FALSE,
  dpar = "epred",
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = 'response',
  .max_cells = 1e6
) {
  ########################
  # ASSERTS AND RECODING #
  ########################
  fit = object
  checkmate::assert_class(fit, "mcpfit")
  warn_custom_jags_code(fit)
  draws_format = rlang::arg_match0(draws_format, c("tidy", "matrix"))
  dpar = assert_dpar(dpar, fit = fit, type = type)
  assert_typescale(type, scale)
  value_col = get_value_col(type, dpar, scale)
  if (is.null(group_effects))
    stop("`group_effects` cannot be NULL. Use `TRUE` for all group-level effects or `FALSE`/`NA` for none.", call. = FALSE)
  if (length(group_effects) == 1 && is.na(group_effects))
    group_effects = FALSE  # brms/lme4-style NA means none
  checkmate::assert_flag(summary)
  assert_width(width, false_ok = TRUE)
  assert_rate(rate, fit)
  checkmate::assert_flag(prior)
  assert_arma(arma, replicate_ok = type == "predict")
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)

  # What data to use
  using_original_data = is.null(newdata) || identical(data.frame(newdata), fit$data)
  if (using_original_data)
    newdata = fit$data

  data_columns = mcp_columns(fit)
  replicate_garma = identical(arma, "replicate") && is_arma(fit)
  arma = !isFALSE(arma)  # Whether AR/MA terms are included. `replicate_garma` says how.
  evaluating_arma = uses_garma(fit, arma, type, dpar)
  if (evaluating_arma)
    assert_arma_series(newdata, data_columns$series)
  if (type == "loglik")
    assert_loglik_garma_history(fit, newdata, arma)

  conditional_garma = evaluating_arma && !replicate_garma &&
    (type %in% c("predict", "residuals") ||
       (type == "fitted" && dpar %in% c("epred", "mu")))
  if (conditional_garma && !using_original_data &&
      data_columns$response %in% names(newdata) && anyNA(newdata[[data_columns$response]]))
    stop(
      "Conditional GARMA evaluation with missing responses is supported only for ",
      "the original fitted data, where retained posterior imputations are available. ",
      "Use `predict(..., arma = \"replicate\")` to generate fresh replicated response series from ",
      "predictor-only `newdata`.",
      call. = FALSE
    )


  ###############
  # FIX NEWDATA #
  ###############
  # Identify grouping columns to exclude based on the `group_effects` argument
  group_info = unpack_group_effects(fit, pars = group_effects)
  model_tables = get_fit_model_tables(fit)
  group_cols = unique(stats::na.omit(model_tables$group_effects$group_col))
  exclude_group_cols = setdiff(group_cols, c(group_info$cols, if (evaluating_arma) data_columns$series, get_predictor_cols(fit)))

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
  series_required = !is.null(data_columns$series) && (
    evaluating_arma ||
    data_columns$series %in% get_predictor_cols(fit) ||
    data_columns$series %in% group_info$cols
  )
  if (!series_required && !is.null(data_columns$series))
    required_cols = required_cols[required_cols != data_columns$series]
  response_not_required = type %in% c("fitted", "predict") && !conditional_garma
  if (response_not_required) {
    required_cols = required_cols[required_cols != data_columns$response]
  } else if (data_columns$response %notin% colnames(newdata)) {
    stop("`newdata` must contain a response column named '", data_columns$response, "' for when `arma == TRUE` and/or `type == 'residuals'`")
  }
  missing_group_cols = setdiff(intersect(required_cols, group_info$cols), colnames(newdata))
  if (length(missing_group_cols) > 0)
    stop("These grouping columns are missing from `newdata`: ", and_collapse(missing_group_cols), ". Add them, or set `group_effects = FALSE` for population-level ", type, " values.", call. = FALSE)
  assert_data_cols(newdata, required_cols)
  assert_response_data(
    fit$family,
    model_tables$segments,
    newdata,
    response_required = !response_not_required,
    aux_required = aux_used
  )

  # Validate against reserved output namespace
  assert_reserved_output_namespace(colnames(newdata), context = "newdata", extra = value_col)

  # Filter newdata columns and attach unique evaluation row index
  newdata = data.frame(newdata[, colnames(newdata) %notin% exclude_group_cols, drop = FALSE])
  newdata$.row = seq_len(nrow(newdata))  # Evaluation key throughout summaries, matrices, plots, and metrics
  newdata_return = newdata


  ########################
  # GET FITS/PREDICTIONS #
  ########################
  data = get_model_data(fit, newdata)
  draws = mcp_draws(
    fit, population = TRUE, group = group_effects, prior = prior, ndraws = ndraws,
    newdata = if (length(group_info$cols) > 0) newdata, pivot = FALSE
  )
  groups = split(group_info$effects$name, group_info$effects$group_col)
  y = if (!replicate_garma) newdata[[data_columns$response]]  # NULL generates fresh AR/MA series

  # Complete the conditional history using the retained JAGS imputations
  impute_response = conditional_garma && anyNA(y)
  if (impute_response) {
    if (prior)
      stop("Missing GARMA histories require posterior draws. Use `arma = \"replicate\"` for prior prediction.", call. = FALSE)
    if (!all(model_tables$group_effects$name %in% group_info$pars))
      stop(
        "This model has group-level effects, and its retained missing-response ",
        "histories are conditional on all of them. GARMA evaluation with missing ",
        "responses therefore currently requires `group_effects = TRUE`.",
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
  }

  # This is the important step! Evaluate the mcp model for all rows of newdata (rows) and draws (columns).
  # Evaluate blocks of about `.max_cells` combinations to keep memory use flat for large data.
  keep_dpars = summary && type == "predict"  # For exact predictive moments and quantiles
  block_size = max(1, floor(.max_cells / nrow(newdata)))
  blocks = split(seq_len(nrow(draws)), ceiling(seq_len(nrow(draws)) / block_size))
  evaluated = lapply(blocks, function(block) {
    block_y = if (impute_response) get_imputed_response(fit, y, draws$.draw[block]) else y
    pars = draws_to_pars(draws[block, ], newdata, groups)
    result = evaluate_model(
      fit, data, pars, length(block), block_y,
      type = if (type == "residuals") "fitted" else type, rate = rate, dpar = dpar, arma = arma, scale = scale
    )
    if (keep_dpars) result else result["value"]
  })
  combine = function(get) do.call(cbind, lapply(evaluated, get))  # Blocks of draws are columns
  values = combine(function(block) block$value)
  dpars = if (keep_dpars) lapply(stats::setNames(nm = fit$family$dpars), function(name) combine(function(block) block$dpars[[name]]))
  if (type == "residuals")
    values = y - values

  # Missing outcomes are latent in the fitted JAGS model, but they are not observed-data likelihood contributions
  if (type == "loglik") {
    observed = !is.na(y)
    if (!any(observed))
      stop("Log-likelihood evaluation requires at least one observed response.")
    values = values[observed, , drop = FALSE]
    newdata_return = newdata_return[observed, , drop = FALSE]
  }


  ##########
  # RETURN #
  ##########
  # Summaries for each row of newdata
  if (summary) {
    if (type == "predict" && !is.null(fit$family$r$variance)) {
      # Exact predictive moments over the draws via the laws of total expectation and variance
      epred = matrix(fit$family$r$epred(dpars, data$response, rate = rate), nrow(values))
      variance = matrix(fit$family$r$variance(dpars, data$response, rate = rate), nrow(values))
      newdata_return$.mean = rowMeans(epred)
      newdata_return$.sd = sqrt(rowMeans(variance) + rowMeans((epred - rowMeans(epred))^2))
    } else {
      newdata_return$.mean = rowMeans(values)
      newdata_return$.sd = apply(values, 1, stats::sd)
    }

    # Central interval. For predictions, the exact quantiles of the predictive distribution.
    if (!isFALSE(width)) {
      probs = c((1 - width) / 2, (1 + width) / 2)
      quantiles = if (type == "predict" && !is.null(fit$family$r$cdf)) {
        get_mixture_quantiles(dpars, data$response, probs, fit$family, rate = rate)
      } else {
        t(apply(values, 1, stats::quantile, probs = probs, names = FALSE, na.rm = type == "residuals"))
      }
      newdata_return$.lower = quantiles[, 1]
      newdata_return$.upper = quantiles[, 2]
    }
    return(tibble::as_tibble(newdata_return[names(newdata_return) != ".row"]))
  }

  # Draws-by-rows matrix, as in brms and rstantools
  if (draws_format == "matrix") {
    values = t(values)
    colnames(values) = newdata_return$.row
    return(values)
  }

  # Tidy draws like tidybayes: all draws for each row of newdata, with the parameter values used for that row
  pars = draws_to_pars(draws, newdata_return, groups)
  pars = pars[intersect(mcp_pars(fit)$parameter, names(pars))]
  n_rows = nrow(newdata_return)
  tibble::as_tibble(c(
    lapply(newdata_return, rep, each = nrow(draws)),
    lapply(draws[c(".chain", ".iteration", ".draw")], rep, times = n_rows),
    lapply(pars, function(value) if (is.matrix(value)) as.vector(t(value)) else rep(value, times = n_rows)),
    stats::setNames(list(as.vector(t(values))), value_col)
  ))
}



# Name of the per-draw value column, following tidybayes: `.epred`, `.linpred`,
# `.prediction`, etc. Other dpars are named after the dpar, e.g., `.sigma`.
get_value_col = function(type, dpar = "epred", scale = "response") {
  if (type == "fitted" && scale == "linear")
    return(".linpred")
  if (type == "fitted" && dpar != "epred")
    return(paste0(".", dpar))
  switch(type,
    fitted = ".epred",
    predict = ".prediction",
    residuals = ".residual",
    loglik = ".loglik"
  )
}


#' Fitted and predicted values of `mcp` models fits
#'
#' Evaluate the model on data, either summarised (per data-row) or per draw. You
#' can use draws from the prior (`prior = TRUE`), select a distributional
#' parameter with `dpar`, and choose the response or linear-predictor scale with
#' `scale` where applicable. See [Fits and predictions](https://lindeloev.github.io/mcp/articles/predict.html)
#' for worked examples, including forecasts with future change points.
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
#' but with fixed arguments for `fitted`: `rate = FALSE, dpar = 'epred'`.
#'
#' `log_lik()` defaults to an unsummarised draws-by-observation matrix, as used
#' by `loo` and other posterior workflows. Non-default `group_effects` and `arma`
#' settings evaluate conditional or counterfactual log-likelihoods (e.g.,
#' without group-level effects or serial correlation); they cannot be used in
#' `loo()` or `waic()` because estimating information criteria for reduced
#' models requires refitting.
#'
#' Missing responses in the original data remain missing in the response column.
#' `fitted()` returns their expected responses, while `predict()` generates fresh
#' posterior predictive response draws (including at missing rows). In GARMA models,
#' retained JAGS imputations supply the history used for later fitted and predicted rows.
#' [Read more](https://lindeloev.github.io/mcp/articles/missing.html).
#'
#' @param group_effects Group-level effects to include. One of:
#'   * `TRUE` All group-level deviations.
#'   * `FALSE` or `NA` No group-level deviations, i.e., population-level
#'     predictions. Like `re_formula = NA` in `{brms}`.
#'   * `"cp"` or `"predictor"`: All group-level deviations belonging to that part of
#'     the model.
#'   * Character vector of group-level parameter names, e.g., `c("cp_1_id", "x_2_id")`.
#'     `r lifecycle::badge("experimental")` `*` is a wildcard, e.g., `"state*_1_id"` selects all levels of a factor
#'     and `"*_id"` selects all effects by `id`. See names using `ranef(fit)`.
#'   Can be combined, e.g., `c("cp", "sigma_*")`.
#' @param object An `mcpfit` object.
#' @param newdata A `tibble` or a `data.frame` containing predictors in the model. 
#' - If `NULL` (default), the original data is used.
#' - Grouping factors must only contain levels from the fitted data. Use `group_effects = FALSE`
#'   to predict without group-level effects, e.g., for new groups.
#' - For models with `ar()` or `ma()`: `fitted()`, `residuals()`, `log_lik()`,
#'   and `predict()` condition on the response history by default,
#'   so `newdata` must include the response. For `fitted()`, `predict()`, and
#'   `residuals()`, missing response histories are supported only in the original
#'   fitted data, using retained posterior imputations as histories. Predictions
#'   are fresh response draws, including at missing rows. With `arma = "replicate"`,
#'   `predict()` and [`posterior_predict()`][rstantools::posterior_predict] generate
#'   fresh response series recursively, so their `newdata` need only contain
#'   predictors and required response auxiliaries. `log_lik()`
#'   is unavailable when a missing response enters a later observed history.
#' - For models with `y | weights()`: Require the weights column except for `fitted()` and `predict()`.
#' @param summary Logical. Summarise the draws for each row in `newdata` (`TRUE`) or return the draws (`FALSE`)?
#' @param width Width of the central interval (strictly between 0 and 1),
#'   returned in the `.lower` and `.upper` columns. `FALSE` omits the interval.
#'   Only in effect when `summary == TRUE`. To get custom quantiles, summarise the draws from
#'   `summary = FALSE`, e.g., using `dplyr::reframe()` (see examples).
#' @param rate Logical scalar. For binomial models, return counts (`rate = FALSE`, the default for `fitted()` and `predict()`) or
#'   the observed or expected success proportion (`rate = TRUE`). Predictions and
#'   count-scale fitted values require a trials column in `newdata`.
#'   Distributional parameters such as `dpar = "mu"` evaluate the parameter itself (e.g., success probability)
#'   and are unaffected by `rate`. Other families do not support `rate = TRUE`. For rates per unit of exposure in models with
#'   `offset(log(exposure))`, set `exposure = 1` in `newdata`.
#' @param prior Logical. Evaluate prior draws (`TRUE`) instead of posterior draws (`FALSE`, default).
#'   The selected draws must be available; prior-only fits require `prior = TRUE`.
#' @param dpar What distributional parameter to evaluate. E.g.,
#'
#'   * `"epred"` (default): Expected response from the full model (or `NULL` for compatibility with brms etc.).
#'   * `"mu"`: The conditional mean (or success probability per trial for binomial/bernoulli models), on the link or response scale.
#'   * `"sigma"`: The standard deviation of the residuals.
#'   * `"ar1"`, `"ar2"`, `"ma1"`, `"ma2"`, etc. depending on which AR or MA
#'     coefficient you want to evaluate.
#' @param arma How to include AR and MA effects.
#'   * `TRUE`: Condition on the observed response history, i.e., one-step-ahead evaluation of
#'     the GARMA residual recurrence. Requires the response variable in `newdata`.
#'   * `FALSE`: Disregard AR and MA effects. For `family = gaussian()`, `predict()` uses only `sigma` for residuals.
#'   * `"replicate"`: Generate fresh response series recursively instead of conditioning on
#'     the observed responses. Only for predictions: `predict()`, `posterior_predict()`, and
#'     `pp_check()` (where it is the default). Applies equally to prior and posterior draws.
#'
#'   For posterior evaluation of the original data, retained JAGS imputations
#'   supply missing GARMA histories. In models with group-level effects, this
#'   currently requires all such effects (`group_effects = TRUE`).
#' @param ndraws Integer or `NULL`. Number of posterior draws to return/summarise.
#'   `NULL` means "all". More draws trade speed for accuracy.
#' @param draws_format One of "tidy" or "matrix". Controls the output format when `summary == FALSE`.
#' @param scale One of
#'   * `"response"`: return on the response scale, i.e., after applying the inverse link function.
#'   * `"linear"`: return on the linear-predictor (link) scale, where the linear
#'     trends are modeled.
#' @param varying Deprecated. Use `group_effects` instead.
#' @param probs Deprecated. Use `width` instead.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @param samples_format Deprecated. Use `draws_format` instead.
#' @param ... Must be empty. Reserved for future use.
#' @return
#'   * If `summary = TRUE`: A tibble with the columns of `newdata` (or
#'     `fit$data` if `newdata` is `NULL`), followed by summaries for each row:
#'
#'      - `.mean`: The mean of the draws.
#'      - `.sd`: The SD of the draws. With posterior draws (the default), this is the
#'        posterior predictive SD for `predict()` and the posterior SD of the
#'        evaluated quantity otherwise.
#'      - `.lower` and `.upper`: The lower and upper ends of the central interval given in `width`.
#'
#'     For `predict()` with built-in families, these are the exact moments and quantiles of the
#'     posterior predictive distribution (the mixture over posterior draws), so they do not vary with the seed.
#'
#'     With `prior = TRUE`, these are the analogous prior summaries. mcp's own
#'     columns are dot-prefixed to avoid collisions with data columns.
#'
#'   * If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes` `tibble` with all the posterior
#'     draws (`Nd`) evaluated at each row in `newdata` (`Nn`), i.e., with `Nd x Nn` rows, ordered by
#'     `.row` and then `.draw`. If there are
#'     group-level effects, the returned data is expanded with the relevant levels for each row.
#'
#'     The return columns are:
#'
#'      - Predictors from `newdata`, plus its response column when supplied.
#'      - Draw descriptors: ".chain", ".iteration", ".draw" (see the `posterior` and `tidybayes` packages), and `.row`, the row number in the evaluated `newdata`.
#'      - Draw values: one column for each parameter in the model.
#'      - The estimate (matching tidybayes conventions): `.epred` for expected responses, `.linpred` with `scale = "linear"`,
#'        the `.dpar` for other dpars (e.g., `.mu` or `.sigma`), `.prediction`, `.residual`, or `.loglik`.
#'
#'   * If `summary = FALSE` and `draws_format = "matrix"`: An `N_draws` X `nrows(newdata)` matrix with fitted/predicted
#'       values (depending on the function). This format is used by `brms` and it's useful as `yrep` in
#'      `bayesplot::ppc_*` functions.
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
#' head(predict(demo_fit, width = 0.8))  # 80% posterior predictive interval
#' head(predict(demo_fit, prior = TRUE))  # Prior predictive
#' head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#' head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#'
#' # Custom quantiles from the draws, here for fitted().
#' # Note that for predict(), this adds more Monte Carlo error than `width`.
#' fitted(demo_fit, summary = FALSE) |>
#'   dplyr::group_by(.row, time) |>
#'   dplyr::reframe(quantile = c(0.1, 0.5, 0.9), .epred = quantile(.epred, quantile)) |>
#'   head()
#'
#' # Evaluate at novel data
#' novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
#' head(predict(demo_fit, newdata = novel_data))
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
#' @export
predict.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  width = 0.95,
  rate = FALSE,
  prior = FALSE,
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  probs = lifecycle::deprecated(),
  ...
) {
  ndraws = legacy_arg(ndraws, nsamples, missing(ndraws), "predict.mcpfit", "ndraws", "nsamples")
  group_effects = legacy_arg(group_effects, varying, missing(group_effects), "predict.mcpfit", "group_effects", "varying", lifecycle::deprecate_warn)
  width = legacy_width(width, probs, missing(width), "predict.mcpfit")
  draws_format = legacy_arg(draws_format, samples_format, missing(draws_format), "predict.mcpfit", "draws_format", "samples_format")
  legacy_which_y(list(...), "predict")

  if (missing(rate))
    legacy_binomial_rate(object, newdata, "predict")

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "predict",
    width = width,
    rate = rate,
    prior = prior,
    dpar = NULL,
    group_effects = group_effects,
    arma = arma,
    ndraws = ndraws,
    draws_format = draws_format
  )
}


#' @aliases fitted fitted.mcpfit
#' @describeIn execute-mcp-model Expected response
#' @export
fitted.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  width = 0.95,
  rate = FALSE,
  prior = FALSE,
  dpar = "epred",
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = "response",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  probs = lifecycle::deprecated(),
  ...
) {
  ndraws = legacy_arg(ndraws, nsamples, missing(ndraws), "fitted.mcpfit", "ndraws", "nsamples")
  group_effects = legacy_arg(group_effects, varying, missing(group_effects), "fitted.mcpfit", "group_effects", "varying", lifecycle::deprecate_warn)
  width = legacy_width(width, probs, missing(width), "fitted.mcpfit")
  draws_format = legacy_arg(draws_format, samples_format, missing(draws_format), "fitted.mcpfit", "draws_format", "samples_format")
  which_y = legacy_which_y(list(...), "fitted")
  if (!is.null(which_y) && missing(dpar))
    dpar = which_y

  if (missing(rate) && (is.null(dpar) || identical(dpar, "epred")) && scale == "response")
    legacy_binomial_rate(object, newdata, "fitted")

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "fitted",
    width = width,
    rate = rate,
    prior = prior,
    dpar = dpar,
    group_effects = group_effects,
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
#' mcp exports these generics, so they work without `{rstantools}`. If `{rstantools}` or
#' `{brms}` masks them, their generics dispatch to the same methods.
#'
#' @param object An `mcpfit` object.
#' @inheritParams predict.mcpfit
#' @param draws,ndraws Number of posterior draws to return. `draws` follows the
#'   `{rstantools}` convention; `ndraws` is the mcp spelling. Supply at most one.
#' @param re.form,re_formula Group-level effects to include. `NULL` includes all
#'   effects and `NA` excludes them.
#' @param dpar Distributional parameter for `posterior_epred()` and
#'   `posterior_linpred()`; `NULL` uses the expected response.
#' @param transform For `posterior_linpred()`, return the inverse-linked
#'   parameter on the response scale instead of the linear predictor. For
#'   binomial models, this returns the success probability rather than expected
#'   counts.
#' @param seed Optional integer seed for draw selection and posterior prediction.
#' @param ... Must be empty. Reserved for future use.
#' @return A numeric `N_draws` by `nrow(newdata)` matrix.
#' @details For GARMA models, `posterior_predict()` conditions on the observed
#'   response history, just like `predict()`. Use `arma = "replicate"` in either
#'   method to generate fresh response histories recursively. Missing responses
#'   in the original data are filled with retained imputations only to supply
#'   histories. With `arma = TRUE`, `posterior_predict()` draws from the response
#'   distributions whose means `posterior_epred()` returns, including at missing
#'   rows. These methods require posterior draws. For prior prediction, use
#'   `predict()` with `prior = TRUE`; `arma` selects the same behavior.
#'
#'   For binomial models, `posterior_epred()` and `posterior_predict()` (and
#'   corresponding `{tidybayes}` workflows such as `add_epred_draws()`) follow
#'   `{brms}` and `{rstantools}` conventions by returning values on the outcome
#'   count scale (`rate = FALSE`), i.e., expected counts \eqn{E[Y] = n\mu} and
#'   simulated counts in \eqn{\{0, \dots, n\}}, matching `fitted()` and `predict()`.
#'   Use `rate = TRUE` for proportions. To obtain the success probability
#'   parameter \eqn{\mu} on the \eqn{[0, 1]} scale regardless of trial counts, pass
#'   `dpar = "mu"`. In contrast, `posterior_linpred(transform = TRUE)` returns the
#'   inverse-linked parameter directly, evaluating success probabilities \eqn{\mu}
#'   in \eqn{[0, 1]} rather than counts.
#' @seealso [fitted.mcpfit()], [predict.mcpfit()]
#' @export
#' @exportS3Method rstantools::posterior_epred
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
#' @export
#' @exportS3Method rstantools::posterior_predict
posterior_predict.mcpfit = function(
  object,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  seed = NULL,
  rate = FALSE,
  arma = TRUE,
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
    arma = arma,
    ...
  )
}


#' @rdname posterior_epred.mcpfit
#' @export
#' @exportS3Method rstantools::posterior_linpred
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
    rate = object$family$response$has_rate,  # Binomial: transformed linpred is mu, not counts
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
  arma = TRUE,
  ...
) {
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

  group_effects = resolve_re_formula(re.form, re_formula)
  if (is.null(dpar))
    dpar = "epred"
  if (!is.null(seed))
    set.seed(seed)
  pp_eval(
    object,
    newdata = newdata,
    summary = FALSE,
    type = type,
    width = FALSE,
    rate = rate,
    prior = FALSE,
    dpar = if (type == "fitted") dpar else NULL,
    group_effects = group_effects,
    arma = arma,
    ndraws = ndraws,
    draws_format = "matrix",
    scale = scale
  )
}


# If mcp masks the rstantools generics below, forward e.g. brmsfit objects to them
#' @rdname posterior_epred.mcpfit
#' @export
posterior_epred = function(object, ...) {
  if (!inherits(object, "mcpfit") && isNamespaceLoaded("rstantools"))
    return(rstantools::posterior_epred(object, ...))
  UseMethod("posterior_epred")
}

#' @rdname posterior_epred.mcpfit
#' @export
posterior_predict = function(object, ...) {
  if (!inherits(object, "mcpfit") && isNamespaceLoaded("rstantools"))
    return(rstantools::posterior_predict(object, ...))
  UseMethod("posterior_predict")
}

#' @rdname posterior_epred.mcpfit
#' @export
posterior_linpred = function(object, ...) {
  if (!inherits(object, "mcpfit") && isNamespaceLoaded("rstantools"))
    return(rstantools::posterior_linpred(object, ...))
  UseMethod("posterior_linpred")
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


#' @rdname execute-mcp-model
#' @export
log_lik = function(object, ...) {
  # If mcp masks rstantools::log_lik(), forward e.g. brmsfit objects to it
  if (!inherits(object, "mcpfit") && isNamespaceLoaded("rstantools"))
    return(rstantools::log_lik(object, ...))
  UseMethod("log_lik")
}

#' @aliases log_lik log_lik.mcpfit
#' @describeIn execute-mcp-model Pointwise log-likelihood
#' @export
#' @exportS3Method rstantools::log_lik
log_lik.mcpfit = function(
  object,
  newdata = NULL,
  summary = FALSE,
  width = 0.95,
  prior = FALSE,
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "matrix",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  probs = lifecycle::deprecated(),
  ...
) {
  ndraws = legacy_arg(ndraws, nsamples, missing(ndraws), "log_lik.mcpfit", "ndraws", "nsamples")
  group_effects = legacy_arg(group_effects, varying, missing(group_effects), "log_lik.mcpfit", "group_effects", "varying", lifecycle::deprecate_warn)
  width = legacy_width(width, probs, missing(width), "log_lik.mcpfit")
  draws_format = legacy_arg(draws_format, samples_format, missing(draws_format), "log_lik.mcpfit", "draws_format", "samples_format")
  legacy_which_y(list(...), "log_lik")

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "loglik",
    width = width,
    prior = prior,
    dpar = NULL,
    group_effects = group_effects,
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
  width = 0.95,
  prior = FALSE,
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  nsamples = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  probs = lifecycle::deprecated(),
  ...
) {
  ndraws = legacy_arg(ndraws, nsamples, missing(ndraws), "residuals.mcpfit", "ndraws", "nsamples")
  group_effects = legacy_arg(group_effects, varying, missing(group_effects), "residuals.mcpfit", "group_effects", "varying", lifecycle::deprecate_warn)
  width = legacy_width(width, probs, missing(width), "residuals.mcpfit")
  legacy_which_y(list(...), "residuals")

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "residuals",
    width = width,
    rate = FALSE,
    prior = prior,
    dpar = NULL,
    group_effects = group_effects,
    arma = arma,
    ndraws = ndraws,
    draws_format = draws_format
  )
}
