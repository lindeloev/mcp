#' Fit Multiple Linear Segments And Their Change Points
#'
#' Given a model (a list of segment formulas), `mcp` infers the posterior
#' distributions of the parameters of each segment as well as the change points
#' between segments. [See more details and worked examples on the mcp website](https://lindeloev.github.io/mcp/).
#' All segments must regress on the same x-variable. You can run
#' `fit = mcp(model, data, sample=FALSE)` to avoid sampling if you just want to
#' inspect the priors (`fit$prior` and [prior_summary()]), the JAGS code
#' `fit$jags_code`, or the R function to simulate data (`fit$simulate`).
#'
#' mcp models ordered change points. The ordering is imposed through the prior.
#'
#' @aliases mcp
#' @param data Table-like data in long format (data.frame, tibble, data.table, etc.)
#'   with syntactic column names.
#'   Missing values in the response variable are imputed using the posterior predictive. 
#'   \code{\link{fitted.mcpfit}} or \code{\link{predict.mcpfit}} details how to see the imputed values.
#' @param model A list of formulas - one for each segment. The many examples 
#'   on the [mcp website](https://lindeloev.github.io/mcp/). But briefly:
#' 
#'   The first formula has the format `response ~ predictors` while the following formulas have
#'   the format `response ~ cp ~ predictors`. Here, `cp` names the change-point
#'   part of the formula rather than a literal variable. The response and
#'   change-point parts can be omitted (`cp ~ predictor` assumes the same
#'   response; `~ predictor` assumes an intercept-only change point). Terms normally carry 
#'   into later segments until redefined (see details).
#'
#'   The following terms can be modeled:
#'
#'   * *Regular formulas:* e.g., `~ 1 + x`. [Read more](https://lindeloev.github.io/mcp/articles/formulas.html).
#'
#'   * *Extended formulas:* e.g., `~ x:group + I(x^2) + exp(z)`. [Read more](https://lindeloev.github.io/mcp/articles/formulas.html).
#'     R-side bases such as `scale()`, `poly()`, and `splines::ns()` are evaluated
#'     before sampling, and their fitted scaling or basis is reused for `newdata`.
#'
#'   * *Group-level effects:* e.g., `~ 1 + (1 | id)` for a group-level
#'     intercept, or `~ 1 + (factor || id)` for independent intercept and
#'     factor-contrast deviations. [Read more](https://lindeloev.github.io/mcp/articles/varying.html).
#'
#'   * *Gaussian residual standard deviation:* e.g., `~sigma(1)` for a simple
#'     standard-deviation change or `~sigma(1 + x + group)` for more advanced
#'     structures. Explicit `sigma()` formulas model log-SD, while the implicit constant `sigma_1` in a
#'     model without `sigma()` remains on the response scale.
#'     [Read more](https://lindeloev.github.io/mcp/articles/variance.html)
#'
#'   * *Time-series residuals:* use `ar(p)` and `ma(q)` separately or together,
#'     e.g., `~ 1 + ar(1) + ma(1)`. Both accept an optional regression formula
#'     and observation `boundary`. GARMA terms support Gaussian, binomial,
#'     Poisson, and negative-binomial families with their default links. They
#'     define a finite conditional recurrence and do not jointly constrain AR
#'     coefficients to stationarity or MA coefficients to invertibility.
#'     [Read more](https://lindeloev.github.io/mcp/articles/arma.html)
#' 
#'   * *Weights:* `y | weights(w) ~ ...` usespositive Gaussian precision weights, 
#'     giving residual SD `sigma / sqrt(w)`.
#' 
#'   * *Binomial:* use `successes | trials(total) ~ ...` with `family = binomial()`.
#'
#' @param prior Named list. Names are parameter names (`cp_i`, `Intercept_i`, `xvar_i`,
#'  `sigma`) and the values are either
#'
#'  * A distribution in mcp's JAGS-string syntax (e.g.,
#'      `Intercept_1 = "dnorm(0, 1) T(0,)"`) indicating a
#'      conventional prior distribution. Data-calibrated, regularizing defaults
#'      are used where priors are not specified. These are designed for stable
#'      estimation and prediction, but should be justified before hypothesis testing.
#'      `mcp` uses conventional distribution scales rather than JAGS precision:
#'      SD for `dnorm()`, scale for `dt()`, `ddexp()`, and `dlogis()`, and
#'      log-SD for `dlnorm()`. See
#'      details. With multiple change points, the default is a regularizing
#'      Student-t prior centered at `min(x)` and sequentially truncated between
#'      the preceding change point and `max(x)`. User-specified `dunif()` priors
#'      and explicit truncations are used as written, so users must supply the
#'      ordering bounds, e.g., `dunif(cp_1, max(time))` for `cp_2`.
#'  * A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed value.
#'  * A model parameter name (e.g., `Intercept_2 = "Intercept_1"`), indicating that this parameter is shared -
#'      typically between segments. If two group-level deviations are shared this way,
#'      they will need to have the same grouping variable.
#'  * A scaled Dirichlet prior is supported for change points if they are all set to
#'      `cp_i = "dirichlet(alpha)"` with the same positive `alpha` for every
#'      change point. `alpha = 1` is most often used. This prior is symmetric over segment spacings,
#'      unlike the regularizing t-tail default for multiple change points, but it
#'      mixes less efficiently, so you will often need to set `iter` higher.
#'      It is recommended for hypothesis testing and for the estimation of more
#'      than 5 change points. [Read more](https://lindeloev.github.io/mcp/articles/priors.html).
#' @param family A supported family: `gaussian()`, `binomial()`, `bernoulli()`,
#'   `poisson()`, or `negbinomial()`, with a supported link function; e.g.,
#'   `gaussian(link = "log")`.
#' @param par_x String (default: `NULL` which is auto-detect).
#' @param sample One of
#'   * `"post"`: Sample the posterior.
#'   * `"prior"`: Sample only the prior. Plots, summaries, etc. will
#'       use the prior. This is useful for prior predictive checks.
#'   * `"both"`: Sample both prior and posterior. Plots, summaries, etc.
#'       will default to using the posterior. The prior only has effect when doing
#'       Savage-Dickey density ratios in \code{\link{hypothesis}}.
#'   * `"none"` or `FALSE`: Do not sample. Returns an mcpfit
#'       object without sample. This is useful if you only want to check
#'       prior strings (fit$prior), the JAGS model (fit$jags_code), etc.
#' @param cores Deprecated and ignored. Configure parallel processing with a
#'   [future][future::plan] plan instead, for example
#'   `future::plan(future::multisession, workers = 3)`. With the default future
#'   plan, chains are sampled sequentially. The argument remains available for
#'   backwards compatibility.
#' @param chains Positive integer. Number of chains to run.
#' @param iter Positive integer. Number of post-warmup draws from each chain.
#'   The total number of draws is `iter * chains`.
#' @param warmup Positive integer. Number of initial iterations per chain which
#'   are discarded before sampling. Set higher if needed for sampler adaptation; 
#'   use diagnostics to assess convergence.
#' @param adapt Deprecated; use `warmup` instead.
#' @param inits A list if initial values for the parameters. This can be useful
#'   if a model fails to converge. Read more in \code{\link[rjags]{jags.model}}.
#'   Defaults to `NULL`, i.e., no inits.
#' @param jags_code String. Pass JAGS code to `mcp` to use directly. This is useful if
#'   you want to tweak the code in `fit$jags_code` and run it within the `mcp`
#'   framework. R-side simulation and prediction methods continue to use the
#'   mcp-default formulas (with warning), so they may no longer match custom JAGS code.
#' @param seed `NULL` or a positive integer. Seed for the JAGS random-number
#'   generators. When supplied, `inits` must be a single named list shared by all chains.
#' @param diagnostics Named list of diagnostic warning thresholds. Available
#'   elements are `rhat = 1.01`, `ess_bulk = 400`, `ess_tail = 400`,
#'   `ar = 0.10`, and `ma = 0.10`. An empty list uses these defaults; a partial
#'   list overrides only the supplied values. Set an element to `NULL` to disable
#'   that diagnostic, or use `FALSE` to disable all configurable diagnostic
#'   warnings. In [summary.mcpfit()], `NULL` inherits the settings used to fit
#'   the model, while a list or `FALSE` overrides the diagnostic footer.
#' @param quiet Logical. Suppress routine JAGS output and mcp sampling-status
#'   messages? Defaults to `FALSE`.
#' @param series Only affects models with `ar()` or `ma()` terms.
#'  * `NULL` (default): one long series.
#'  * character: data column name identifying independent AR/MA series.
#' @details
#'   * Terms normally carries into the next segment unless otherwise defined.
#'     `~ 0 + x` adds an increment to the current slope and therefore
#'     maintains continuity, whereas `~ 1 + x` resets the segment level and slope.
#'     Its coefficients are absolute segment parameters, so the fitted mean can
#'     be disjoined at the change point.

#'   * Response evaluation: [fitted.mcpfit()] returns the expected response, while
#'     [predict.mcpfit()] returns posterior-predictive response draws. For
#'     binomial models, `rate = TRUE` returns success proportions and
#'     `rate = FALSE` returns counts; count-scale fitted values and predictions
#'     require trial counts in `newdata`. See [execute-mcp-model] for all
#'     response-evaluation methods.
#'
#'   Notes on priors:
#'   * Default population-level `cp_\*` priors are ordered. For user priors,
#'       `mcp` adds truncation (e.g., `T(cp_1, )`) only when the prior has neither
#'       explicit truncation nor an inherently bounded form such as `dunif()` or
#'       `dirichlet()`. After sampling, all population- and group-level change
#'       points are checked for strict ordering; population-level change points
#'       are also checked against the observed x-range. This includes numerically
#'       fixed change points.
#'   * Data-dependent prior values can be written directly, for example
#'       `min(time)`, `max(time)`, `median(response)`, `mad(response)`,
#'       `max(time) - min(time)`, `segment_width(time)`, `n_segments()`, and `n_cp()`.
#'       They are resolved from the model data before JAGS code is generated.
#'       The older constants `MINX`, `MAXX`, `MEANX`, `SDX`, `MINY`, `MAXY`,
#'       `MEANY`, `SDY`, and `N_CP` remain accepted with a deprecation warning.
#'   * Prior strings use conventional scale parameterizations. `mcp` converts
#'       these to the parameterization required by JAGS when generating code:
#'       inverse variance for `dnorm()`, `dt()`, and `dlnorm()`, and inverse
#'       scale for `ddexp()` and `dlogis()`. Use
#'       `prior_summary(fit)` for resolved priors and
#'       `prior_summary(fit, verbose = TRUE)` for their rules and descriptions.
#' @return An \code{\link{mcpfit}} object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \donttest{
#' # Define the segments using formulas. A change point is estimated between each formula.
#' model = list(
#'   response ~ 1,  # Plateau in the first segment (Intercept_1)
#'   ~ 0 + time,    # Joined slope (time_2) at cp_1
#'   ~ 1 + time     # Disjoined slope (Intercept_3, time_3) at cp_2
#' )
#'
#' # Fit it and sample the prior too.
#' # future::plan(future::multisession, workers = 3)  # Uncomment for parallel sampling
#' data = mcp_example_data("demo")  # Simulated data example
#' demo_fit = mcp(model, data = data, sample = "both", seed = 42)
#'
#' # See parameter estimates
#' summary(demo_fit)
#'
#' # Visual inspection of the results
#' plot(demo_fit)  # Visualization of model fit/predictions
#' plot_pars(demo_fit)  # Parameter distributions
#' pp_check(demo_fit)  # Prior/Posterior predictive checks
#'
#' # Test a hypothesis
#' hypothesis(demo_fit, "cp_1 > 10")
#'
#' # Make predictions
#' fitted(demo_fit)
#' predict(demo_fit)
#' predict(demo_fit, newdata = data.frame(time = c(55.545, 80, 132)))
#'
#' # Compare to a one-intercept-only model (no change points) with default prior
#' model_null = list(response ~ 1)
#' fit_null = mcp(model_null, data = data, par_x = "time")  # fit another model here
#' demo_loo = loo(demo_fit)
#' null_loo = loo(fit_null)
#' loo::loo_compare(demo_loo, null_loo)
#'
#' # Inspect the prior. Useful for prior predictive checks.
#' summary(demo_fit, prior = TRUE)
#' plot(demo_fit, prior = TRUE)
#'
#' # Show all priors. Default priors are added where you don't provide any
#' print(demo_fit$prior)
#'
#' # Set priors and re-run
#' prior = list(
#'   Intercept_1 = 15,
#'   time_2 = "dt(0, 2, 1) T(0, )",  # t-dist slope. Truncated to positive.
#'   cp_2 = "dunif(cp_1, 80)",    # change point to segment 2 > cp_1 and < 80.
#'   Intercept_3 = "Intercept_1"           # Shared intercept between segment 1 and 3
#' )
#'
#' fit3 = mcp(model, data = data, prior = prior)
#'
#' # Show the JAGS model
#' demo_fit$jags_code
#' }
mcp = function(model,
               data,
               prior = list(),
               family = gaussian(),
               par_x = NULL,
               sample = "post",
               cores = NULL,
               chains = 3,
               iter = 3000,
               warmup = 1500,
               adapt = lifecycle::deprecated(),
               inits = NULL,
               jags_code = NULL,
               seed = NULL,
               diagnostics = list(),
               quiet = FALSE,
               series = NULL) {
  custom_jags_code = !is.null(jags_code)

  matched_call = match.call()

  ################
  # CHECK INPUTS #
  ################
  # Check model
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")
  assert_rel(model)
  assert_no_offsets(model)

  # Check data and data-model correspondence
  if (missing(data) || is.null(data) || !is.data.frame(data))
    stop("`data` is required in mcp() since mcp v0.4.0. Passing data = NULL or omitting data is no longer supported.", call. = FALSE)

  checkmate::assert_data_frame(data)
  non_syntactic = names(data)[make.names(names(data)) != names(data)]
  if (length(non_syntactic) > 0)
    stop(
      "`data` has non-syntactic column name(s): ",
      paste0("`", non_syntactic, "`", collapse = ", "), ". Rename them before fitting.",
      call. = FALSE
    )
  data = data.frame(data)
  assert_arma_series(data, series)

  checkmate::assert_string(par_x, null.ok = TRUE)
  par_x = get_par_x(model, data, par_x)
  rhs_vars = get_rhs_vars(model)
  model_vars = unique(c(get_model_vars(model), par_x, series))
  assert_model_data(data, par_x, rhs_vars)
  assert_data_cols(data, cols = model_vars, fail_funcs = c(is.infinite))
  data = data[, model_vars]  # Remove unused data

  # Plain named lists are the v0.4 JAGS-string prior format. Future classed
  # sampler-agnostic prior objects can dispatch before this legacy path.
  validate_prior_v1(prior)

  # Transform family to mcpfamily
  if (!is.family(family) && !is.mcpfamily(family))
    stop("`family` is not a valid family or mcpfamily. Should be gaussian(), binomial(), mcpfamily(guassian(link = 'log')), etc.")

  if (is.mcpfamily(family) == FALSE)
    family = mcpfamily(family)

  # More checking...
  checkmate::assert(
    checkmate::check_choice(sample, c("post", "prior", "both", "none")),
    checkmate::check_false(sample),
    .var.name = "sample"
  )
  if (!is.null(cores)) {
    if (!identical(cores, "all"))
      checkmate::assert_int(cores, lower = 1)

    cores_details = paste0(
      "`cores` is ignored. Parallel processing is now controlled by the active ",
      "future plan. For example, call ",
      "`future::plan(future::multisession, workers = 3)` before `mcp()`, and ",
      "`future::plan(future::sequential)` when the workers are no longer needed."
    )
    if (identical(cores, "all") || (is.numeric(cores) && cores > 1))
      cores_details = paste0(
        "Setting `cores` above one no longer enables parallel processing. ",
        cores_details
      )

    lifecycle::deprecate_warn(
      when = "0.4.0",
      what = "mcp(cores)",
      details = cores_details
    )
  }

  if (lifecycle::is_present(adapt)) {
    if (!missing(warmup))
      stop("Supply only one of `warmup` and deprecated `adapt`.", call. = FALSE)

    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = "mcp(adapt)",
      with = "mcp(warmup)"
    )
    warmup = adapt
  }

  checkmate::assert_int(chains, lower = 1)
  checkmate::assert_int(iter, lower = 1)
  checkmate::assert_int(warmup, lower = 1)
  checkmate::assert_list(inits, null.ok = TRUE)
  checkmate::assert_int(seed, lower = 1, null.ok = TRUE)
  diagnostics = resolve_diagnostics(diagnostics)
  checkmate::assert_flag(quiet)

  # jags_code
  if(!is.null(jags_code))
    if (!is.character(jags_code) || !stringr::str_detect(gsub(" ", "", jags_code), "model\\{"))
      stop("`jags_code` must be NULL or a string with a JAGS model, including 'model {...}'.")


  ##################
  # MODEL BUILDING #
  ##################
  # Build model metadata.
  segment_tables = get_segment_tables(model, data, family, par_x)
  segments = segment_tables$segments
  cps = segment_tables$cps
  assert_model_data(data, par_x, group_cols = stats::na.omit(cps$group_col))
  predictor_tables = get_predictor_tables(model, data, family, par_x)
  predictors = predictor_tables$predictors
  family = resolve_dpar_specs(family, predictors, model)
  group_effects = get_group_effects(cps, predictor_tables$group_effects)

  # Make prior
  prior = get_prior(segments, cps, predictors, group_effects, family, prior, data)
  prior_table = attr(prior, "prior_table")
  prior_context = attr(prior, "prior_context")
  attr(prior, "prior_table") = NULL
  attr(prior, "prior_context") = NULL
  assert_fixed_sigma(prior_table, predictors, family)

  # Assemble model metadata used by fitted-model methods
  parameters = get_pars_table(predictors, cps, group_effects, family)
  data_columns = c(
    list(par_x = par_x, response = unique(segments$y), series = series),
    lapply(get_family_aux_columns(family, segments), function(column) {
      if (is.na(column)) NULL else column
    })
  )
  model_tables = list(
    data_columns = data_columns,
    segments = segments,
    cps = cps,
    predictors = predictors,
    group_effects = group_effects,
    parameters = parameters,
    design_specs = predictor_tables$design_specs
  )

  # Validate AR/MA configuration
  has_arma = any(predictors$dpar %in% c("ar", "ma"))
  if (has_arma) {
    if (is.null(family$garma))
      stop(
        "family = ", family$family, "(link = \"", family$link,
        "\") does not define the GARMA behavior required by ar() or ma()."
      )

    x_by_series = split(data[[par_x]], if (is.null(series)) 1 else data[[series]])
    x_unordered = any(vapply(
      x_by_series,
      function(x) is.unsorted(x) && is.unsorted(rev(x)),
      logical(1)
    ))
    if (x_unordered)
      message("'", par_x, "' is unordered. Please note that ar() and ma() apply in data-frame row order, not the values of '", par_x, "'.")
  } else if (!is.null(series)) {
    stop("`series` is only used by models containing ar() or ma().")
  }

  # Make formulas
  formula_jags = get_formula_jags(segments, predictors, group_effects, par_x, family)
  formula_r = get_formula_r(formula_jags, predictors, group_effects, cps, par_x)

  # Make jags code if it is not provided by the user
  if (is.null(jags_code)) {
    ar_order = get_arma_order(predictors, "ar")
    ma_order = get_arma_order(predictors, "ma")
    jags_code = get_jags_code(
      prior, segments, group_effects, formula_jags, ar_order, ma_order, family, par_x,
      prior_table, prior_context, series = !is.null(series)
    )
  }


  ##########
  # SAMPLE #
  ##########
  jags_data = get_jags_data(
    data, family, segments, predictors, group_effects, jags_code, series
  )

  # Monitor model parameters and, for generated JAGS code, latent responses
  all_pars = names(prior)
  missing_response_rows = which(is.na(data[[segments$y[1]]]))
  imputed_response_nodes = if (custom_jags_code || length(missing_response_rows) == 0) character() else
    paste0(segments$y[1], "[", missing_response_rows, "]")

  # Sample posterior
  if (sample %in% c("post", "both")) {
    mcmc_post = run_jags(
      data = data,
      jags_code = jags_code,
      jags_data = jags_data,
      pars = c(all_pars, imputed_response_nodes),
      sample = "post",
      n.chains = chains,
      n.iter = iter,
      n.adapt = warmup,
      inits = inits,
      seed = seed,
      quiet = quiet
    ) %>%
      recover_levels(data, group_effects)

    class(mcmc_post) = "mcmc.list"
    if (length(imputed_response_nodes) > 0) {
      mcmc_imputed = mcmc_post[, imputed_response_nodes, drop = FALSE]
      retained_parameter_nodes = setdiff(colnames(mcmc_post[[1]]), imputed_response_nodes)
      mcmc_post = mcmc_post[, retained_parameter_nodes, drop = FALSE]
    } else {
      mcmc_imputed = NULL
    }
    assert_ordered_cp_draws(mcmc_post, cps, data[[par_x]])
    warn_nonconvergence(mcmc_post, diagnostics)
  } else {
    mcmc_post = NULL
    mcmc_imputed = NULL
  }

  # Sample prior
  if (sample %in% c("prior", "both")) {
    # Set response = NA if we only sample prior
    jags_data_prior = jags_data
    jags_data_prior[[segments$y[1]]] = rep(NA, nrow(data))

    mcmc_prior = run_jags(
      data = data,
      jags_code = jags_code,
      jags_data = jags_data_prior,
      pars = all_pars,
      sample = "prior",
      n.chains = chains,
      n.iter = iter,
      n.adapt = warmup,
      inits = inits,
      seed = seed,
      quiet = quiet
    ) %>%
      recover_levels(data, group_effects)

    class(mcmc_prior) = "mcmc.list"
    assert_ordered_cp_draws(mcmc_prior, cps, data[[par_x]])
  } else {
    mcmc_prior = NULL
  }


  ##########
  # RETURN #
  ##########
  # Return normalized formulas without discarding user-defined environments.
  model = Map(stats::as.formula, segments$form, env = segments$form_env)
  class(model) = c("mcplist", "list")
  class(prior) = c("mcplist", "list")
  class(jags_code) = c("mcptext", "character")  # for nicer printing

  # Make mrpfit object
  mcpfit = list(
    # By user (same order as mcp argument)
    model = model,
    data = data,
    prior = prior,
    family = family,
    call = matched_call,

    # Results
    mcmc_post = mcmc_post,
    mcmc_prior = mcmc_prior,

    # Extracted model
    jags_code = jags_code,
    simulate = get_fitsimulate(cps, predictors, group_effects),

    # Pass info to *.mcpfit() functions.
    # Not meant to be used by the end user.
    .internal = list(
      model_tables = model_tables,
      formula_jags = formula_jags,
      formula_r = formula_r,
      prior_table = prior_table,
      prior_context = prior_context,
      prior_format = "jags_string_v1",
      diagnostics = diagnostics,
      custom_jags_code = custom_jags_code,
      imputed_response = mcmc_imputed,
      imputed_response_rows = if (is.null(mcmc_imputed)) integer() else missing_response_rows
    )
  )
  class(mcpfit) = "mcpfit"

  if (!is.null(mcmc_post) && has_arma)
    warn_arma_fit(mcpfit, diagnostics = diagnostics)

  # Return it
  mcpfit
}
