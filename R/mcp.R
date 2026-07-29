#' Fit Multiple Linear Segments And Their Change Points
#'
#' Given a model (a list of segment formulas), `mcp` infers the posterior
#' distributions of the parameters of each segment as well as the change points
#' between segments. [See more details and worked examples on the mcp website](https://lindeloev.github.io/mcp/).
#' All segments must regress on the same x-variable. Change
#' points are forced to be ordered using truncation of the priors. You can run
#' `fit = mcp(model, data, sample=FALSE)` to avoid sampling if you just want to
#' inspect the priors (`fit$prior` and [prior_summary()]), the JAGS code
#' `fit$jags_code`, or the R function to simulate data (`fit$simulate`).
#'
#' @aliases mcp
#' @param data Table-like data in long format (data.frame, tibble, data.table, etc.)
#' @param model A list of formulas - one for each segment. The first formula
#'   has the format `response ~ predictors` while the following formulas
#'   have the format `response ~ changepoint ~ predictors`. The response
#'   and change points can be omitted (`changepoint ~ predictors` assumes same
#'   response. `~ predictors` assumes an intercept-only change point). The
#'   following can be modeled:
#'
#'   * *Regular formulas:* e.g., `~ 1 + x`). [Read more](https://lindeloev.github.io/mcp/articles/formulas.html).
#'
#'   * *Extended formulas:*, e.g., `~ I(x^2) + exp(z)`. [Read more](https://lindeloev.github.io/mcp/articles/formulas.html).
#'
#'   * *Variance:* e.g., `~sigma(1)` for a simple variance change or
#'     `~sigma(1 + I(x^2))`) for more advanced variance structures. Explicit
#'     sigma formulas model log-SD, while the implicit constant `sigma_1` in a
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
#' @param prior Named list. Names are parameter names (`cp_i`, `Intercept_i`, `xvar_i`,
#'  `sigma``) and the values are either
#'
#'  * A JAGS distribution (e.g., `Intercept_1 = "dnorm(0, 1) T(0,)"`) indicating a
#'      conventional prior distribution. Uninformative priors based on data
#'      properties are used where priors are not specified. This ensures good
#'      parameter estimations, but it is a questionable for hypothesis testing.
#'      `mcp` uses SD (not precision) for dnorm, dt, dlogis, etc. See
#'      details. Change points are forced to be ordered through the priors using
#'      truncation, except for uniform priors where the lower bound should be
#'      greater than the previous change point, `dunif(cp_1, max(time))`.
#'  * A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed value.
#'  * A model parameter name (e.g., `Intercept_2 = "Intercept_1"`), indicating that this parameter is shared -
#'      typically between segments. If two varying effects are shared this way,
#'      they will need to have the same grouping variable.
#'  * A scaled Dirichlet prior is supported for change points if they are all set to
#'      `cp_i = "dirichlet(N)` where `N` is the alpha for this change point and
#'      `N = 1` is most often used. This prior is less informative about the
#'      location of the change points than the default uniform prior, but it
#'      samples less efficiently, so you will often need to set `iter` higher.
#'      It is recommended for hypothesis testing and for the estimation of more
#'      than 5 change points. [Read more](https://lindeloev.github.io/mcp/articles/priors.html).
#' @param family One of `gaussian()`, `binomial()`, `bernoulli()`, `poisson()`,
#'   or `negbinomial()`.
#'   with a supported link function, e.g., `gaussian(link = "log")`.
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
#' @param adapt Positive integer. Also sometimes called "burnin", this is the
#'   number of samples used to reach convergence. Set lower for greater speed.
#'   Set higher if the chains haven't converged yet or look at [tips, tricks, and debugging](https://lindeloev.github.io/mcp/articles/tips.html).
#' @param inits A list if initial values for the parameters. This can be useful
#'   if a model fails to converge. Read more in \code{\link[rjags]{jags.model}}.
#'   Defaults to `NULL`, i.e., no inits.
#' @param jags_code String. Pass JAGS code to `mcp` to use directly. This is useful if
#'   you want to tweak the code in `fit$jags_code` and run it within the `mcp`
#'   framework.
#' @details Notes on priors:
#'   * Order restriction is automatically applied to cp_\* parameters using
#'       truncation (e.g., `T(cp_1, )`) so that they are in the correct order on the
#'       x-axis UNLESS you do it yourself. The one exception is for dunif
#'       distributions where you have to do it as above.
#'   * Data-dependent prior values can be written directly, for example
#'       `min(time)`, `max(time)`, `median(response)`, `mad(response)`,
#'       `max(time) - min(time)`, `segment_width(time)`, `n_segments()`, and `n_cp()`.
#'       They are resolved from the model data before JAGS code is generated.
#'       The older constants `MINX`, `MAXX`, `MEANX`, `SDX`, `MINY`, `MAXY`,
#'       `MEANY`, `SDY`, and `N_CP` remain accepted with a deprecation warning.
#'   * Use SD when you specify priors for dt, dlogis, etc. JAGS uses precision
#'       but `mcp` converts to precision under the hood via the sd_to_prec()
#'       function. So you will see SDs in `fit$prior` but precision ($1/SD^2)
#'       in `fit$jags_code`. Use `prior_summary(fit)` for resolved priors and
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
#' demo_fit = mcp(model, data = data, sample = "both")
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
#' demo_fit$loo = loo(demo_fit)
#' fit_null$loo = loo(fit_null)
#' loo::loo_compare(demo_fit$loo, fit_null$loo)
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
               adapt = 1500,
               inits = NULL,
               jags_code = NULL) {

  ################
  # CHECK INPUTS #
  ################
  # Check model
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")
  assert_rel(model)

  # Check data and data-model correspondence
  checkmate::assert_data_frame(data)
  data = data.frame(data)

  checkmate::assert_string(par_x, null.ok = TRUE)
  par_x = get_par_x(model, data, par_x)
  rhs_vars = get_rhs_vars(model)
  assert_data_cols(data, cols = rhs_vars, fail_funcs = c(is.na, is.nan))

  model_vars = unique(c(get_model_vars(model), par_x))
  assert_data_cols(data, cols = model_vars, fail_funcs = c(is.infinite))
  data = data[, model_vars]  # Remove unused data

  # Check prior
  checkmate::assert_list(prior)

  which_duplicated = duplicated(names(prior))
  if (any(which_duplicated))
    stop("`prior` has duplicated entries for the same parameter: ", and_collapse(names(prior)[which_duplicated]))

  # Transform family to mcpfamily
  if (is.family(family) == FALSE & is.mcpfamily(family) == FALSE)
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

  checkmate::assert_int(chains, lower = 1)
  checkmate::assert_list(inits, null.ok = TRUE)

  # jags_code
  if(!is.null(jags_code))
    if (!is.character(jags_code) || !stringr::str_detect(gsub(" ", "", jags_code), "model\\{"))
      stop("`jags_code` must be NULL or a string with a JAGS model, including 'model {...}'.")


  ##################
  # MODEL BUILDING #
  ##################
  # Make an abstract table representing the segments and their relations.
  par_x = get_par_x(model, data, par_x)
  segment_table = get_segment_table(model, data, family, par_x)  #"ST" for "segment table", "CP" for "change points"
  ST = segment_table$ST
  CP = segment_table$CP
  rhs_table = get_rhs_table(model, data, family, par_x)
  family = resolve_dpar_specs(family, rhs_table, model)

  # Make prior
  prior = get_prior(ST, CP, rhs_table, family, prior, data)
  prior_table = attr(prior, "prior_table")
  prior_context = attr(prior, "prior_context")
  attr(prior, "prior_table") = NULL
  attr(prior, "prior_context") = NULL

  # Make lists of parameters
  all_pars = names(prior)  # There is a prior for every parameter
  family_dpars = family$dpar_specs$dpar
  pars_table = get_pars_table(rhs_table, CP, family)
  pars = list(
    x = par_x,
    y = unique(ST$y),
    cp = paste0("cp_", 1:nrow(ST))[seq_len(nrow(ST)-1)],  # N_cp = N_segments - 1
    fixed = pars_table$name[pars_table$dpar == "mu"],
    population = c(),
    varying = logical0_to_null(c(stats::na.omit(ST$cp_group))),
    sigma = pars_table$name[pars_table$dpar %in% setdiff(family_dpars, "mu")],
    arma = pars_table$name[pars_table$dpar %notin% c("cp", family_dpars)],
    trials = logical0_to_null(stats::na.omit(unique(ST$trials))),
    weights = logical0_to_null(stats::na.omit(unique(ST$weights)))
  )
  cp_population = pars_table$name[pars_table$dpar == "cp" & pars_table$name %notin% pars$varying]
  pars$population = pars_table$name[pars_table$name %in% c(cp_population, pars$fixed, pars$sigma, pars$arma)]

  # Check parameters
  # ARMA models
  if (length(pars$arma) > 0) {
    if (is.null(family$garma))
      stop(
        "family = ", family$family, "(link = \"", family$link,
        "\") does not define the GARMA behavior required by ar() or ma()."
      )

    if (is.unsorted(data[, par_x]) && is.unsorted(rev(data[, par_x])))
      message("'", par_x, "' is unordered. Please note that ar() and ma() apply in data-frame row order, not the values of '", par_x, "'.")
  }

  # Make formulas
  formula_jags = get_formula_jags(ST, rhs_table, par_x, family)
  formula_r = get_formula_r(formula_jags, rhs_table, pars)

  ar_order = get_arma_order(rhs_table, "ar")
  ma_order = get_arma_order(rhs_table, "ma")

  # Make jags code if it is not provided by the user
  if (is.null(jags_code)) {
    jags_code = get_jags_code(
      prior, ST, formula_jags, ar_order, ma_order, family, par_x,
      prior_table, prior_context
    )
  }
  if (sample %in% c("post", "both"))
    warn_high_order_arma(ar_order, ma_order, "fit")


  ##########
  # SAMPLE #
  ##########
  jags_data = get_jags_data(data, family, ST, rhs_table, jags_code)

  # Sample posterior
  if (sample %in% c("post", "both")) {
    mcmc_post = run_jags(
      data = data,
      jags_code = jags_code,
      jags_data = jags_data,
      pars = all_pars,  # Monitor log-likelihood for loo/waic
      sample = "post",
      n.chains = chains,
      n.iter = iter,
      n.adapt = adapt,
      inits = inits
    ) %>%
      recover_levels(data, ST)

    class(mcmc_post) = "mcmc.list"
    warn_nonconvergence(mcmc_post)
  } else {
    mcmc_post = NULL
  }

  # Sample prior
  if (sample %in% c("prior", "both")) {
    # Set response = NA if we only sample prior
    jags_data_prior = jags_data
    jags_data_prior[[ST$y[1]]] = rep(NA, nrow(data))

    mcmc_prior = run_jags(
      data = data,
      jags_code = jags_code,
      jags_data = jags_data_prior,
      pars = all_pars,
      sample = "prior",
      n.chains = chains,
      n.iter = iter,
      n.adapt = adapt,
      inits = inits
    ) %>%
      recover_levels(data, ST)

    class(mcmc_prior) = "mcmc.list"
  } else {
    mcmc_prior = NULL
  }


  ##########
  # RETURN #
  ##########
  model = lapply(ST$form, stats::as.formula, env = globalenv())
  class(model) = c("mcplist", "list")
  class(prior) = c("mcplist", "list")
  class(pars) = c("mcplist", "list")  # for nicer printing
  class(jags_code) = c("mcptext", "character")  # for nicer printing

  # Make mrpfit object
  mcpfit = list(
    # By user (same order as mcp argument)
    model = model,
    data = data,
    prior = prior,
    family = family,

    # Results
    mcmc_post = mcmc_post,
    mcmc_prior = mcmc_prior,
    loglik = NULL,
    loo = NULL,
    waic = NULL,

    # Extracted model
    pars = pars,
    jags_code = jags_code,
    simulate = get_fitsimulate(pars),

    # Pass info to *.mcpfit() functions.
    # Not meant to be used by the end user.
    .internal = list(
      ST = ST,
      CP = CP,
      rhs_table = rhs_table,
      pars_table = pars_table,
      formula_jags = formula_jags,
      formula_r = formula_r,
      prior_table = prior_table,
      prior_context = prior_context,
      mcp_version = utils::packageVersion("mcp")  # For helpful messages about backwards compatibility
    )
  )
  class(mcpfit) = "mcpfit"

  # Return it
  mcpfit
}
