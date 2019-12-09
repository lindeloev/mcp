#' Fit Multiple Linear Segments And Their Change Points
#'
#' Given a list of linear segments, `mcp` infers the posterior
#' distributions of the parameters of each segment as well as the change points
#' between segments. [See more details and worked examples on the mcp website](https://lindeloev.github.io/mcp/).
#' All segments must regress on the same x-variable. Change
#' points are forced to be ordered using truncation of the priors. You can run
#' `fit = mcp(segments, sample=FALSE)` to avoid sampling and the need for
#' data if you just want to get the priors (`fit$prior`), the JAGS code
#' `fit$jags_code`, or the R function to simulate data (`fit$simulate`).
#'
#' @aliases mcp
#' @param data Data.frame or tibble in long format.
#' @param segments A list of formulas - one for each segment. The first formula
#'   has the format `response ~ predictors` while the following formulas
#'   have the format `response ~ changepoint ~ predictors`. The response
#'   and change points can be omitted (`changepoint ~ predictors` assumes same
#'   response. `~ predictors` assumes an intercept-only change point). Segments
#'   can model
#'
#'   * *Regular formulas:* e.g., `~ 1 + x`). [Read more](https://lindeloev.github.io/mcp/articles/formulas.html).
#'
#'   * *Extended formulas:*, e.g., `~ I(x^2) + exp(x) + sin(x)`. [Read more](https://lindeloev.github.io/mcp/articles/formulas.html).
#'
#'   * *Variance:* e.g., `~sigma(1)` for a simple variance change or
#'     `~sigma(rel(1) + I(x^2))`) for more advanced variance structures. [Read more](https://lindeloev.github.io/mcp/articles/variance.html)
#'
#'   * *Autoregression:* e.g., `~ar(1)` for a simple onset/change in AR(1) or
#'     `ar(2, 0 + x`) for an AR(2) increasing by `x`. [Read more](https://lindeloev.github.io/mcp/articles/arma.html)
#'
#' @param prior Named list. Names are parameter names (`cp_i`, `int_i`, `xvar_i`,
#'  `sigma``) and the values are either
#'
#'  * A JAGS distribution (e.g., `int_1 = "dnorm(0, 1) T(0,)"`) indicating a
#'      conventional prior distribution. Uninformative priors based on data
#'      properties are used where priors are not specified. This ensures good
#'      parameter estimations, but it is a questionable for hypothesis testing.
#'      `mcp` uses SD (not precision) for dnorm, dt, dlogis, etc. See
#'      details. Change points are forced to be ordered through the priors using
#'      truncation, except for uniform priors where the lower bound should be
#'      greater than the previous change point, `dunif(cp_1, MAXX)`.
#'  * A numerical value (e.g., `int_1 = -2.1`) indicating a fixed value.
#'  * A model parameter name (e.g., `int_2 = "int_1"`), indicating that this parameter is shared -
#'      typically between segments. If two varying effects are shared this way,
#'      they will need to have the same grouping variable.
#' @param family One of `gaussian()`, `binomial()`, `bernoulli()`, or `poission()`.
#'   Only default link functions are currently supported.
#' @param par_x String (default: NULL). Only relevant if no segments contains
#'   slope (no hint at what x is). Set this, e.g., par_x = "time".
#' @param sample One of
#'   * `"post"` (default): Sample the posterior.
#'   * `"prior"`: Sample only the prior. Plots, summaries, etc. will
#'       use the prior. This is useful for prior predictive checks.
#'   * `"both"`: Sample both prior and posterior. Plots, summaries, etc.
#'       will default to using the posterior. The prior only has effect when doing
#'       Savage-Dickey density ratios in \code{\link{hypothesis}}.
#'   * `"none"` or `FALSE`: Do not sample. Returns an mcpfit
#'       object without sample. This is useful if you only want to check
#'       prior strings (fit$prior), the JAGS model (fit$jags_code), etc.
#' @param cores Positive integer or "all". Number of cores.
#'
#'   * `1`: serial sampling
#'   * `>1`: parallel sampling on this number of cores. Ideally set `chains`
#'     to the same value. Note: `cores > 1` takes a few extra seconds the first
#'     time it's called but subsequent calls will start sampling immediately.
#'   * `"all"`: use all cores but one and sets `chains` to the same value. This is
#'     a convenient way to maximally use your computer's power.
#' @param chains Positive integer. Number of chains to run.
#' @param iter Positive integer. Number of post-warmup samples to draw.
#' @param adapt Positive integer. Also sometimes called "burnin", this is the
#'   number of samples used to reach convergence. Set lower for greater speed.
#'   Set higher if the chains haven't converged yet or look at [tips, tricks, and debugging](https://lindeloev.github.io/mcp/articles/tips.html).
#' @param inits A list if initial values for the parameters. This can be useful
#'   if a model fails to converge. Read more in \code{\link[rjags]{jags.model}}.
#'   Defaults to `NULL`, i.e., no inits.
#' @param jags_code Pass JAGS code to `mcp` to use directly. This is useful if
#'   you want to tweak the code in `fit$jags_code` and run it within the `mcp`
#'   framework.
#' @details Notes on priors:
#'   * Order restriction is automatically applied to cp_\* parameters using
#'       truncation (e.g., `T(cp_1, )`) so that they are in the correct order on the
#'       x-axis UNLESS you do it yourself. The one exception is for dunif
#'       distributions where you have to do it as above.
#'   * In addition to the model parameters, `MINX` (minimum x-value), `MAXX`
#'       (maximum x-value), `SDX` (etc...), `MINY`, `MAXY`, and `SDY`
#'       are also available when you set priors. They are used to set uninformative
#'       default priors.
#'   * Use SD when you specify priors for dt, dlogis, etc. JAGS uses precision
#'       but `mcp` converts to precision under the hood via the sd_to_prec()
#'       function. So you will see SDs in `fit$prior` but precision ($1/SD^2)
#'       in `fit$jags_code`
#' @return An \code{\link{mcpfit}} object.
#' @seealso \link{get_segment_table}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @importFrom stats gaussian binomial
#' @export
#' @examples
#' \donttest{
#' # Define the segments using formulas. A change point is estimated between each formula.
#' segments = list(
#'   response ~ 1,  # Plateau in the first segment (int_1)
#'   ~ 0 + time,    # Joined slope (time_2) at cp_1
#'   ~ 1 + time     # Disjoined slope (int_3, time_3) at cp_2
#' )
#'
#' # Fit it. The `ex_demo` dataset is included in mcp. Sample the prior too.
#' # options(mc.cores = 3)  # Run in parallel to speed up. Disable if this fails.
#' ex_fit = mcp(segments, data = ex_demo, sample = "both")
#' }
#'
#' # See parameter estimates
#' summary(ex_fit)
#'
#' # Visual inspection of the results
#' plot(ex_fit)
#' plot_pars(ex_fit)
#'
#' # Test a hypothesis
#' hypothesis(ex_fit, "cp_1 > 10")
#'
#' \donttest{
#' # Compare to a one-intercept-only model (no change points) with default prior
#' segments_null = list(response ~ 1)
#' fit_null = mcp(segments_null, data = ex_demo, par_x = "time")  # fit another model here
#' ex_fit$loo = loo(ex_fit)
#' fit_null$loo = loo(fit_null)
#' loo::loo_compare(ex_fit$loo, fit_null$loo)
#' }
#'
#' # Inspect the prior. Useful for prior predictive checks.
#' summary(ex_fit, prior = TRUE)
#' plot(ex_fit, prior = TRUE)
#'
#' # Show all priors. Default priors are added where you don't provide any
#' print(ex_fit$prior)
#'
#' # Set priors and re-run
#' prior = list(
#'   int_1 = "dt(10, 5, 1) T(0, )", # t-dist intercept. Truncated to positive.
#'   year_1 = "dnorm(0, 3)",      # slope of segment 1. Mean = 0, SD = 5.
#'   cp_2 = "dunif(cp_1, 80)",    # change point to segment 2 > cp_1 and < 80.
#'   year_2 = "year_1",           # Shared slope between segment 2 and 1
#'   int_3 = 15                   # Fixed intercept of segment 3
#' )
#'
#' \donttest{
#' fit3 = mcp(segments, data = ex_demo, prior = prior)
#' }
#'
#' # Show the JAGS model
#' cat(ex_fit$jags_code)

mcp = function(segments,
               data = NULL,
               prior = list(),
               family = gaussian(),
               par_x = NULL,
               sample = "post",
               cores = 1,
               chains = 3,
               iter = 3000,
               adapt = 1000,
               inits = NULL,
               jags_code = NULL) {

  ################
  # CHECK INPUTS #
  ################

  # Check data
  if (is.null(data) & sample == TRUE)
    stop("Cannot sample without data.")

  if (!is.null(data)) {
    if (!is.data.frame(data) & !tibble::is_tibble(data))
      stop("`data` must be a data.frame or a tibble.")

    data = data.frame(data)  # Force into data frame
  }

  # Check segments
  if (!is.list(segments))
    stop("`segments` must be a list")

  if (length(segments) == 0)
    stop("At least one segment is needed")

  for (segment in segments) {
    if (!inherits(segment, "formula"))
      stop("all segments must be formulas.")
  }

  # Check prior
  if (!is.list(prior))
    stop("`prior` must be a named list.")

  which_duplicated = duplicated(names(prior))
  if (any(which_duplicated))
    stop("`prior` has duplicated entries for the same parameter: ", paste0(names(prior)[which_duplicated]), collapse = ", ")

  # Check and recode family
  if (class(family) != "family")
    stop("`family` must be one of gaussian() or binomial()")

  if (!family$family %in% c("gaussian", "binomial", "bernoulli", "poisson"))
    stop("`family` must be one of gaussian(), binomial(), or bernoulli()")

  if (family$family == "gaussian" & !family$link %in% c("identity"))
    stop("'identity' is currently the only supported link function for gaussian().")

  if (family$family %in% c("binomial", "bernoulli") & !family$link %in% c("logit"))
    stop("'logit' is currently the only supported link function for binomial() and bernoulli().")

  if (family$family == "poisson" & !family$link %in% c("log"))
    stop("'log' is currently the only supported link function for poisson().")

  # if (family$family %in% c("bernoulli", "binomial") & !is.null(autocor))
  #   stop("`autocor` is not yet supported for Bernoulli and binomial families.")

  family = get_family(family$family, family$link)  # convert to mcp family


  # Check other stuff
  if (!is.null(par_x) & !is.character(par_x))
    stop("`par_x` must be NULL or a string.")

  # Sampler settings
  if (!sample %in% c("post", "prior", "both") & !is.logical(sample))
    stop("`sample` must be 'post', 'prior', 'both', or 'none'/FALSE")

  if (cores < 1 | !check_integer(cores, "cores"))
    stop("`cores` has to be 1 or greater (parallel sampling).")

  if (chains < 1 | !check_integer(chains, "chains"))
    stop("`chains` has to be 1 or greater.")

  if (cores > chains)
    message("`cores` is greater than `chains`. Not all cores will be used.")

  # Parallel fails on R version 3.6.0 and lower (sometimes at least).
  if (cores > 1 & getRversion() < "3.6.1")
    message("Parallel sampling (`cores` > 1) sometimes err on R versions below 3.6.1. You have ", R.Version()$version.string, ". Consider upgrading if it fails or hangs.")


  ##################
  # MODEL BUILDING #
  ##################
  # Make an abstract table representing the segments and their relations.
  # ("ST" for "segment table").
  ST = get_segment_table(segments, data, family, par_x)

  # Make prior
  prior = get_prior(ST, family, prior)

  # Make lists of parameters
  all_pars = names(prior)  # There is a prior for every parameter
  pars = list(
    x = unique(ST$x),
    y = unique(ST$y),
    trials = logical0_to_null(stats::na.omit(unique(ST$trials))),
    varying = logical0_to_null(c(stats::na.omit(ST$cp_group))),
    sigma = all_pars[stringr::str_detect(all_pars, "^sigma_")],
    arma = all_pars[stringr::str_detect(all_pars, "(^ar|^ma)[0-9]")]
  )
  pars$reg = all_pars[!all_pars %in% c(pars$varying, pars$sigma, pars$arma)]
  pars$population = c(pars$reg, pars$sigma, pars$arma)

  # Check parameters
  # ARMA models
  if (length(pars$arma) > 0) {
    if (family$link %in% c("logit", "probit"))
      message("The current implementation of autoregression can be fragile for link='logit'. In particular, if there are any all-success trials (e.g., 10/10), the only solution is for 'ar' to be 0.00. If fitting succeeds, do a proper assessment of model convergence.")

    if (length(pars$sigma) > 1)
    message("You are using ar() together with sigma(). Autoregression usually assumes homoskedasticity (equal variance at all x). This is not a problem if if intercepts in ar() are joined by intercepts in sigma() like `~ [formula] + sigma(1) + ar(N)`. It may not be a problem for slopes on `sigma` either, but this has not been assessed thoroughly yet. So this is a note to be cautious about interpreting the sigma- and ar-parameters for now.")

    if (is.unsorted(data[, pars$x]) & is.unsorted(rev(data[, pars$x])))
      message("'", pars$x, "' is unordered. Please note that ar() applies in the order of data of the data frame - not the values.")
  }

  # Make formula_str and simulate
  formula_str_sim = get_all_formulas(ST, prior, pars$x, ytypes = c("ct", "sigma", "arma"))
  simulate = get_simulate(formula_str_sim, pars, nrow(ST), family)

  # Make jags code
  max_arma_order = get_arma_order(pars$arma)
  formula_str_jags = get_all_formulas(ST, prior, pars$x)
  jags_code_generated = get_jagscode(prior, ST, formula_str_jags, max_arma_order, family, sample)
  jags_code = ifelse(is.null(jags_code), yes = jags_code_generated, no = jags_code)  # Get from user?


  ##########
  # SAMPLE #
  ##########
  # Sample posterior
  if (sample %in% c("post", "both")) {
    samples = run_jags(
      data = data,
      jags_code = jags_code,
      pars = c(all_pars, "loglik_"),  # Monitor log-likelihood for loo/waic
      ST = ST,
      cores = cores,
      sample = "post",
      n.chains = chains,
      n.iter = iter,
      n.adapt = adapt,
      inits = inits
    )

    # Move loglik columns out to it's own list, keeping parameters and loglik apart
    loglik_cols = stringr::str_starts(colnames(samples[[1]]), 'loglik_')  # detect loglik cols
    mcmc_loglik = lapply(samples, function(x) x[, loglik_cols])
    mcmc_post = lapply(samples, function(x) x[, !loglik_cols])
  }

  # Sample prior
  if (sample %in% c("prior", "both")) {
    samples = run_jags(
      data = data,
      jags_code = jags_code,
      pars = all_pars,  # Not loglik
      ST = ST,
      cores = cores,
      sample = "prior",
      n.chains = chains,
      n.iter = iter,
      n.adapt = adapt,
      inits = inits
    )

    # Move loglik columns out to it's own list, keeping parameters and loglik apart
    loglik_cols = stringr::str_starts(colnames(samples[[1]]), 'loglik_')  # detect loglik cols
    mcmc_prior = lapply(samples, function(x) x[, !loglik_cols])
  }


  ##########
  # RETURN #
  ##########
  # Fill in the missing samples
  if (exists("mcmc_post")) class(mcmc_post) = "mcmc.list"
  if (exists("mcmc_prior")) class(mcmc_prior) = "mcmc.list"
  if (exists("mcmc_loglik")) class(mcmc_loglik) = "mcmc.list"
  if (!exists("mcmc_post")) mcmc_post = NULL
  if (!exists("mcmc_prior")) mcmc_prior = NULL
  if (!exists("mcmc_loglik")) mcmc_loglik = NULL

  # Make mrpfit object
  mcpfit = list(
    # By user (same order as mcp argument)
    segments = lapply(ST$form, stats::as.formula, env=globalenv()),  # with explicit response and cp
    data = data,
    prior = prior,
    family = family,

    # Results
    mcmc_post = mcmc_post,
    mcmc_prior = mcmc_prior,
    mcmc_loglik = mcmc_loglik,
    loo = NULL,
    waic = NULL,

    # Extracted model
    pars = pars,
    jags_code = jags_code,
    simulate = simulate,

    # Pass info to *.mcpfit() functions.
    # Not meant to be used by the end user.
    .other = list(
      ST = ST
    )
  )
  class(mcpfit) = "mcpfit"

  # Return it
  mcpfit
}
