#' Fit Multiple Linear Segments And Their Change Points
#'
#' Given a list of linear segments, `mcp` infers the posterior
#' distributions of the parameters of each segment as well as the change points
#' between segments. [See more details and worked examples on the mcp website](https://lindeloev.github.io/mcp/).
#' All segments must regress on the same x-variable. Change
#' points are forced to be ordered using truncation of the priors. You can run
#' `fit = mcp(segments, sample=FALSE)` to avoid sampling and the need for
#' data if you just want to get the priors (`fit$prior`), the JAGS code
#' `fit$jags_code`, or the R function to simulate data (`fit$func_y`).
#'
#' @aliases mcp
#' @param data Data.frame or tibble in long format.
#' @param segments A list of formulas - one for each segment. The first formula
#'   has the format `response ~ predictors` while the following formulas
#'   have the format `response ~ changepoint ~ predictors`. The response
#'   and change points can be omitted (`changepoint ~ predictors` assumes same
#'   response. `~ predictors` assumes an intercept-only change point).
#'
#'   See examples for more details.
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
#'   * 1: serial sampling
#'   * >1: parallel sampling on this number of cores. Ideally set `chains`
#'     to the same value.
#'   * "all": use all cores but one.
#' @param chains Positive integer. Number of chains to run.
#' @param iter Positive integer. Number of post-warmup samples to draw.
#' @param adapt Positive integer. Number of iterations to find sampler settings.
#' @param update Positive integer. Also sometimes called "burnin", this is the
#'   number of regular samples before anything is recorded. Use to reach
#'   convergence.
#' @param jags_explicit Pass JAGS code to `mcp` to use directly. Useful if
#'   you want to make small tweaks, but mostly used for the development of mcp.
#' @param ... Further parameters for \code{\link[dclone]{jags.fit}}.
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
#' @return An `mcpfit` object.
#' @seealso \link{get_segment_table}
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @importFrom stats gaussian binomial
#' @export
#' @examples
#' \dontrun{
#' # Define the segments that are separated by change points
#' segments = list(
#'   score ~ 1 + year,  # intercept + slope
#'    ~ 0 + year,  # joined slope
#'    ~ 0,  # joined plateau
#'    ~ 1  # disjoined plateau
#' )
#'
#' # Sample and see results
#' fit = mcp(segments, data)
#' summary(fit)
#'
#' # Visual inspection of the results
#' plot(fit)
#' plot(fit, "combo")
#'
#' # Test a hypothesis
#' hypothesis(fit, "cp_1 > 10")
#'
#' # Compare to a one-intercept-only model (no change points) with default prior
#' segments2 = list(score ~ 1)
#' fit2 = mcp(segments2, data)  # fit another model here
#' fit$loo = loo(fit)
#' fit2$loo = loo(fit)
#' loo_compare(fit, fit2)
#'
#' # Sample the prior and inspect it using all the usual methods (prior predictive checks)
#' fit_prior = mcp(segments, data, sample = "prior")
#' summary(fit_prior)
#' plot(fit_prior)
#'
#' # Show all priors. Default priors are added where you don't provide any
#' fit$prior
#'
#' # Set priors and re-run
#' # cp_i are change points.
#' # int_i are intercepts.
#' # x_i are slopes.
#' # i is the segment number (change points are to the right of the segment)
#' prior = list(
#'   int_1 = "dt(10, 30) T(0, )",  # t-dist intercept. Truncated to > 0
#'   year_1 = "dnorm(0, 5)",  # slope of segment 1. Mean = 0, SD = 5.
#'   cp_2 = "dunif(cp_1, 40)",  # change point to segment 2 > cp_1.
#'   year_2 = "year_1",  # Shared slope between segment 2 and 1
#'   int_3 = 15  # Fixed intercept of segment 3
#' )
#' fit3 = mcp(segments, data, prior)
#'
#' # Show JAGS model
#' cat(fit$jags_code)
#' }

mcp = function(segments,
               data = NULL,
               prior = list(),
               family = gaussian(),
               par_x = NULL,
               sample = "post",
               cores = 1,
               chains = 3,
               iter = 3000,
               adapt = 1500,
               update = 1500,
               jags_explicit = NULL,
               ...) {

  # Check data
  if (is.null(data) & sample == TRUE)
    stop("Cannot sample without data.")

  if (!is.null(data)) {
    if (!is.data.frame(data) & !tibble::is_tibble(data))
      stop("`data` must be a data.frame or a tibble.")
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

  # Check family
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
  # Throw a warning
  major = as.numeric(R.Version()$major)
  minor = as.numeric(R.Version()$minor)
  fails_parallel = (major < 3 | (major == 3 & minor < 6.1))
  if (cores > 1 & fails_parallel == TRUE)
    warning("Parallel sampling (`cores` > 1) has been shown to err on R versions below 3.6.1. You have ", R.Version()$version.string, ". Consider upgrading if it fails or hangs.")


  # Get an abstract segment table ("ST")
  ST = get_segment_table(segments, data, family$family, par_x)

  par_x = unique(ST$x)
  par_y = unique(ST$y)
  par_trials = unique(ST$trials)

  # Get prior and lists of parameters
  prior = get_prior(ST, family$family, prior, ar_order)
  pars_varying = logical0_to_null(c(stats::na.omit(ST$cp_group)))
  pars_population = names(prior)[!names(prior) %in% pars_varying]  # Simply the absence of varying pars

  # Make formula_str and func_y
  formula_str = get_formula_str(ST, par_x)
  if (family$family == "gaussian") {
    formula_str_sigma = get_formula_str(ST, par_x, sigma = TRUE)
    formula_str = paste0(formula_str, "\n\n", formula_str_sigma)
  }

  pars_funcy = pars_population[!pars_population %in% ST$cp_sd]
  func_y = get_func_y(formula_str, par_x, par_trials, pars_funcy, pars_varying, nrow(ST), family$family)

  # Make jags code and sample it.
  jags_code = get_jagscode(prior, ST, formula_str, family$family, sample)

  # Sample posterior
  if (sample %in% c("post", "both")) {
    samples = run_jags(
      data = data,
      jags_code = ifelse(is.null(jags_explicit), jags_code, jags_explicit),
      pars = c(pars_population, pars_varying, "loglik_"),  # population-level, varying, and loglik for loo/waic
      ST = ST,
      cores = cores,
      sample = "post",
      n.chains = chains,
      n.iter = iter,
      n.adapt = adapt,
      n.update = update,
      ...
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
      jags_code = ifelse(is.null(jags_explicit), jags_code, jags_explicit),
      pars = c(pars_population, pars_varying),  # population-level, varying, but NOT loglik
      ST = ST,
      cores = cores,
      sample = "prior",
      n.chains = chains,
      n.iter = iter,
      n.adapt = adapt,
      n.update = update,
      ...
    )

    # Move loglik columns out to it's own list, keeping parameters and loglik apart
    loglik_cols = stringr::str_starts(colnames(samples[[1]]), 'loglik_')  # detect loglik cols
    mcmc_prior = lapply(samples, function(x) x[, !loglik_cols])
  }

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
    segments = lapply(ST$form, as.formula, env=globalenv()),  # with explicit response and cp
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
    pars = list(
      population = pars_population,
      varying = pars_varying,
      x = par_x,
      y = par_y,
      trials = par_trials
    ),

    jags_code = jags_code,
    func_y = func_y,

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


#' Converts logical(0) to null. Returns x otherwise
#'
#'@aliases logical0_to_null
#'@param x Anything
#'@return NULL or x

logical0_to_null = function(x) {
  if (length(x) > 0)
    return(x)
  else return(NULL)
}
