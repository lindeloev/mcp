source("R/make_jagscode.R")
source("R/run_jags.R")
source("R/unpack_segments.R")

#' Fit Multiple Linear Segments And Their Change Points
#'
#' Change points are forced to be ordered using truncation.
#'
#' @param data Data.frame or tibble in long format.
#' @param segments List of formulas. Break points are estimated in between. The left-hand side specifies the chainge points and the right-hand side specifies the linear formula.
#' @param prior Named list of parameters and associated priors in JAGS code. Uninformative default priors are used where priors are not specified.
#' @param param_x String (default: NULL). Only relevant if no segments contains slope (no hint at what x is). Set this, e.g., param_x = "time".
#' @param sample Boolean (default: TRUE). Set to FALSE if you only want to check priors, the JAGS model, etc.
#' @param ... Parameters for `jags.parfit` which channels them to `jags.fit`.
#' @keywords mcmc, jags, mct
#' @export
#' @examples
#' # Define the segments that are separated by change points
#' segments = list(
#'   score ~ 1 + year,  # intercept + slope
#'    1 ~ 0 + year,  # joined slope
#'    1 ~ 0,  # joined plateau
#'    1 ~ 1  # disjoined plateau
#' )
#'
#' # Set priors.
#' # cp_i are change points.
#' # int_i are intercepts.
#' # x_i are slopes.
#' # i is the segment number (change points are to the right of the segment)
#' prior = list(
#'   int_1 = "dunif(10, 30)",  # intercept of segment 1
#'   cp_2 = "dunif(cp_1, 40),  # change point between segment 1 and 2. Must be greater than cp_1. Order restriction is applied automatically for everything but dunif (a JAGS limitation).
#'   year_2 = "dnorm(0, 1/5^2)  # slope of segment 1. Mean = 0, SD = 5.
#' )
#'
#' # Start sampling
#' fit = mcp(segments, data, prior)
#'
#' # Visual inspection of the results
#' plot(fit)
#' plot(fit, "combo")
#'
#' # Compare to an one-intercept-only model (no change points) with default prior
#' segments2 = list(1 ~ 1)
#' fit2 = mcp(segments2, data)  # fit another model here
#' fit$loo = loo(fit)
#' fit2$loo = loo(fit)
#' loo_compare(fit, fit2)
#'
#' # Show all priors (not just those specified manually)
#' fit$prior
#'
#' # Do stuff with the parameter estimates
#' fit$pars$model  # check out which parameters are inferred.
#' library(tidybayes)
#' spread_draws(fit$samples, cp_1, cp_2, int_1, year_1, year_2) %>%
#'    # tidybayes stuff here
#'
#' # Show JAGS model
#' cat(fit$jags_code)


mcp = function(segments, data, prior = list(), param_x = NULL, sample = TRUE, ...) {

  # Check input values
  if(sample) {
    if(!is.data.frame(data) & !tibble::is_tibble(data)) {
      stop("`data` must be a data.frame or a tibble.")
    }
  } else data = NULL  # define variable but nothing more
  if(!is.list(segments)) {
    stop("`segments` must be a list")
  }
  for(segment in segments) {
    if(!inherits(segment, "formula")) stop("all segments must be formulas.")
  }
  if(!is.list(prior)) {
    stop("`prior` must be a named list.")
  }
  if(!is.null(param_x) & !is.character(param_x)) {
    stop("`param_x` must be NULL or a string.")
  }
  if(!is.logical(sample)) {
    stop("`sample` must be TRUE or FALSE")
  }

  # Get prior, func_y, formula_jags, and param_x/param_y
  unpacked = unpack_segments(segments, prior, param_x)
  prior = unpacked$prior


  # Build model
  jags_code = make_jagscode(
    data = data,
    prior = prior,
    formula_jags = unpacked$formula_jags,
    nsegments = length(segments),
    sample = sample,
    param_x = unpacked$param_x,
    param_y = unpacked$param_y)

  # If samples should drawn. If not, just do everything else.
  if(sample) {
    # Sample it
    samples = run_jags(
      data = data,
      model = jags_code,
      params = c(names(prior), "loglik_"),
      ...
    )

    # Move loglik columns out to it's own list, keeping parameters and loglik apart
    loglik_cols = stringr::str_starts(colnames(samples[[1]]), 'loglik_')  # detect loglik cols
    loglik = lapply(samples, function(x) x[, loglik_cols])
    samples = lapply(samples, function(x) x[, !loglik_cols])
  } else {
    samples = NULL
    loglik = NULL
  }

  # Make mrpfit object
  mcpfit = list(
    samples = samples,
    loglik = loglik,

    loo = NULL,
    waic = NULL,

    data = data,
    prior = prior,

    pars = list(
      model = names(prior),
      x = unpacked$param_x,
      y = unpacked$param_y
    ),

    segments = segments,
    jags_code = jags_code,
    func_y = unpacked$func_y
  )
  class(mcpfit) = "mcpfit"

  # Return it
  mcpfit
}
