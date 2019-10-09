source("R/make_jagscode.R")
source("R/run_jags.R")

#' Fit Multiple Linear Segments And Their Change Points
#'
#' Change points are forced to be ordered using truncation.
#'
#' @param data A data.frame or tibble in long format.
#' @param segments A list of formulas. Break points are estimated in between. The left-hand side specifies the chainge points and the right-hand side specifies the linear formula.
#' @param prior A named list of parameters and associated priors in JAGS code. Uninformative default priors are used where priors are not specified.
#' @param param_x A string. Only relevant if no segments contains slope (no hint at what x is). Set this, e.g., param_x = "time".
#' @param ... Parameters for `jags.parfit` which channels them to `jags.fit`.
#' @keywords mcmc, jags, mct
#' @import stringr
#' @export
#' @examples
#' segments = list(
#'   score ~ 1 + year,  # intercept + slope
#'    1 ~ 0 + year,  # joined slope
#'    1 ~ 0,  # joined plateau
#'    1 ~ 1  # disjoined plateau
#' )
#'
#' prior = list(
#'   int_1 = "dunif(10, 30)",  # intercept of segment 1
#'   cp_2 = "dunif(cp_1, 40),  # change point between segment 1 and 2. Must be greater than cp_1. Order restriction is applied automatically for everything but dunif (a JAGS limitation).
#'   year_2 = "dnorm(0, 1/5^2)  # slope of segment 1. Mean = 0, SD = 5.
#' )
#'
#' fit = mcp(data, segments, prior)
#'
#' # See results
#' plot(fit)
#' plot(fit, "combo")
#'
#' # Compare models predictive performance
#' fit2 = ...  # fit another model here
#' fit$loo = loo(fit)
#' fit2$loo = loo(fit)
#' loo_compare(fit, fit2)
#'
#' # Show all priors (not just those specified manually)
#' fit$prior
#'
#' # Show JAGS model
#' cat(fit$model_jags)
#'
#'


mcp = function(data, segments, prior = list(), param_x = NULL, ...) {

  # Check input values
  if(!is.data.frame(data) & !is.tibble(data)) {
    stop("`data` must be a data.frame or a tibble.")
  }
  if(!is.list(segments)) {
    stop("`segments` must be a list")
  }
  for(segment in segments) {
    if(!inherits(segment, "formula")) stop("all segments must be formulas.")
  }
  if(!is.list(prior)) {
    stop("`prior` must be a named list.")
  }

  # Build model
  model_obj = make_jagscode(data, segments, prior, param_x = param_x)

  # Sample it
  samples = run_jags(
    data = data,
    model = model_obj$model,
    params = c(unlist(model_obj$all_pars), model_obj$param_name_x, "y_", "sigma", "loglik_"),
    ...
  )

  # Split loglik columns way
  loglik_cols = str_starts(colnames(samples[[1]]), 'loglik_')  # detect loglik cols
  loglik = lapply(samples, function(x) x[, loglik_cols])
  samples = lapply(samples, function(x) x[, !loglik_cols])

  # Make mrpfit object
  mcpfit = list(
    samples = samples,
    loglik = loglik,

    loo = NULL,
    waic = NULL,

    data = data,
    prior = prior,

    pars = list(
      model = model_obj$all_pars,
      x = model_obj$param_name_x,
      y = model_obj$param_name_y
    ),

    segments = segments,
    model_jags = model_obj$model
  )
  class(mcpfit) = "mcpfit"

  # Return it
  mcpfit
}
