source("R/make_jagscode.R")
source("R/run_jags.R")
source("R/unpack_segments.R")

#' Fit Multiple Linear Segments And Their Change Points
#'
#' Given a list of linear segments, \code{mcp} infers the posterior
#' distributions of the parameters of each segment as well as the change points
#' between segments. All segments must regress on the same x-axis. Change points
#' are forced to be ordered using truncation of the priors. You can run
#' \code{fit = mcp(segments, sample=FALSE)} to avoid sampling and the need for
#' data if you just want to get the priors (\code{fit$prior}), the JAGS code
#' \code{fit$jags_code}, or the R function to simulate data (\code{fit$func_y}).
#'
#' @aliases mcp
#' @param data Data.frame or tibble in long format.
#' @param segments A list of formulas - one for each segment. The right-hand
#'   side specifices the form of intercepts and slopes. For the first segment,
#'   the left-hand side is the response variable. In the following segments, the
#'   left-hand side is the change point (on x). See examples for more details.
#' @param prior Named list. Names are parameter names (cp_i, int_i, [x_var]_i,
#'   sigma) and the values are the associated priors in JAGS code. Uninformative
#'   default priors are used where priors are not specified.
#'   \code{mct} uses SD (not precision) for dnorm, dt, dlogis, etc. See details.
#'   Change points are forced to be ordered through the priors using truncation,
#'   \code{dnorm(0, 1) T(cp_1, )}, except for uniform priors where the lower
#'   bound should be greater than the previous change point, \code{dunif(cp_1, )}.
#' @param param_x String (default: NULL). Only relevant if no segments contains
#'   slope (no hint at what x is). Set this, e.g., param_x = "time".
#' @param sample Boolean (default: TRUE). Set to FALSE if you only want to check
#'   priors, the JAGS model, etc.
#' @param family WORK IN PROGRESS. One of "gauss" (default), "binomial"
#' @param ... Parameters for \code{jags.parfit} which channels them to \code{jags.fit}.
#' @details Noites on priors:
#'   * Order restriction is automatically applied to cp_\* parameters using
#'     truncation (e.g., T(cp_1, )) so that they are in the correct order on the
#'     x-axis UNLESS you do it yourself. The one exception is for dunif
#'     distributions where you have to do it as above.
#'   * In addition to the model parameters, \code{MINX} (minimum x-value), \code{MAXX}
#'     (maximum x-value), \code{SDX} (etc...), \code{MINY}, \code{MAXY}, and \code{SDY}
#'     are also available when you set priors. They are used to set uninformative
#'     default priors.
#'   * Use SD when you specify priors for dt, dlogis, etc. JAGS uses precision
#'     but mct converts to precision under the hood via the sd_to_prec()
#'     function. So you will see SDs in \code{fit$prior} but precision ($1/SD^2)
#'     in \code{fit$jags_model}
#' @return An \code{mcpfit} object.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \dontrun{
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
#'   int_1 = "dt(10, 30) T(0, )",  # t-dist intercept. Truncated to > 0
#'   cp_2 = "dunif(cp_1, 40),  # change point to segment 2 > cp_1.
#'   year_2 = "dnorm(0, 5)  # slope of segment 1. Mean = 0, SD = 5.
#' )
#'
#' # Start sampling
#' fit = mcp(segments, data, prior)
#'
#' # Visual inspection of the results
#' plot(fit)
#' plot(fit, "combo")
#'
#' # Compare to a one-intercept-only model (no change points) with default prior
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
#' }


mcp = function(segments, data = NULL, prior = list(), family = "gaussian", param_x = NULL, sample = TRUE, ...) {

  # Check input values
  if(is.null(data) & sample == TRUE) {
    stop("Cannot sample without data.")
  }
  if(sample == TRUE) {
    if(!is.data.frame(data) & !tibble::is_tibble(data)) {
      stop("`data` must be a data.frame or a tibble.")
    }
  } else {
    data = NULL  # define variable but nothing more
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
  if(!is.null(param_x) & !is.character(param_x)) {
    stop("`param_x` must be NULL or a string.")
  }
  if(!is.logical(sample)) {
    stop("`sample` must be TRUE or FALSE")
  }

  # Get prior, func_y, formula_jags, and param_x/param_y
  unpacked = unpack_segments(segments, prior, param_x)
  prior = unpacked$prior

  # Check variables in data
  if(sample == TRUE) {
    if(!unpacked$param_x %in% colnames(data)) {
      stop(paste0("The slope variable ", unpacked$param_x, " is not a column in data."))
    }
    if(!unpacked$param_y %in% colnames(data)) {
      stop(paste0("The response variable", unpacked$param_y, " is not a column in data."))
    }
  }


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
      jags_code = jags_code,
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
