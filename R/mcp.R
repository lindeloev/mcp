# source("R/get_segment_table.R")
# source("R/get_prior.R")
# source("R/get_formula.R")
# source("R/get_jagscode.R")
# source("R/run_jags.R")

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
#'   sigma) and the values are either
#'
#'    * A JAGS distribution (e.g., \code{int_1 = "dnorm(0, 1) T(0,)"}) indicating a
#'      conventional prior distribution. Uninformative priors based on data
#'      propertiesare used where priors are not specified. This ensures good
#'      parameter estimations, but it is a questionable for hypothesis testing.
#'      \code{mcp} uses SD (not precision) for dnorm, dt, dlogis, etc. See
#'      details. Change points are forced to be ordered through the priors using
#'      truncation, except for uniform priors where the lower bound should be
#'      greater than the previous change point, \code{dunif(cp_1, MAXX)}.
#'    * A numerical value (e.g., \code{int_1 = -2.1}) indicating a fixed value.
#'    * A model parameter name (e.g., \code{int_2 = "int_1"}), indicating that this parameter is shared -
#'      typically between segments. If two varying effects are shared this way,
#'      they will need to have the same grouping variable.
#' @param par_x String (default: NULL). Only relevant if no segments contains
#'   slope (no hint at what x is). Set this, e.g., par_x = "time".
#' @param sample Boolean (default: TRUE). Set to FALSE if you only want to check
#'   priors, the JAGS model, etc.
#' @param family WORK IN PROGRESS. One of "gauss" (default), "binomial"
#' @param jags_explicit Pass JAGS code to \code{mcp} to use directly. Useful if
#'   you want to make small tweaks, but mostly used for the development of mcp.
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
#'     but \code{mcp} converts to precision under the hood via the sd_to_prec()
#'     function. So you will see SDs in \code{fit$prior} but precision ($1/SD^2)
#'     in \code{fit$jags_code}
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
#' # Start sampling
#' fit = mcp(segments, data)
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
#' # Set priors and re-run
#' # cp_i are change points.
#' # int_i are intercepts.
#' # x_i are slopes.
#' # i is the segment number (change points are to the right of the segment)
#' prior = list(
#'   int_1 = "dt(10, 30) T(0, )",  # t-dist intercept. Truncated to > 0
#'   year_1 = "dnorm(0, 5)",  # slope of segment 1. Mean = 0, SD = 5.
#'   cp_2 = "dunif(cp_1, 40),  # change point to segment 2 > cp_1.
#'   year_2 = "year_1",  # Shared slope between segment 2 and 1
#'   int_3 = 15  # Fixed intercept of segment 3
#' )
#' fit3 = mcp(segments, data, prior)
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


mcp = function(segments, data = NULL, prior = list(), family = "gaussian", par_x = NULL, sample = TRUE, jags_explicit = NULL, ...) {

  # Check input values
  if (is.null(data) & sample == TRUE)
    stop("Cannot sample without data.")

  if (!is.null(data)) {
    if (!is.data.frame(data) & !tibble::is_tibble(data))
      stop("`data` must be a data.frame or a tibble.")
  }

  if (!is.list(segments))
    stop("`segments` must be a list")

  if (length(segments) == 0)
    stop("At least one segment is needed")

  for (segment in segments) {
    if (!inherits(segment, "formula"))
      stop("all segments must be formulas.")
  }

  if (!is.list(prior))
    stop("`prior` must be a named list.")

  if (any(duplicated(names(prior))))
    stop("`prior` has duplicated entries for the same parameter.")

  if (!is.null(par_x) & !is.character(par_x))
    stop("`par_x` must be NULL or a string.")

  if (!is.logical(sample))
    stop("`sample` must be TRUE or FALSE")


  # Get an abstract segment table ("ST")prior, func_y, formula_jags, and par_x/param_y
  ST = get_segment_table(segments, data, par_x)
  par_x = unique(ST$x)
  param_y = unique(ST$y)

  # Get prior and lists of parameters
  prior = get_prior(ST, prior)
  #params_population = names(prior)
  params_population = stats::na.omit(unique(c("sigma", ST$int_name, ST$slope_name, ST$cp_name[-1], ST$cp_sd)))
  params_varying = logical0_to_null(c(stats::na.omit(ST$cp_group)))

  # Make formula_str and func_y
  formula_str = get_formula_str(ST)

  params_funcy = params_population[!params_population %in% ST$cp_sd]
  func_y = get_func_y(formula_str, par_x, params_funcy, params_varying, nrow(ST))

  # Make jags code and sample if sample==TRUE. If not, just skip it.
  jags_code = get_jagscode(data, prior, ST, formula_str)

  if (sample) {
    # Sample it
    samples = run_jags(
      data = data,
      jags_code = ifelse(is.null(jags_explicit), jags_code, jags_explicit),
      params = c(params_population, params_varying, "loglik_"),  # population-level, varying, and loglik for loo/waic
      ST = ST,
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
      population = params_population,
      varying = params_varying,
      x = par_x,
      y = param_y
    ),

    segments = segments,
    jags_code = jags_code,
    func_y = func_y,

    # Not really meant to be used by the end user.
    # But useful to handle the class.
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
