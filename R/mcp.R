source("R/make_jagscode.R")
source("R/run_jags.R")

#' Fit Multiple Linear Segments And Their Change Points
#'
#' @param data A data.frame or tibble containing the variables expressed in `model` in long format.
#' @param ... Parameters for `jags.parfit` which channels them to `jags.fit`.
#' @keywords mcmc, jags, mct
#' @export
#' @examples
#' run_jags(data, model, params)
#' 

mcp = function(data, segments, prior = list(), ...) {
  
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
  model_obj = make_jagscode(data, segments, prior)
  
  # Sample it
  mcmc = run_jags(
    data = data, 
    model = model_obj$model, 
    params = c(unlist(model_obj$all_pars), model_obj$par_name_x, 'y_', 'sigma'),
    ...
  )
  
  # Return it
  mcmc
}