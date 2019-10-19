#' Make JAGS code for Multiple Change Point model
#'
#' @aliases get_jagscode
#' @inheritParams mcp
#' @param formula_str String. The formula string returned by \code{build_formula_str}.
#' @param nsegments Number of segments
#' @param param_y String. Name of response parameter.
#' @return String. A JAGS model.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_jagscode = function(data, prior, formula_str, nsegments, param_x, param_y) {
  # Transform priors from SD to precision
  for(i in 1:length(prior)) {
    if(stringr::str_detect(prior[[i]], "dnorm|dt|dcauchy|ddexp|dlogis|dlnorm")) {
      prior[[i]] = sd_to_prec(prior[[i]])
    }
  }

  # Transform formula_str into JAGS format
  formula_jags = gsub("PARAM_X", paste0(param_x, "[i_]"), formula_str)

  # Begin building JAGS model. `mm` is short for "mcp model".
  # Add fixed variables.
  mm = paste0("
data {
  # X values
  MINX = ", ifelse(!is.null(data), min(data[, param_x]), "[requires data]"), "
  MAXX = ", ifelse(!is.null(data), max(data[, param_x]), "[requires data]"), "
  MEANX = ", ifelse(!is.null(data), mean(unlist(data[, param_x])), "[requires data]"), "
  SDX = ", ifelse(!is.null(data), sd(unlist(data[, param_x])), "[requires data]"), "

  # Y values
  MINY = ", ifelse(!is.null(data), min(data[, param_y]), "[requires data]"), "
  MAXY = ", ifelse(!is.null(data), max(data[, param_y]), "[requires data]"), "
  MEANY = ", ifelse(!is.null(data), mean(unlist(data[, param_y])), "[requires data]"), "
  SDY = ", ifelse(!is.null(data), sd(unlist(data[, param_y])), "[requires data]"), "
}

model {
  # Priors
  ")

  # Add all priors
  # First as vector and code whether it's fixed or not
  prior_vector = unlist(prior)
  is_fixed = stringr::str_detect(prior_vector, "^[-0-9.]+$")

  # Distributed for non-fixed. Equal sign for fixed.
  if(sum(!is_fixed) > 0) {
    prior_str_dist = paste0(names(prior_vector[!is_fixed]), " ~ ", prior_vector[!is_fixed], collapse="\n    ")
  } else prior_str_dist = ""
  if(sum(is_fixed) > 0) {
    prior_str_fix = paste0(names(prior_vector[is_fixed]), " = ", prior_vector[is_fixed], collapse="  # Fixed\n    ")
  } else prior_str_fix = ""

  # Now add them in!
  mm = paste0(mm, prior_str_dist, "\n    ", prior_str_fix)
  mm = paste0(mm, "
  cp_0 = -10^100  # mcp helper value; minus infinity
  cp_", nsegments, " = 10^100  # mcp helper value; plus infinity\n
    ")


  # * Build list of parameters involved in each segment
  mm = paste0(mm, "\n
  # Model and likelihood
  for(i_ in 1:length(", param_x, ")) {

    # Fitted value
    y_[i_] = \n")

  # Add JAGS code for fitted values and indent it
  mm = paste0(mm, gsub("\n", "\n      ", formula_jags))

  # Finally the likelihood
  mm = paste0(mm, "\n\n    # Likelihood
    ", param_y, "[i_] ~ dnorm(y_[i_], 1 / sigma^2)")

  # Log-density.
  mm = paste0(mm, "\n\n    # Log-density for LOO/WAIC computation
    loglik_[i_] = logdensity.norm(", param_y, "[i_], y_[i_], 1 / sigma^2)")

  # Finish up
  mm = paste0(mm, "
  }
}")

  # Return the model
  mm
}
