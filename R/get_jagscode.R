#' Make JAGS code for Multiple Change Point model
#'
#' @aliases get_jagscode
#' @inheritParams mcp
#' @param formula_str String. The formula string returned by \code{build_formula_str}.
#' @param ST Segment table. Returned by \code{get_segment_table()}.
#' @return String. A JAGS model.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_jagscode = function(data, prior, ST, formula_str) {
  # Begin building JAGS model. `mm` is short for "mcp model".
  # Add fixed variables.
  mm = paste0("
data {
  # X values
  MINX = ", ifelse(!is.null(data), paste0("min(", ST$x[1], ")"), "[requires data]"), "
  MAXX = ", ifelse(!is.null(data), paste0("max(", ST$x[1], ")"), "[requires data]"), "
  MEANX = ", ifelse(!is.null(data), paste0("mean(", ST$x[1], ")"), "[requires data]"), "
  SDX = ", ifelse(!is.null(data), paste0("sd(", ST$x[1], ")"), "[requires data]"), "

  # Y values
  MINY = ", ifelse(!is.null(data), paste0("min(", ST$y[1], ")"), "[requires data]"), "
  MAXY = ", ifelse(!is.null(data), paste0("max(", ST$y[1], ")"), "[requires data]"), "
  MEANY = ", ifelse(!is.null(data), paste0("mean(", ST$y[1], ")"), "[requires data]"), "
  SDY = ", ifelse(!is.null(data), paste0("sd(", ST$y[1], ")"), "[requires data]"), "
}

model {
  # Priors for population-level effects\n")
  ##########
  # PRIORS #
  ##########
  # Transform priors from SD to precision
  all_dprec = "dnorm|dt|dcauchy|ddexp|dlogis|dlnorm"  # JAGS precision dists
  for (i in seq_len(length(prior))) {
    if (stringr::str_detect(prior[[i]], all_dprec)) {
      prior[[i]] = sd_to_prec(prior[[i]])
    }
  }

  # Split up priors into population and varying
  prior_pop = prior[!names(prior) %in% ST$cp_group]
  prior_varying = prior[names(prior) %in% ST$cp_group]

  # Use get_prior_str() to add population-level priors
  for (i in 1:length(prior_pop)) {
    mm = paste0(mm, get_prior_str(prior_pop, i))
  }

  # Helpers for change points:
  mm = paste0(mm, "  cp_0 = -10^100  # mcp helper value; minus infinity\n")
  mm = paste0(mm, "  cp_", nrow(ST), " = 10^100  # mcp helper value; plus infinity\n")

  # Use get_prior_str() to add varying priors
  mm = paste0(mm, "\n  # Priors for varying effects\n")
  if (length(prior_varying) > 0) {
    for (i in 1:length(prior_varying)) {
      mm = paste0(mm, get_prior_str(
        prior = prior_varying,
        i = i,
        varying_group = na.omit(ST$cp_group_col[ST$cp_group == names(prior_varying[i])])
      ))
    }
  }


  ###########
  # FORMULA #
  ###########
  # Transform formula_str into JAGS format. Insert par_x and varying indices
  formula_jags = gsub("PAR_X", paste0(ST$x[1], "[i_]"), formula_str)
  for (i in seq_len(nrow(ST))) {
    formula_jags = gsub(paste0("CP_", i, "_INDEX"), paste0("[", ST$cp_group_col[i], "[i_]]"), formula_jags)
  }

  # Insert formula_jags
  mm = paste0(mm, "

  # Model and likelihood
  for (i_ in 1:length(", ST$x[1], ")) {

    # Fitted value
    y_[i_] = \n")

  # Add JAGS code for fitted values and indent it
  mm = paste0(mm, gsub("\n", "\n      ", formula_jags))

  # Finally the likelihood
  mm = paste0(mm, "\n\n    # Likelihood
    ", ST$y[1], "[i_] ~ dnorm(y_[i_], 1 / sigma^2)")


  ###############
  # OTHER STUFF #
  ###############

  # Log-density.
  mm = paste0(mm, "\n\n    # Log-density for LOO/WAIC computation
    loglik_[i_] = logdensity.norm(", ST$y[1], "[i_], y_[i_], 1 / sigma^2)")

  # Finish up
  mm = paste0(mm, "
  }
}")

  # Return the model
  mm
}


#' Get JAGS code for a prior
#'
#' @aliases get_prior_str
#' @inheritParams mcp
#' @param i The index in \code{prior} to get code for
#' @param varying_group String or NULL (default). Null indicates a population-
#'   level prior. String indicates a varying-effects prior (one for each group
#'   level).
#' @return A string
#'
get_prior_str = function(prior, i, varying_group = NULL) {
  # Helpers
  value = prior[i]
  name = names(prior[i])

  # Is this fixed?
  all_d = "dunif|dbern|dbeta|dbin|dchisqr|ddexp|dexp|df|dgamma|dgen.gamma|dhyper|dlogis|dlnorm|dnegbin|dnchisqr|dnorm|dpar|dpois|dt|dweib"  # All JAGS distributions
  is_fixed = stringr::str_detect(value, "^[-0-9.]+$") |
    value %in% names(prior)

  # If not either number or known parameter, it should be a known distribution.
  if (!is_fixed & !stringr::str_detect(value, all_d))
    stop("The prior \"", name, " = ", value, "\" is not a known distribution, a number, nor a model parameter.")

  # If it is a known distribution
  if (!is_fixed) {
    # ... and this is a population-level effect
    if (is.null(varying_group)) {
      return(paste0("  ", name, " ~ ", value, "\n"))
    } else {
      # It is a varying effect!
      return(paste0("  for (", varying_group, "_ in 1:n_unique_", varying_group, ") {
    ", name, "_uncentered[", varying_group, "_] ~ ", value, "
  }
  ", name, " = ", name, "_uncentered - mean(", name, "_uncentered)  # vectorized zero-centering\n"))
    }
  } else {
    # Fixed value. Just equate name and value
    return(paste0("  ", name, " = ", value, "  # Fixed\n"))
  }
}



#'¨Transform a prior from SD to precision.
#'
#' JAGS uses precision rather than SD. This function converts
#' \code{dnorm(4.2, 1.3)} into \code{dnorm(4.2, 1/1.3^2)}. It allows users to specify
#' priors using SD and then it's transformed for the JAGS code. It works for the
#' following distributions: dnorm|dt|dcauchy|ddexp|dlogis|dlnorm. In all of
#' these,
#' tau/sd is the second parameter.
#'
#'@aliases sd_to_prec
#'@param prior_str String. A JAGS prior. Can be truncated, e.g.
#'  \code{dt(3, 2, 1) T(my_var, )}.
#'@return A string
#'@author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'@export
#'

sd_to_prec = function(prior_str) {
  if (stringr::str_detect(prior_str, "dnorm|dt|dcauchy|ddexp|dlogis|dlnorm")) {
    # Identify trunc. If present, cut it out of the prior_str
    trunc_start = gregexpr("T\\(", prior_str)[[1]]
    if (trunc_start != -1) {
      trunc = substr(prior_str, trunc_start, 1000)
      prior_str = substr(prior_str, 0, trunc_start-1)
    } else trunc = ""

    # Split by commas to get distribution arguments.
    # Then convert the second to precision
    pieces = stringr::str_trim(strsplit(prior_str, ",")[[1]])
    pieces = gsub(" ", "", pieces)  # Remove spaces. Just so we know what we have.

    # In two-parameter distributions, tau is followed by a parenthesis.
    # Remove it so we just have the "value"
    is_dt = stringr::str_starts(pieces[1], "dt\\(")
    if (!is_dt) {
      pieces[2] = substr(pieces[2], 1, nchar(pieces[2])-1)
    }
    pieces[2] = paste0("1/(", pieces[2], ")^2")  # convert to precision

    # Stitch together again and return
    new_prior = paste0(paste0(pieces, collapse = ", "), ifelse(!is_dt, ") ", " "), trunc)
    return(new_prior)
  }
  else return(prior_str)
}
