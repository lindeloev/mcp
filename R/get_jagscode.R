#' Make JAGS code for Multiple Change Point model
#'
#' @aliases get_jagscode
#' @keywords internal
#' @inheritParams mcp
#' @param formula_str String. The formula string returned by `build_formula_str`.
#' @param ST Segment table. Returned by `get_segment_table()`.
#' @param arma_order Positive integer. The autoregressive order.
#' @return String. A JAGS model.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_jagscode = function(prior, ST, formula_str, arma_order, family, sample) {
  # Begin building JAGS model. `mm` is short for "mcp model".
  # Add fixed variables.
  mm = paste0("
model {
  # Priors for population-level effects\n")
  ##########
  # PRIORS #
  ##########
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
  if (length(prior_varying) > 0) {
    mm = paste0(mm, "\n  # Priors for varying effects\n")
    for (i in 1:length(prior_varying)) {
      mm = paste0(mm, get_prior_str(
        prior = prior_varying,
        i = i,
        varying_group = stats::na.omit(ST$cp_group_col[ST$cp_group == names(prior_varying[i])])
      ))
    }
  }

  # Autocorrelation: detect if there is an intercept or slope on AR
  y_code = "y_[i_]"
  has_ar = !all(is.na(unlist(ST$ar_code))) | !all(is.na(unlist(ST$ar_int)))
  if (has_ar) {
    # Add computation of autocorrelated residuals
    mm = paste0(mm, get_ar_code(arma_order, ST, family))
    y_code = paste0(y_code, " + sum(ar_[i_, ])")
  }

  # Add inverse link function to back-transform to observed metric
  if (family$link != "identity")  # not identity
    y_code = paste0(family$linkinv_jags, "(", y_code, ")")




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
  for (i_ in 1:length(", ST$x[1], ")) {")

  # Add JAGS code for fitted values and indent it
  mm = paste0(mm, gsub("\n", "\n    ", formula_jags))

  # Finally the likelihood
  mm = paste0(mm, "\n\n    # Likelihood and log-density for family = ", family$family, "()
    ")

  if (family$family == "gaussian") {
    mm = paste0(mm, ST$y[1], "[i_] ~ dnorm(", y_code, ", 1 / sigma_[i_]^2)
    loglik_[i_] = logdensity.norm(", ST$y[1], "[i_], ", y_code, ", 1 / sigma_[i_]^2)")

  } else if (family$family == "binomial") {
    mm = paste0(mm, ST$y[1], "[i_] ~ dbin(", y_code, ", ", ST$trials[1], "[i_])
    loglik_[i_] = logdensity.bin(", ST$y[1], "[i_], ", y_code, ", ", ST$trials[1], "[i_])")
  } else if (family$family == "bernoulli") {
    mm = paste0(mm, ST$y[1], "[i_] ~ dbern(", y_code, ")
    loglik_[i_] = logdensity.bern(", ST$y[1], "[i_], ", y_code, ")")
  } else if (family$family == "poisson") {
    mm = paste0(mm, ST$y[1], "[i_] ~ dpois(", y_code, ")
    loglik_[i_] = logdensity.pois(", ST$y[1], "[i_], ", y_code, ")")
  }

  # If only the prior is sampled, remove the loglik_[i_] line
  if (sample == "prior")
    mm = gsub("loglik.*?$","", mm)


  ###############
  # OTHER STUFF #
  ###############

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
#' @keywords internal
#' @inheritParams mcp
#' @param i The index in `prior` to get code for
#' @param varying_group String or NULL (default). Null indicates a population-
#'   level prior. String indicates a varying-effects prior (one for each group
#'   level).
#' @return A string
#'
get_prior_str = function(prior, i, varying_group = NULL) {
  # Helpers
  value = prior[[i]]
  name = names(prior[i])

  # Is this fixed?
  all_d = "dunif|dbern|dbeta|dbin|dchisqr|ddexp|dexp|df|dgamma|dgen.gamma|dhyper|dlogis|dlnorm|dnegbin|dnchisqr|dnorm|dpar|dpois|dt|dweib"  # All JAGS distributions
  is_fixed = stringr::str_detect(value, "^[-0-9.]+$") |
    value %in% names(prior)

  # If not either number or known parameter, it should be a known distribution.
  if (!is_fixed & !stringr::str_detect(value, all_d))
    stop("The prior '", name, " = ", value, "' is not a known distribution, a number, nor a model parameter.")

  # If it is a known distribution
  if (!is_fixed) {
    # Convert to precision
    value = sd_to_prec(value)

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



#' Transform a prior from SD to precision.
#'
#' JAGS uses precision rather than SD. This function converts
#' `dnorm(4.2, 1.3)` into `dnorm(4.2, 1/1.3^2)`. It allows users to specify
#' priors using SD and then it's transformed for the JAGS code. It works for the
#' following distributions: dnorm|dt|dcauchy|ddexp|dlogis|dlnorm. In all of
#' these,
#' tau/sd is the second parameter.
#'
#'@aliases sd_to_prec
#'@param prior_str String. A JAGS prior. Can be truncated, e.g.
#'  `dt(3, 2, 1) T(my_var, )`.
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



#' Get JAGS code for residuals for each AR-order
#'
#' @aliases ar_code
#' @keywords internal
#' @inheritParams get_jagscode
#' @param ar_order The order of the autoregressive component
get_ar_code = function(ar_order, ST, family) {
  code = ""
  for (i in seq_len(ar_order)) {
    # Get code for link(y[i - order])
    if (family$family != "binomial") {
      y_obs = paste0(family$link_jags, "(", ST$y[1], "[i_ - ", i, "])")
    } else {
      y_obs = paste0(family$link_jags, "(", ST$y[1], "[i_ - ", i, "] / ", ST$trials[1], "[i_ - ", i, "])")
    }

    code = paste0(code, "

  # AR(", i, ") on residuals:
  ar_[1:", i, ", ", i, "] = c(", paste0(rep("0", i), collapse = ","), ")
  for (i_ in ", i + 1, ":length(", ST$x[1], ")) {ar_[i_, ", i, "] = ar", i, "_[i_] * (", y_obs, " - y_[i_-", i, "])}")
  }
  return(code)
}
