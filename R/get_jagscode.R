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
  # Priors
  ")
  ###########################
  # POPULATION-LEVEL PRIORS #
  ###########################
  # Transform priors from SD to precision
  for (i in seq_len(length(prior))) {
    if (stringr::str_detect(prior[[i]], "dnorm|dt|dcauchy|ddexp|dlogis|dlnorm")) {
      prior[[i]] = sd_to_prec(prior[[i]])
    }
  }

  # Add all priors
  # First as vector and code whether it's fixed or not
  prior_vector = unlist(prior)
  is_fixed = stringr::str_detect(prior_vector, "^[-0-9.]+$")

  # Distributed for non-fixed. Equal sign for fixed.
  if (sum(!is_fixed) > 0) {
    prior_str_dist = paste0(names(prior_vector[!is_fixed]), " ~ ", prior_vector[!is_fixed], collapse = "\n    ")
  } else prior_str_dist = ""
  if (sum(is_fixed) > 0) {
    prior_str_fix = paste0(names(prior_vector[is_fixed]), " = ", prior_vector[is_fixed], collapse = "  # Fixed\n    ")
  } else prior_str_fix = ""

  # Now add them in!
  mm = paste0(mm, prior_str_dist, "\n  ", prior_str_fix)
  mm = paste0(mm, "
  cp_0 = -10^100  # mcp helper value; minus infinity
  cp_", nrow(ST), " = 10^100  # mcp helper value; plus infinity\n
    ")


  ##################
  # VARYING PRIORS #
  ##################
  for (i in seq_len(nrow(ST))) {
    S = ST[i, ]

    # Varying change points (intercepts)
    if (S$cp_ran_int == TRUE & !is.na(S$cp_group)) {
      # Special truncate for first and last change point to stay within MINX-MAXX
      trunc = paste0("T(", ifelse(i == 2, "MINX", ST$cp_code_prior[i-1]), ", ", ifelse(i == nrow(ST), "MAXX", ST$cp_code_prior[i+1]), ")")
      mm = paste0(mm, "
  # Zero-centered varying change point for segment ", i, "
  for (", S$cp_group_col, "_ in 1:n_unique_", S$cp_group_col, ") {
    ", S$cp_group, "_uncentered[", S$cp_group_col, "_] ~ dnorm(0, 1/", S$cp_sd, "^2) ", trunc, "
  }
  ", S$cp_group, " = ", S$cp_group, "_uncentered - mean(", S$cp_group, "_uncentered)  # vectorized zero-centering
  ")
    }
  }



  ###########
  # FORMULA #
  ###########
  # Transform formula_str into JAGS format. Insert par_x and varying indices
  formula_jags = gsub("PAR_X", paste0(ST$x[1], "[i_]"), formula_str)
  for(i in seq_len(nrow(ST))) {
    formula_jags = gsub(paste0("CP_", i, "_INDEX"), paste0("[", ST$cp_group_col[i], "[i_]]"), formula_jags)
  }

  # Insert formula_jags
  mm = paste0(mm, "\n

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
