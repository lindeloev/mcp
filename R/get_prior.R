##########################
# LIST OF DEFAULT PRIORS #
##########################

# Generic default priors that applies to many families
cp_prior = list(
  cp_1 = "dunif(MINX, MAXX)",  # If there is only one change point
  cp = "dt(MINX, (MAXX - MINX) / N_CP, N_CP - 1)",
  cp_rel = "dunif(0, MAXX - %s)",
  sd = "dnorm(0, 2 * (MAXX - MINX) / N_CP) T(0, )"
)

sigma_prior_identity = list(
  sigma_int = "dnorm(0, SDY) T(0, )",
  sigma_slope = "dt(0, SDY / (MAXX - MINX), 3)"
)

sigma_prior_log = list(
  sigma_int = "dnorm(log(SDY) / 2, log(SDY) / 2)",
  sigma_slope = "dt(0, log(SDY) / (MAXX - MINX), 3)"
)

arma_prior = list(
  arma_int = "dunif(-1, 1)",
  arma_slope = "dnorm(0, 1 / (MAXX - MINX))"  # 68% of changing 1 over the observed x values
)


# Per-family priors, mixing in the generic priors
priors = list(
  gaussian_identity = c(cp_prior, sigma_prior_identity, arma_prior, list(
    ct_int = "dt(0, 3 * SDY, 3)",
    ct_slope = "dt(0, SDY / (MAXX - MINX), 3)"
  )),

  gaussian_log = c(cp_prior, sigma_prior_log, arma_prior, list(
    ct_int = "dt(log(MEANY), log(SDY), 3)",  # The SD is inherently an overshoot here, so embodies the "3 * SDY" of the identity
    ct_slope = "dt(0, log(SDY) / (3 * (MAXX - MINX)), 3)"
  )),

  # A logit/probit of +/- 5 is quite extreme. Very compatible with 3
  binomial_logit = c(cp_prior, arma_prior, list(
    ct_int = "dnorm(0, 3)",
    ct_slope = "dnorm(0, 3 / (MAXX - MINX))"
  )),

  binomial_identity = c(cp_prior, arma_prior, list(
    ct_int = "dunif(0, 1)",
    ct_slope = "dnorm(0, 1 / (MAXX - MINX))"
  )),

  poisson_log = c(cp_prior, arma_prior, list(
    ct_int = "dnorm(0, 10)",
    ct_slope = "dnorm(0, 10)"
  )),

  poisson_identity = c(cp_prior, arma_prior, list(
    ct_int = "dt(MEANY, MEANY, 3)",  # Double-dipping data (bad). Try not to go below zero.
    ct_slope = "dt(0, SDY / (MAXX - MINX), 3)"
  )),

  exponential_identity = c(cp_prior, list(
    ct_int = "dgamma(0.1, 0.1)",
    ct_slope = "dnorm(0, log(SDY))"
  ))
)

# Logit and probit are so similar that the same prior applies
priors$binomial_probit = priors$binomial_logit

# Identical priors for binomial and bernoulli.
priors$bernoulli_logit = priors$binomial_logit
priors$bernoulli_probit = priors$binomial_probit
priors$bernoulli_identity = priors$binomial_identity




######################
# GET DEFAULT PRIORS #
######################
# Change point
get_default_prior_cp = function(ST, i, default_brutto) {
  if (i < 2)
    stop("Only i >= 2 allowed.")

  # An absolute change point intercept
  if (i >= 2)
    #return(sprintf(priors[[family$family]]$cp, ST$cp_code_prior[i - 1]))
    return(truncate_prior_cp(ST, i, default_brutto$cp))
}


# Varying-by-group change point
get_default_prior_cp_group = function(ST, i) {
  if (i < 2)
    stop("Only i >= 2 allowed.")

  # Truncate between last change point and next change point, including their
  # varying effects, but keep in the observed range (MINX, MAXX).
  if (i == 2)
    trunc_from = paste0("MINX - ", ST$cp_code_prior[i])
  if (i > 2)
    trunc_from = paste0(ST$cp_code_prior[i-1], " - ", ST$cp_code_prior[i])
  if (i == nrow(ST))
    trunc_to = paste0("MAXX - ", ST$cp_code_prior[i])
  if (i < nrow(ST))
    trunc_to = paste0(ST$cp_code_prior[i + 1], " - ", ST$cp_code_prior[i])
  trunc = paste0("T(", trunc_from, ", ", trunc_to, ")")
  return(paste0("dnorm(0, ", ST$cp_sd[i], ") ", trunc))
}



#######################
# USER-DEFINED PRIORS #
#######################
truncate_prior_cp = function(ST, i, prior_str) {
  # Helper: Current segment.
  S = ST[i,]

  # User provided prior. Truncate it to ensure correct order
  is_bounded = stringr::str_detect(prior_str, "dunif|dirichlet")
  is_truncated = stringr::str_detect(prior_str, "T\\(")
  is_fixed = is.numeric(prior_str)

  # OK, we need to add truncation ourselves.
  if (!is_bounded & !is_truncated & !is_fixed) {
    if (S$cp_int_rel != 0) {
      # Relative: be positive (greater than former cp) and within observed range
      return(paste0(prior_str, " T(0, MAXX - ", ST$cp_code_prior[i - 1], ")"))
    }
    else {
      # Absolute: be greater than the former change point and within observed range
      return(paste0(prior_str, " T(", ST$cp_code_prior[i - 1], ", MAXX)"))
    }
  } else {
    # Return unaltered
    return(prior_str)
  }
}



####################
# ALL TOGETHER NOW #
####################

#' Get priors for all parameters in a segment table.
#'
#' Starts by finding all default priors. Then replace them with user priors.
#' User priors for change points are truncated appropriately using
#' `truncate_prior_cp``, if not done manually by the user already.
#'
#' @aliases get_prior
#' @keywords internal
#' @inheritParams mcp
#' @param ST Tibble. A segment table returned by `get_segment_table`.
#' @param prior A list of user-defined priors. Will overwrite the relevant
#'   default priors.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A named list of strings. The names correspond to the parameter names
#'   and the strings are the JAGS code for the prior (before converting SD to
#'   precision).
#'

get_prior = function(ST, family, prior = list()) {
  # Get ready to populate this list
  default_this = list()
  brutto_name = paste0(family$family, "_", family$link)
  if (!exists(brutto_name, where = priors))
    stop("mcp does not currently support `family = ", family$family, "(link = '", family$link, "')`. Please raise an issue on GitHub if you need this.")
  default_brutto = priors[[brutto_name]]  # Get the full set of priors belong to this family-link combo

  # Add model-specific paramters
  for (i in seq_len(nrow(ST))) {
    # Helper: Current segment.
    S = ST[i, ]

    # Intercept
    if (!is.na(S$ct_int))
      default_this[[S$ct_int[[1]]$name]] = default_brutto$ct_int

    # Each slope
    if (!is.na(S$ct_slope)) {
      for (name in S$ct_slope[[1]]$name) {
        default_this[[name]] = default_brutto$ct_slope
      }
    }

    # Change point.
    if (i > 1) {
      if (nrow(ST) == 2) {
        default_this[[S$cp_name]] = default_brutto$cp_1
      } else {
        default_this[[S$cp_name]] = get_default_prior_cp(ST, i, default_brutto)
      }
    }

    # Change point varying effects
    if (!is.na(S$cp_sd)) {
      default_this[[S$cp_sd]] = default_brutto$sd
      default_this[[S$cp_group]] = get_default_prior_cp_group(ST, i)
    }

    # Sigma intercept
    if (!is.na(S$sigma_int)) {
      for (name in S$sigma_int[[1]]$name) {
        default_this[[name]] = default_brutto$sigma_int
      }
    }

    # Sigma slope
    if (!is.na(S$sigma_slope)) {
      for (name in S$sigma_slope[[1]]$name) {
        default_this[[name]] = default_brutto$sigma_slope
      }
    }

    # MA intercept
    for (order in seq_len(sum(!is.na(S$ma_int[[1]])))) {  # Number of entries in int
      for (name in S$ma_int[[1]][[order]]$name) {
        default_this[[name]] = default_brutto$arma_int
      }
    }

    # MA slope
    for (order in seq_len(length(S$ma_slope[[1]]))) {  # Number of entries in slope
      if (!all(is.na(S$ma_slope[[1]][[order]]) == TRUE)) {  # If this slope exists...
        for (name in S$ma_slope[[1]][[order]]$name) {
          default_this[[name]] = default_brutto$arma_slope
        }
      }
    }

    # AR intercept
    for (order in seq_len(sum(!is.na(S$ar_int[[1]])))) {  # Number of entries in int
      for (name in S$ar_int[[1]][[order]]$name) {
        default_this[[name]] = default_brutto$arma_int
      }
    }

    # AR slope
    for (order in seq_len(length(S$ar_slope[[1]]))) {  # Number of entries in slope
      if (!all(is.na(S$ar_slope[[1]][[order]]) == TRUE)) {  # If this slope exists...
        for (name in S$ar_slope[[1]][[order]]$name) {
          default_this[[name]] = default_brutto$arma_slope
        }
      }
    }

    # Truncate change point prior if supplied by user
    if (i > 1 & ST$cp_name[i] %in% names(prior)) {
      prior[[S$cp_name]] = truncate_prior_cp(ST, i, prior[[S$cp_name]])
    }
  }

  # A check
  name_matches = names(prior) %in% names(default_this)
  if (!all(name_matches))
    stop("Prior(s) were specified for the following parmameter name(s) that are not part of the model: '", paste0(names(prior)[!name_matches], collapse = "', '"), "'")

  # Replace default priors with user prior and return
  this_prior = utils::modifyList(default_this, prior)

  # Sort according to type
  i_cp = stringr::str_starts(names(this_prior), "cp_")
  i_ints = stringr::str_starts(names(this_prior), "int_")
  i_slopes = stringr::str_starts(names(this_prior), paste0(ST$x[1], "_"))
  i_sigma = stringr::str_starts(names(this_prior), "sigma_")
  i_ar = stringr::str_starts(names(this_prior), "ar[0-9]+")
  i_ma = stringr::str_starts(names(this_prior), "ma[0-9]+")

  this_prior = c(this_prior[i_cp], this_prior[i_ints], this_prior[i_slopes], this_prior[i_sigma], this_prior[i_ar], this_prior[i_ma])
  class(this_prior) = "mcpprior"
  return(this_prior)
}
