##########################
# LIST OF DEFAULT PRIORS #
##########################

# Generic default priors that applies to many families
cp_prior = list(
  cp_1 = "dunif(MINX, MAXX)",
  cp = "dunif(%s, MAXX)",
  cp_rel = "dunif(0, MAXX - %s)",
  sd = "dnorm(0, (MAXX - MINX) / 2) T(0, )"
)

sigma_prior = list(
  sigma_int = "dnorm(0, SDY) T(0, )",
  sigma_slope = "dt(0, SDY / (MAXX - MINX), 3)"
)

arma_prior = list(
  arma_int = "dunif(-1, 1)",
  arma_slope = "dnorm(0, 1)"
)


# Per-family priors, mixing in the generic priors
priors = list(
  gaussian = c(cp_prior, sigma_prior, arma_prior, list(
    ct_slope = "dt(0, SDY / (MAXX - MINX), 3)",
    ct_int = "dt(0, 3 * SDY, 3)"
  )),

  # Identical priors for binomial and bernoulli.
  # A logit of +/- 5 is quite extreme. Very compatible with 3
  binomial = c(cp_prior, list(
    ct_slope = "dnorm(0, 3 / (MAXX - MINX))",
    ct_int = "dnorm(0, 3)"
  )),

  bernoulli = c(cp_prior, list(
    ct_slope = "dnorm(0, 3 / (MAXX - MINX))",
    ct_int = "dnorm(0, 3)"
  )),

  poisson = c(cp_prior, arma_prior, list(
    ct_slope = "dnorm(0, 10)",
    ct_int = "dnorm(0, 10)"
  ))
)



######################
# GET DEFAULT PRIORS #
######################
# Change point
get_default_prior_cp = function(ST, i, family) {
  if (i < 2)
    stop("Only i >= 2 allowed.")

  # First change point
  if (i == 2)
    return(priors[[family$family]]$cp_1)

  # A relative change point intercept
  if (i > 2 & ST$cp_int_rel[i] != 0)
    return(sprintf(priors[[family$family]]$cp_rel, ST$cp_code_prior[i - 1]))

  # An absolute change point intercept
  if (i > 2 & ST$cp_int_rel[i] == 0)
    return(sprintf(priors[[family$family]]$cp, ST$cp_code_prior[i - 1]))
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
  #trunc = paste0("T(", ifelse(i == 2, "MINX", ST$cp_code_prior[i-1]), ", ", ifelse(i == nrow(ST), "MAXX", ST$cp_code_prior[i+1]), ")")
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
  is_dunif = stringr::str_detect(prior_str, "dunif")
  is_truncated = stringr::str_detect(prior_str, "T\\(")
  is_fixed = is.numeric(prior_str)

  # OK, we need to add truncation ourselves.
  if (!is_dunif & !is_truncated & !is_fixed) {
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
#' @param ST Tibble. A segment table returned by `get_segment_table`.
#' @param prior A list of user-defined priors. Will overwrite the relevant
#'   default priors.
#' @inheritParams mcp
#'

get_prior = function(ST, family, prior = list()) {
  # Populate this list
  default_prior = list()

  # # Add MA priors
  # if (is.numeric(ar_order)) {
  #   for(i in seq_len(ar_order)) {
  #     default_prior[[paste0("phi", i)]] = priors[[family$family]]$phi
  #   }
  # }

  # Add model-specific paramters
  for (i in seq_len(nrow(ST))) {
    # Helper: Current segment.
    S = ST[i, ]

    # Intercept
    if (!is.na(S$ct_int))
      default_prior[[S$ct_int[[1]]$name]] = priors[[family$family]]$ct_int

    # Each slope
    if (!is.na(S$ct_slope)) {
      for (name in S$ct_slope[[1]]$name) {
        default_prior[[name]] = priors[[family$family]]$ct_slope
      }
    }

    # Change point
    if (i > 1)
      default_prior[[S$cp_name]] = get_default_prior_cp(ST, i, family)

    # Change point varying effects
    if (!is.na(S$cp_sd)) {
      default_prior[[S$cp_sd]] = priors[[family$family]]$sd
      default_prior[[S$cp_group]] = get_default_prior_cp_group(ST, i)
    }

    # Sigma intercept
    if (!is.na(S$sigma_int)) {
      for (name in S$sigma_int[[1]]$name) {
        default_prior[[name]] = priors[[family$family]]$sigma_int
      }
    }

    # Sigma slope
    if (!is.na(S$sigma_slope)) {
      for (name in S$sigma_slope[[1]]$name) {
        default_prior[[name]] = priors[[family$family]]$sigma_slope
      }
    }

    # MA
    # We know that S$ma_int and S$ma_slope have the same order
    n_ma = sum(!is.na(S$ma_int[[1]]))  # Number of ma intercepts
    for (order in seq_len(n_ma)) {
      # Intercept
      if (!all(is.na(S$ma_int[[1]][[order]]) == TRUE)) {  # If this intercept exists...
        for (name in S$ma_int[[1]][[order]]$name) {
          default_prior[[name]] = priors[[family$family]]$arma_int
        }
      }

      # Slope
      if (!all(is.na(S$ma_slope[[1]][[order]]) == TRUE)) {  # If this slope exists...
        for (name in S$ma_slope[[1]][[order]][["slope"]]$name) {
          default_prior[[name]] = priors[[family$family]]$arma_slope
        }
      }
    }

    # AR
    # We know that S$ar_int and S$ar_slope have the same order
    n_ar = sum(!is.na(S$ar_int[[1]]))  # Number of ma intercepts
    for (order in seq_len(n_ar)) {
      # Intercept
      if (!all(is.na(S$ar_int[[1]][[order]]) == TRUE)) {  # If this intercept exists...
        for (name in S$ar_int[[1]][[order]]$name) {
          default_prior[[name]] = priors[[family$family]]$arma_int
        }
      }

      # Slope
      if (!all(is.na(S$ar_slope[[1]][[order]]) == TRUE)) {  # If this slope exists...
        for (name in S$ar_slope[[1]][[order]][["slope"]]$name) {
          default_prior[[name]] = priors[[family$family]]$arma_slope
        }
      }
    }

    # Truncate change point prior if supplied by user
    if (i > 1 & ST$cp_name[i] %in% names(prior)) {
      prior[[S$cp_name]] = truncate_prior_cp(ST, i, prior[[S$cp_name]])
    }
  }

  # Replace default priors with user prior and return
  prior = utils::modifyList(default_prior, prior)

  # Sort by name and return
  prior[order(names(prior))]
}
