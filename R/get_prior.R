##########################
# LIST OF DEFAULT PRIORS #
##########################

# Generic default priors
common_prior = list(
  cp_1 = "dunif(MINX, MAXX)",
  cp = "dunif(%s, MAXX)",
  cp_rel = "dunif(0, MAXX - %s)",
  sd = "dnorm(0, (MAXX - MINX) / 2) T(0, )",
  phi = "dunif(-1, 1)"
)

# Per-family priors
priors = list(
  gaussian = c(common_prior, list(
    slope = "dt(0, SDY / (MAXX - MINX), 3)",
    int = "dt(0, 3 * SDY, 3)",
    sigma_int = "dnorm(0, SDY) T(0, )",
    sigma_slope = "dt(0, SDY / (MAXX - MINX), 3)"  # TO DO: this is probably very wide, leading to bad inits?
  )),

  # Identical priors for binomial and bernoulli.
  # A logit of +/- 5 is quite extreme. Very compatible with 3
  binomial = c(common_prior, list(
    slope = "dnorm(0, 3 / (MAXX - MINX))",
    int = "dnorm(0, 3)"
  )),

  bernoulli = c(common_prior, list(
    slope = "dnorm(0, 3 / (MAXX - MINX))",
    int = "dnorm(0, 3)"
  )),

  poisson = c(common_prior, list(
    slope = "dnorm(0, 10)",
    int = "dnorm(0, 10)"
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
    return(priors[[family]]$cp_1)

  # A relative change point intercept
  if (i > 2 & ST$cp_int_rel[i] != 0)
    return(sprintf(priors[[family]]$cp_rel, ST$cp_code_prior[i - 1]))

  # An absolute change point intercept
  if (i > 2 & ST$cp_int_rel[i] == 0)
    return(sprintf(priors[[family]]$cp, ST$cp_code_prior[i - 1]))
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

get_prior = function(ST, family, prior = list(), ar_order) {
  # Populate this list
  default_prior = list()

  # Add autocorrelation priors
  if (is.numeric(ar_order)) {
    for(i in seq_len(ar_order)) {
      default_prior[[paste0("phi_", i)]] = priors[[family]]$phi
    }
  }

  # Add sigma (may be deleted if specific variance change points are detected below)
  #if (family == "gaussian")
  #  default_prior[["sigma"]] = priors[[family]]$sigma

  # Add model-specific paramters
  for (i in seq_len(nrow(ST))) {
    # Helper: Current segment.
    S = ST[i, ]

    # Intercept
    if (!is.na(S$int))
      default_prior[[S$int[[1]]$name]] = priors[[family]]$int

    # Each slope
    if (!is.na(S$slope)) {
      for (name in S$slope[[1]]$name) {
        default_prior[[name]] = priors[[family]]$slope
      }
    }

    # Change point
    if (i > 1)
      default_prior[[S$cp_name]] = get_default_prior_cp(ST, i, family)

    # Change point varying effects
    if (!is.na(S$cp_sd)) {
      default_prior[[S$cp_sd]] = priors[[family]]$sd
      default_prior[[S$cp_group]] = get_default_prior_cp_group(ST, i)
    }

    # Sigma intercept
    if (!is.na(S$sigma_int)) {
      # Delete the default_named and insert segment-named instead.
      #default_prior[["sigma"]] = NULL
      #default_prior[["sigma_1"]] = priors[[family]]$sigma

      for (name in S$sigma_int[[1]]$name) {
        default_prior[[name]] = priors[[family]]$sigma_int
      }
    }

    # Sigma slope
    if (!is.na(S$sigma_slope)) {
      for (name in S$sigma_slope[[1]]$name) {
        default_prior[[name]] = priors[[family]]$sigma_slope
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
