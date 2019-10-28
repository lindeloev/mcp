##########################
# LIST OF DEFAULT PRIORS #
##########################

# Generic default priors
common_prior = list(
  cp_1 = "dunif(MINX, MAXX)",
  cp = "dunif(%s, MAXX)",
  cp_rel = "dunif(0, MAXX - %s)"
)

# Per-family priors
priors = list(
  gaussian = c(common_prior, list(
    slope = "dnorm(0, 3 * SDY / (MAXX - MINX))",
    int = "dnorm(0, 3 * SDY)",
    sigma = "dnorm(0, 3 * SDY) T(0, )",
    sd = "dnorm(0, 3 * SDY) T(0, )"
  )),

  binomial = c(common_prior, list(
    # a logit of +/- 10 is quite extreme
    slope = "dnorm(0, 10 / (MAXX - MINX))",
    int = "dnorm(0, 10)",
    sd = "dnorm(0, 10) T(0, )"
  ))
)

# # Default priors for NORMAL (prior_gaussian)
# prior_gaussian = c(prior_common, list(
#   slope = "dnorm(0, SDY / (MAXX - MINX) * 3)",
#   int = "dnorm(0, SDY * 3)",
#   sigma = "dnorm(0, SDY * 3) T(0, )",
#   sd = "dnorm(0, SDY * 3) T(0, )"
# ))



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

  S = ST[i, ]
  trunc = paste0("T(", ifelse(i == 2, "MINX", ST$cp_code_prior[i-1]), ", ", ifelse(i == nrow(ST), "MAXX", ST$cp_code_prior[i+1]), ")")
  return(paste0("dnorm(0, ", S$cp_sd, ") ", trunc))
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
#' \code{truncate_prior_cp}, if not done manually by the user already.
#'
#' @aliases get_prior
#' @param ST Tibble. A segment table returned by \code{get_segment_table}.
#' @param prior A list of user-defined priors. Will overwrite the relevant
#'   default priors.
#' @inheritParams mcp
#'

get_prior = function(ST, family, prior = list()) {
  # Populate this list
  default_prior = list()

  # Add model-agnostic parameters
  if (family == "gaussian")
    default_prior[["sigma"]] = priors[[family]]$sigma

  # Add model-specific paramters
  for (i in seq_len(nrow(ST))) {
    # Helper: Current segment.
    S = ST[i, ]

    # Intercept
    if (!is.na(S$int_name))
      default_prior[[S$int_name]] = priors[[family]]$int

    # Slope
    if (!is.na(S$slope_name))
      default_prior[[S$slope_name]] = priors[[family]]$slope

    # Change point
    if (i > 1)
      default_prior[[S$cp_name]] = get_default_prior_cp(ST, i, family)

    # Change point varying effects
    if (!is.na(S$cp_sd)) {
      default_prior[[S$cp_sd]] = priors[[family]]$sd
      default_prior[[S$cp_group]] = get_default_prior_cp_group(ST, i)
    }

    # Truncate change point prior if supplied by user
    if (i > 1 & ST$cp_name[i] %in% names(prior)) {
      prior[[S$cp_name]] = truncate_prior_cp(ST, i, prior[[S$cp_name]])
    }
  }

  # Replace default priors with user prior and return
  utils::modifyList(default_prior, prior)
}
