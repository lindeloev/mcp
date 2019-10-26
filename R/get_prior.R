##########################
# LIST OF DEFAULT PRIORS #
##########################

# Default priors for NORMAL (dp_n)
dp_n = list(
  cp_1 = "dunif(MINX, MAXX)",
  cp = "dunif(%s, MAXX)",
  cp_rel = "dunif(0, MAXX - %s)",
  slope = "dnorm(0, SDY / (MAXX - MINX) * 3)",
  int = "dnorm(0, SDY * 3)",
  sigma = "dnorm(0, SDY * 3) T(0, )",
  sd = "dnorm(0, SDY * 3) T(0, )"
)

# Default priors for BERNOULLI (dp_b)
# PLACEHOLDER. NOT IN USE YET
dp_b = list(
  # empty
)



######################
# GET DEFAULT PRIORS #
######################
# Change point
get_default_prior_cp = function(ST, i) {
  if (i < 2)
    stop("Only i >= 2 allowed.")

  # First change point
  if (i == 2)
    return(dp_n$cp_1)

  # A relative change point intercept
  if (i > 2 & ST$cp_int_rel[i] != 0)
    return(sprintf(dp_n$cp_rel, ST$cp_code_prior[i - 1]))

  # An absolute change point intercept
  if (i > 2 & ST$cp_int_rel[i] == 0)
    return(sprintf(dp_n$cp, ST$cp_code_prior[i - 1]))
}


# Varying-by-group change point
get_default_prior_cp_group = function(ST, i) {
  if(i < 2)
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

get_prior = function(ST, prior = list()) {
  # Get default priors
  default_prior = list(sigma = dp_n$sigma)  # Populate this list
  for (i in seq_len(nrow(ST))) {
    # Helper: Current segment.
    S = ST[i, ]

    # Intercept
    if (!is.na(S$int_name))
      default_prior[[S$int_name]] = dp_n$int

    # Slope
    if (!is.na(S$slope_name))
      default_prior[[S$slope_name]] = dp_n$slope

    # Change point
    if (i > 1)
      default_prior[[S$cp_name]] = get_default_prior_cp(ST, i)

    # Change point varying effects
    if (!is.na(S$cp_sd)) {
      default_prior[[S$cp_sd]] = dp_n$sd
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
