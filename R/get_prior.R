##########################
# LIST OF DEFAULT PRIORS #
##########################

# Default priors for NORMAL (dp_n)
dp_n = list(
  cp_1 = "dunif(MINX, MAXX)",
  cp = "dunif(%s, MAXX)",
  cp_rel = "dunif(0, MAXX - %s)",
  slope = "dt(0, SDY / (MAXX - MINX) * 3, 1)",
  int = "dt(0, SDY * 3, 1)",
  sigma = "dnorm(0, SDY) T(0, )"
)

# Default priors for BERNOULLI (dp_b)
# PLACEHOLDER. NOT IN USE YET
dp_b = list(
  # empty
)



######################
# GET DEFAULT PRIORS #
######################

get_default_prior_cp = function(ST, i) {
  if(i < 2)
    stop("Only i >= 2 allowed.")

  # First change point
  if(i == 2)
    return(dp_n$cp_1)

  # A relative change point intercept
  if(i > 2 & ST$cp_int_rel[i] != 0)
    return(sprintf(dp_n$cp_rel, ST$cp_code[i - 1]))

  # An absolute change point intercept
  if(i > 2 & ST$cp_int_rel[i] == 0)
    return(sprintf(dp_n$cp, ST$cp_code[i - 1]))
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
  if(!is_dunif & !is_truncated & !is_fixed) {
    if(S$cp_int_rel != 0) {
      # Relative: be positive (greater than former cp) and within observed range
      return(paste0(prior_str, " T(0, MAXX - ", ST$cp_code[i-1], ")"))
    }
    else {
      # Absolute: be greater than the former change point and within observed range
      return(paste0(prior_str, " T(", ST$cp_code[i-1], ", MAXX)"))
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
#' Inserts default priors where none are specified by the user.
#'
#' @alias get_prior
#' @param ST Tibble. A segment table returned by \code{get_segment_table}.
#' @inheritParams mcp
#' @importFrom utils modifyList
#'
# Extracts default priors for model and replace them with user-priors when provided

get_prior = function(ST, prior) {
  # Get default priors
  default_prior = list(sigma = dp_n$sigma)  # Populate this list
  for(i in 1:nrow(ST)) {
    # Helper: Current segment.
    S = ST[i, ]

    # Intercept
    if(!is.na(S$int_name))
      default_prior[[S$int_name]] = dp_n$int

    # Slope
    if(!is.na(S$slope_name))
      default_prior[[S$slope_name]] = dp_n$slope

    # Change point
    if(i > 1)
      default_prior[[S$cp_name]] = get_default_prior_cp(ST, i)

    # Truncate change point prior if supplied by user
    if(i > 1 & ST$cp_name[i] %in% names(prior)) {
      prior[[S$cp_name]] = truncate_prior_cp(ST, i, prior[[S$cp_name]])
    }
  }

  # Replace default priors with user prior and return
  modifyList(default_prior, prior)
}
