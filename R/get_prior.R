##########################
# LIST OF DEFAULT PRIORS #
##########################

default_common_priors = tibble::tribble(
  ~family, ~link, ~dpar, ~par_type, ~prior,
  # Changepoint (common for all families)
  "changepoint", "identity", "cp_1", "Intercept", "dunif(MINX, MAXX)",  # If there is only one change point
  "changepoint", "identity", "cp", "Intercept", "dt(MINX, (MAXX - MINX) / N_CP, N_CP - 1)",  # For 2+ change points. Read the mcp article for more details on this "t-prior".
  "changepoint", "identity", "cp_sd", "dummy", "dnorm(0, 2 * (MAXX - MINX) / N_CP) T(0, )",  # Spread of varying effects

  # AR
  NA, "identity", "ar", "Intercept", "dbeta(0, 1)",
  NA, "identity", "ar", "dummy", "dt(0, 1, 3)",
  NA, "identity", "ar", "slope", "dt(0, 1 / (MAXX - MINX), 3)"
)

default_dpar_priors = tibble::tribble(
  ~family, ~link, ~dpar, ~par_type, ~prior,

  # Gaussian mu
  "gaussian", "identity", "mu", "Intercept", "dt(MEANLINKY, SDLINKY, 3)",
  "gaussian", "identity", "mu", "dummy", "dt(0, SDLINKY, 3)",
  "gaussian", "identity", "mu", "slope", "dt(0, N_CP * SDLINKY / (MAXX - MINX), 3)",

  # Gaussian sigma
  "gaussian", "identity", "sigma", "Intercept", "dt(0, SDLINKY, 3) T(0, )",  # Will always be <= observed SDLINKY
  "gaussian", "identity", "sigma", "dummy", "dt(0, SDLINKY, 3)",
  "gaussian", "identity", "sigma", "slope", "dt(0, N_CP * SDLINKY / (MAXX - MINX), 3)",

  # Binomial
  "binomial", "logit", "mu", "Intercept", "dt(0, 2.5, 3)",  # TO DO: center this on data
  "binomial", "logit", "mu", "dummy", "dt(0, 2.5, 3)",
  "binomial", "logit", "mu", "slope", "dt(0, N_CP * 2.5 / (MAXX - MINX), 3)",

  "binomial", "identity", "mu", "Intercept", "dbeta(0, 1)",  # TO DO: center this on data
  "binomial", "identity", "mu", "dummy", "dunif(-1, 1)",
  "binomial", "identity", "mu", "slope", "dt(0, N_CP / (MAXX - MINX), 3)",

  # Poisson
  "poisson", "log", "mu", "Intercept", "dt(0, 10)",  # TO DO: center this on data
  "poisson", "log", "mu", "dummy", "dt(0, 10)",
  "poisson", "log", "mu", "slope", "dt(0, 10)",

  "poisson", "identity", "mu", "Intercept", "dt(MEANLINKY, MEANLINKY, 3)",
  "poisson", "identity", "mu", "dummy", "dt(0, MEANLINKY, 3)",
  "poisson", "identity", "mu", "slope", "dt(0, N_CP * MEANLINKY / (MAXX - MINX), 3)",

  # Negbinomial shape
  "negbinomial", "identity", "shape", "Intercept", "dt(0, 10, 3) T(0, )",  # TO DO: center this on data
  "negbinomial", "identity", "shape", "dummy", "dt(0, 10, 3)",
  "negbinomial", "identity", "shape", "slope", "dt(0, N_CP * 10 / (MAXX - MINX), 3) T(0, )"
)

# Copies sections to other link functions and families
default_dpar_priors = dplyr::bind_rows(
  default_dpar_priors,

  # Gaussian
  default_dpar_priors %>% dplyr::filter(family == "gaussian", link == "identity", dpar == "mu") %>% dplyr::mutate(link = "log"),
  default_dpar_priors %>% dplyr::filter(family == "gaussian", link == "identity", dpar == "sigma") %>% dplyr::mutate(link = "log"),

  # Bernoulli/binomial
  default_dpar_priors %>% dplyr::filter(family == "binomial", link == "logit", dpar == "mu") %>% dplyr::mutate(link = "probit"),
  default_dpar_priors %>% dplyr::filter(family == "binomial", link == "logit", dpar == "mu") %>% dplyr::mutate(family = "bernoulli"),
  default_dpar_priors %>% dplyr::filter(family == "binomial", link == "logit", dpar == "mu") %>% dplyr::mutate(family = "bernoulli", link = "probit"),
  default_dpar_priors %>% dplyr::filter(family == "binomial", link == "identity", dpar == "mu") %>% dplyr::mutate(family = "bernoulli"),

  # Negative binomial
  default_dpar_priors %>% dplyr::filter(family == "poisson", link == "log", dpar == "mu") %>% dplyr::mutate(family = "negbinomial"),
  default_dpar_priors %>% dplyr::filter(family == "poisson", link == "identity", dpar == "mu") %>% dplyr::mutate(family = "negbinomial"),
  default_dpar_priors %>% dplyr::filter(family == "negbinomial", link == "identity", dpar == "mu") %>% dplyr::mutate(link = "log")
)


######################
# GET DEFAULT PRIORS #
######################
# Change point
get_default_prior_cp = function(ST, i, cp_prior) {
  assert_integer(i, lower = 2, len = 1)
  return(truncate_prior_cp(ST, i, cp_prior$cp))
}


# Varying-by-group change point
get_default_prior_cp_group = function(ST, i) {
  assert_integer(i, lower = 2, len = 1)

  # Truncate between last change point and next change point, including their
  # varying effects, but keep in the observed range (MINX, MAXX).
  if (i == 2)
    trunc_from = paste0("MINX - ", ST$cp_name[i])
  if (i > 2)
    trunc_from = paste0(ST$cp_name[i-1], " - ", ST$cp_name[i])
  if (i == nrow(ST))
    trunc_to = paste0("MAXX - ", ST$cp_name[i])
  if (i < nrow(ST))
    trunc_to = paste0(ST$cp_name[i + 1], " - ", ST$cp_name[i])
  trunc = paste0("T(", trunc_from, ", ", trunc_to, ")")
  return(paste0("dnorm(0, ", ST$cp_sd[i], ") ", trunc))
}



####################################
# USER-DEFINED CHANGE POINT PRIORS #
####################################
truncate_prior_cp = function(ST, i, prior_str) {
  assert_integer(i, lower = 2, len = 1)
  assert_types(prior_str, "character", len = 1)

  # Helper: Current segment.
  S = ST[i,]

  # User provided prior. Truncate it to ensure correct order
  is_bounded = stringr::str_detect(prior_str, "dunif|dirichlet")
  is_truncated = stringr::str_detect(prior_str, "T\\(")
  is_fixed = is.numeric(prior_str)

  # Absolute: be greater than the former change point and within observed range
  if (is_bounded == FALSE && is_truncated == FALSE && is_fixed == FALSE) {
    return(paste0(prior_str, " T(", ST$cp_name[i - 1], ", MAXX)"))
  } else {
    # Return unaltered
    return(prior_str)
  }
}



####################
# ALL TOGETHER NOW #
####################

#' Get priors for all parameters in the model
#'
#' Starts by finding all default priors. Then replace them with user priors.
#' User priors for change points are truncated appropriately using
#' `truncate_prior_cp``, if not done manually by the user already.
#'
#' @aliases get_prior
#' @keywords internal
#' @inheritParams mcp
#' @param ST Tibble. A segment table as returned by `get_segment_table`.
#' @param rhs_table Tibble as returned by `get_rhs()`.
#' @param family An `mcpfamily` object as returned by `mcpfamily()`.
#' @param prior A list of user-defined priors. Will overwrite the relevant
#'   default priors.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A named list of strings. The names correspond to the parameter names
#'   and the strings are the JAGS code for the prior (before converting SD to
#'   precision).
get_prior = function(ST, rhs_table, family, prior = list()) {
  assert_types(family, "mcpfamily")

  # Get ready to populate this list
  default_prior = list()

  # Priors for change points
  cp_prior = default_common_priors %>%
    dplyr::filter(family == "changepoint", link == "identity") %>%
    dplyr::select(dpar, prior) %>%
    tibble::deframe() %>%
    as.list()

  for (i in seq_len(nrow(ST))) {
    # Helper: Current segment.
    S = ST[i, ]

    # Change point
    if (i > 1) {
      if (nrow(ST) == 2) {
        default_prior[[S$cp_name]] = cp_prior$cp_1
      } else {
        default_prior[[S$cp_name]] = get_default_prior_cp(ST, i, cp_prior)
      }
    }

    # Change point varying effects
    if (!is.na(S$cp_sd)) {
      default_prior[[S$cp_sd]] = cp_prior$sd
      default_prior[[S$cp_group]] = get_default_prior_cp_group(ST, i)
    }

    # Truncate change point prior if supplied by user
    if (i > 1 && ST$cp_name[i] %in% names(prior)) {
      prior[[S$cp_name]] = truncate_prior_cp(ST, i, prior[[S$cp_name]])
    }
  }

  # Priors for RHS parameters
  rhs_prior = rhs_table %>%
    dplyr::left_join(family$default_prior, by = c("dpar", "par_type")) %>%
    dplyr::select(code_name, prior) %>%
    tibble::deframe() %>%
    as.list()

  if (any(is.na(rhs_prior)))
    stop_github("mcp could not find a default prior for ", and_collapse(names(rhs_prior[is.na(rhs_prior)])))

  # Merge
  default_prior = c(default_prior, rhs_prior)

  # A check
  name_matches = names(prior) %in% names(default_prior)
  if (any(name_matches == FALSE))
    stop("Prior(s) were specified for the following parmameter name(s) that are not part of the model: ", and_collapse(names(prior)[!name_matches]))

  # Replace default priors with user prior and return
  default_prior = utils::modifyList(default_prior, prior)
  return(default_prior)
}
