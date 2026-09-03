# Family specifications ----------------------------------------------------

#' Bernoulli Family for mcp
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
#' @examples
#' # Fit a binary-response model with a probit link
#' data = data.frame(time = 1:6, y = c(0, 0, 0, 1, 1, 1))
#' fit = mcp(list(y ~ 1), data, family = bernoulli(link = "probit"), par_x = "time", sample = FALSE)
#' mcp_pars(fit)  # Show the parameters of the fitted Bernoulli model
bernoulli = function(link = "logit") {
  link = rlang::arg_match0(link, c("identity", "logit", "probit"))

  family = stats::binomial(link = link)
  family$family = "bernoulli"
  mcpfamily(family)
}


#' Negative Binomial for mcp
#'
#' Parameterized as `mu` (the conditional mean) and `shape` (the same quantity
#' as `size` in `rnbinom()`). Thus `Var(y) = mu + mu^2 / shape`, which approaches
#' the Poisson variance as `shape` approaches infinity.
#'
#' @aliases negbinomial
#' @param link Link function for `mu`.
#' @param link_shape Link function for `shape`.
#' @details `shape(1)` is added implicitly and is constant across segments unless
#'   a `shape()` formula is supplied. For example, `y ~ 1 + x + shape(1 + x)`
#'   models both the mean and shape. Regression coefficients for both dpars are
#'   on their link scales.
#' @export
#' @examples
#' # Fit an overdispersed count model with the default log links
#' data = data.frame(time = 1:6, count = c(1, 2, 8, 3, 12, 5))
#' fit = mcp(list(count ~ 1), data, family = negbinomial(), par_x = "time", sample = FALSE)
#' mcp_pars(fit)  # Show the mean and shape parameters
negbinomial = function(link = "log", link_shape = "log") {
  link = rlang::arg_match0(link, "log")
  link_shape = rlang::arg_match0(link_shape, "log")

  family = list(family = "negbinomial", link = link, link_shape = link_shape)
  class(family) = "family"
  mcpfamily(family)
}


#' Create or Test Objects of Class "mcpfamily"
#'
#' \lifecycle{experimental}
#'
#' Converts standard R family objects into `mcpfamily` objects used internally by
#' `mcp`. Supported family and link combinations include:
#' * `gaussian(link = "identity")` or `gaussian(link = "log")`
#' * `binomial(link = "logit")`, `binomial(link = "probit")`, or `binomial(link = "identity")`
#' * `bernoulli(link = "logit")`, `bernoulli(link = "probit")`, or `bernoulli(link = "identity")`
#' * `poisson(link = "log")` or `poisson(link = "identity")`
#' * `negbinomial(link = "log", link_shape = "log")`
#'
#' Note: `mcpfamily` objects are shipped with mcp - there is not (yet) support for user-supplied
#' families.
#'
#' @aliases mcpfamily
#' @param x A family object, e.g., `binomial(link = "identity")`.
#' @seealso \code{\link{family}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' # mcp() converts supported standard family objects automatically
#' data = data.frame(time = 1:6, y = exp(seq(0, 1, length.out = 6)))
#' my_normal = mcpfamily(stats::gaussian(link = "log"))
#' fit = mcp(list(y ~ 1), data, family = my_normal, par_x = "time", sample = FALSE)
#' family(fit)  # Show the mcp family retained in the fit
#'
#' # The converted object can also be inspected directly
#' mcpfamily(stats::binomial())$dpars  # Show its distributional parameters
mcpfamily = function(x) {
  checkmate::assert_true(is.family(x), .var.name = "x")

  family = switch(
    x$family,
    gaussian = mcpfamily_gaussian(x),
    binomial = mcpfamily_binomial(x),
    bernoulli = mcpfamily_bernoulli(x),
    poisson = mcpfamily_poisson(x),
    negbinomial = mcpfamily_negbinomial(x),
    if (inherits(x, "mcpfamily")) x else new_mcpfamily(x)
  )

  family
}


new_mcpfamily = function(family, dpar_specs = family$dpar_specs,
                         default_prior = family$default_prior,
                         response = family$response, 
                         r = family$r,
                         backends = family$backends, 
                         garma = family$garma) {
  if (is.null(dpar_specs) || is.null(default_prior) || is.null(response) ||
      is.null(r) || is.null(backends$jags)) {
    stop(
      "family = ", family$family, "() does not have a complete internal mcp ",
      "family specification."
    )
  }

  # Default values and funcs for the dependent variable (the response)
  response_defaults = list(
    auxiliary = list(weights = list(
      required = FALSE,  # Is this required in the formula? FALSE = can be omitted.
      operations = "log_lik"  # For which operations is this auxiliary column needed?
    )),
    validate = function(y, data, response_columns) invisible(TRUE),  # Default to no validation
    observed = function(y, data, rate) y,
    probability = function(rate) FALSE,
    point_size = "weights"
  )

  # Fill in missing top-level response elements
  missing_response = setdiff(names(response_defaults), names(response))
  response[missing_response] = response_defaults[missing_response]

  # Merge missing auxiliary elements (e.g., weights) from response_defaults
  missing_aux = setdiff(names(response_defaults$auxiliary), names(response$auxiliary))
  response$auxiliary[missing_aux] = response_defaults$auxiliary[missing_aux]

  if (is.null(response$point_size))
    response$point_size = response_defaults$point_size

  # Each family will define a validate(); this applies that.
  family_validate = response$validate
  response$validate = function(y, data, response_columns) {
    # Check weights if present
    if (!is.null(data$weights)) {
      if (!is.numeric(data$weights) || anyNA(data$weights) || any(data$weights <= 0))
        stop("All weights must be numeric and greater than zero.")
    }

    # Apply the family-specific validation
    family_validate(y, data, response_columns)
  }

  # R-side log-likelihood. Weights are defined as applying directly here
  family_log_lik = r$log_lik
  r$log_lik = function(y, dpars, data) {
    weights = if (is.null(data$weights)) 1 else data$weights
    weights * family_log_lik(y, dpars, data)
  }

  family$dpar_specs = dpar_specs
  family$default_prior = normalize_family_default_priors(default_prior)
  family$response = response
  family$r = r
  family$backends = backends
  family$garma = garma
  family = add_dpar_specs(family)

  family$linkfun_str = get_link_str(family$link)
  family$linkinv_str = get_link_str(family$link, inverse = TRUE)
  if (!rlang::has_name(family, "linkfun"))
    family$linkfun = get_link_function(family$link)
  if (!rlang::has_name(family, "linkinv"))
    family$linkinv = get_link_function(family$link, inverse = TRUE)

  class(family) = c("mcpfamily", "family")
  family
}


validate_count_mean = function(mu) {
  if (any(mu > 2146275819, na.rm = TRUE))
    stop("Modelled extremely large count mean (> 2146275819).")
  invisible(TRUE)
}


normalize_family_default_priors = function(defaults) {
  checkmate::assert_data_frame(defaults)
  assert_data_cols(defaults, c("dpar", "par_type", "prior"))
  defaults = tibble::as_tibble(defaults)
  if (!"group_sd_prior" %in% names(defaults))
    defaults$group_sd_prior = NA_character_
  if (!"description" %in% names(defaults))
    defaults$description = NA_character_
  if (!"condition" %in% names(defaults))
    defaults$condition = "always"
  defaults[, c("dpar", "par_type", "prior", "group_sd_prior", "description", "condition")]
}


mcpfamily_gaussian = function(family) {
  if (family$link %notin% c("identity", "log"))
    stop("mcp has no default priors for gaussian(link = \"", family$link, "\") so it's likely not supported.")

  response_location = "round(median(.y), 1)"
  response_scale = "max(2.5, round(mad(.y), 1))"
  if (family$link == "identity") {
    mu_location = response_location
    mu_scale = response_scale
  } else {
    mu_location = "log_response_location(.y)"
    mu_scale = "log_response_scale(.y)"
  }

  default_prior = tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
    "mu", "Intercept", paste0("dt(", mu_location, ", ", mu_scale, ", 3)"), paste0("dt(0, ", mu_scale, ", 3) T(0, )"), "Robustly centered mean intercept with a minimum scale of 2.5", "always",
    "mu", "dummy", paste0("dt(0, ", mu_scale, ", 3)"), paste0("dt(0, ", mu_scale, ", 3) T(0, )"), "Regularizing mean contrast on the link scale", "always",
    "mu", "slope", paste0("dt(0, ", mu_scale, " / predictor_scale(), 3)"), paste0("dt(0, ", mu_scale, " / predictor_scale(), 3) T(0, )"), "Regularizing mean coefficient scaled to a reference predictor change", "always",
    "sigma", "Intercept", paste0("dt(0, ", response_scale, ", 3) T(0, )"), paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Positive residual SD calibrated on the response scale", "constant",
    "sigma", "Intercept", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Weakly regularizing modeled log-SD intercept", "modeled",
    "sigma", "dummy", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Regularizing log-SD contrast", "always",
    "sigma", "slope", "dt(0, 2.5 / predictor_scale(), 3)", "dt(0, 2.5 / predictor_scale(), 3) T(0, )", "Regularizing log-SD coefficient scaled to a reference predictor change", "always"
  )

  response = list(
    auxiliary = list(weights = list(
      required = FALSE,
      operations = "log_lik"
    )),
    point_size = "weights"
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    log_lik = function(y, dpars, data) stats::dnorm(y, dpars$mu, dpars$sigma, log = TRUE),
    rng = function(n, dpars, data, rate = FALSE) {
      stats::rnorm(n, dpars$mu, dpars$sigma)
    },
    cdf = function(q, dpars, data, rate = FALSE) {
      stats::pnorm(q, dpars$mu, dpars$sigma)
    }
  )
  jags = list(
    likelihood = function(context) {
      weights = context$aux("weights", "1")

      if (identical(weights, "1")) {
        # No weight
        paste0(context$y, " ~ dnorm(", context$dpar("mu"), ", 1 / ", context$dpar("sigma"), "^2)")
      } else {
        # The weighted normal contributes sigma^-1 from its normalizing term.
        # Observing zero from dexp(sigma^(1-w)) contributes sigma^(1-w), so
        # their product is Normal(y | mu, sigma)^w up to a weight-only constant.
        c(
          "# Gaussian likelihood raised to the observation weight",
          paste0("likelihood_weight_[i_] = 1 + response_observed_[i_] * (", weights, " - 1)  # Ensures weight 1 if missing data"),
          paste0(context$y, " ~ dnorm(", context$dpar("mu"), ", likelihood_weight_[i_] / ", context$dpar("sigma"), "^2)"),
          paste0("likelihood_zero_[i_] ~ dexp(pow(", context$dpar("sigma"), ", 1 - likelihood_weight_[i_]))")
        )
      }
    }
  )
  garma = if (family$link == "identity") {
    list(
      observed_r = function(y, data, boundary) y,
      observed_jags = function(context) context$y,
      generate_message = "Generating residuals for AR(N) model since the response column/argument was not provided."
    )
  }

  new_mcpfamily(
    family,
    dpar_specs = dplyr::bind_rows(
      new_dpar_spec("mu", family$link),
      new_dpar_spec("sigma", "identity", implicit = TRUE, link_modeled = "log", lower = 0.001)
    ),
    default_prior = default_prior,
    response = response,
    r = r,
    backends = list(jags = jags),
    garma = garma
  )
}


mcpfamily_binomial = function(family) {
  if (family$link == "identity") {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dbeta(1, 1)", "dt(0, 1, 3) T(0, )", "Uniform success-probability intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "dt(0, 1, 3) T(0, )", "Success-probability difference between levels", "always",
      "mu", "slope", "dt(0, 1 / predictor_scale(), 3)", "dt(0, 1 / predictor_scale(), 3) T(0, )", "Regularizing success-probability coefficient scaled to a reference predictor change", "always"
    )
  } else if (family$link %in% c("logit", "probit")) {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dt(0, 1.5, 3)", "dt(0, 1.5, 3) T(0, )", "Weakly regularizing link-scale intercept", "always",
      "mu", "dummy", "dt(0, 1.5, 3)", "dt(0, 1.5, 3) T(0, )", "Weakly regularizing categorical contrast on the link scale", "always",
      "mu", "slope", "dt(0, 1.5 / predictor_scale(), 3)", "dt(0, 1.5 / predictor_scale(), 3) T(0, )", "Weakly regularizing link-scale coefficient scaled to a reference predictor change", "always"
    )
  } else {
    stop("mcp has no default priors for binomial(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    auxiliary = list(trials = list(
      required = TRUE,
      operations = c("epred", "log_lik", "rng", "garma")
    )),
    validate = function(y, data, response_columns) {
      checkmate::assert_integerish(y, lower = 0, .var.name = response_columns$y)
      checkmate::assert_integerish(data$trials, lower = 1, any.missing = FALSE, .var.name = response_columns$trials)
      invalid = which(!is.na(y) & !is.na(data$trials) & y > data$trials)
      if (length(invalid) > 0)
        stop(
          "For family = binomial(), responses in '", response_columns$y,
          "' cannot exceed trials in '", response_columns$trials,
          "'. Found invalid data in row(s): ", paste(invalid, collapse = ", "), "."
        )
      invisible(TRUE)
    },
    observed = function(y, data, rate) if (rate) y / data$trials else y,
    probability = function(rate) isTRUE(rate),
    is_discrete = TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) if (rate) dpars$mu else data$trials * dpars$mu,
    log_lik = function(y, dpars, data) stats::dbinom(y, data$trials, dpars$mu, log = TRUE),
    rng = function(n, dpars, data, rate = FALSE) {
      y = stats::rbinom(n, data$trials, dpars$mu)
      if (rate) y / data$trials else y
    },
    cdf = function(q, dpars, data, rate = FALSE) {
      if (rate) q = round(q * data$trials)
      stats::pbinom(q, data$trials, dpars$mu)
    }
  )
  jags = list(
    likelihood = function(context) {
      weights = context$aux("weights", "1")

      if (identical(weights, "1")) {
        paste0(context$y, " ~ dbin(", context$dpar("mu"), ", ", context$aux("trials"), ")")
      } else {
        c(
          "# Binomial likelihood raised to the observation weight",
          paste0("likelihood_weight_[i_] = 1 + response_observed_[i_] * (", weights, " - 1)  # Ensures weight 1 if missing data"),
          paste0(context$y, " ~ dbin(", context$dpar("mu"), ", ", context$aux("trials"), ")"),
          paste0("likelihood_zero_[i_] ~ dexp(pow(", context$dpar("mu"), ", (likelihood_weight_[i_] - 1) * ", context$y, ") * pow(1 - ", context$dpar("mu"), ", (likelihood_weight_[i_] - 1) * (", context$aux("trials"), " - ", context$y, ")))")
        )
      }
    }
  )
  garma = if (family$link == "logit") {
    list(
      observed_r = function(y, data, boundary) pmin(pmax(y, boundary), data$trials - boundary) / data$trials,
      observed_jags = function(context) paste0(
        "min(max(", context$y, ", ", context$boundary, "), ",
        context$aux("trials"), " - ", context$boundary, ") / ", context$aux("trials")
      )
    )
  }

  new_mcpfamily(
    family,
    dpar_specs = new_dpar_spec("mu", family$link),
    default_prior = default_prior,
    response = response,
    r = r,
    backends = list(jags = jags),
    garma = garma
  )
}


mcpfamily_bernoulli = function(family) {
  if (family$link == "identity") {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dbeta(1, 1)", "dt(0, 1, 3) T(0, )", "Uniform P(y = TRUE) intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "dt(0, 1, 3) T(0, )", "P(y = TRUE) difference between levels", "always",
      "mu", "slope", "dt(0, 1 / predictor_scale(), 3)", "dt(0, 1 / predictor_scale(), 3) T(0, )", "Regularizing P(y = TRUE) coefficient scaled to a reference predictor change", "always"
    )
  } else if (family$link %in% c("logit", "probit")) {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dt(0, 1.5, 3)", "dt(0, 1.5, 3) T(0, )", "Weakly regularizing link-scale intercept", "always",
      "mu", "dummy", "dt(0, 1.5, 3)", "dt(0, 1.5, 3) T(0, )", "Weakly regularizing categorical contrast on the link scale", "always",
      "mu", "slope", "dt(0, 1.5 / predictor_scale(), 3)", "dt(0, 1.5 / predictor_scale(), 3) T(0, )", "Weakly regularizing link-scale coefficient scaled to a reference predictor change", "always"
    )
  } else {
    stop("mcp has no default priors for bernoulli(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    validate = function(y, data, response_columns) {
      invalid = which(!is.na(y) & y %notin% c(0, 1))
      if (length(invalid) > 0)
        stop("Only responses 0 and 1 are allowed for family = bernoulli() in column '", response_columns$y, "'")
      invisible(TRUE)
    },
    probability = function(rate) TRUE,
    is_discrete = TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    log_lik = function(y, dpars, data) stats::dbinom(y, 1, dpars$mu, log = TRUE),
    rng = function(n, dpars, data, rate = FALSE) stats::rbinom(n, 1, dpars$mu),
    cdf = function(q, dpars, data, rate = FALSE) stats::pbinom(q, 1, dpars$mu)
  )
  jags = list(
    likelihood = function(context) {
      weights = context$aux("weights", "1")

      if (identical(weights, "1")) {
        paste0(context$y, " ~ dbern(", context$dpar("mu"), ")")
      } else {
        c(
          "# Bernoulli likelihood raised to the observation weight",
          paste0("likelihood_weight_[i_] = 1 + response_observed_[i_] * (", weights, " - 1)  # Ensures weight 1 if missing data"),
          paste0(context$y, " ~ dbern(", context$dpar("mu"), ")"),
          paste0("likelihood_zero_[i_] ~ dexp(pow(", context$dpar("mu"), ", (likelihood_weight_[i_] - 1) * ", context$y, ") * pow(1 - ", context$dpar("mu"), ", (likelihood_weight_[i_] - 1) * (1 - ", context$y, ")))")
        )
      }
    }
  )

  new_mcpfamily(
    family,
    dpar_specs = new_dpar_spec("mu", family$link),
    default_prior = default_prior,
    response = response,
    r = r,
    backends = list(jags = jags),
    garma = NULL
  )
}


mcpfamily_poisson = function(family) {
  response_location = "round(median(.y), 1)"
  response_scale = "max(2.5, round(mad(.y), 1))"
  if (family$link == "identity") {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", paste0("dt(", response_location, ", ", response_scale, ", 3) T(0, )"), paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Positive count intercept calibrated on the response scale", "always",
      "mu", "dummy", paste0("dt(0, ", response_scale, ", 3)"), paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Count contrast calibrated on the response scale", "always",
      "mu", "slope", paste0("dt(0, ", response_scale, " / predictor_scale(), 3)"), paste0("dt(0, ", response_scale, " / predictor_scale(), 3) T(0, )"), "Count coefficient scaled to a reference predictor change", "always"
    )
  } else if (family$link == "log") {
    count_y = "log(pmax(.y, 0.1))"
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", paste0("dt(round(median(", count_y, "), 1), max(2.5, round(mad(", count_y, "), 1)), 3)"), "dt(0, 2.5, 3) T(0, )", "Robustly centered log-count intercept with a minimum scale of 2.5", "always",
      "mu", "dummy", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Regularizing categorical contrast on the log scale", "always",
      "mu", "slope", "dt(0, 2.5 / predictor_scale(), 3)", "dt(0, 2.5 / predictor_scale(), 3) T(0, )", "Regularizing log-count coefficient scaled to a reference predictor change", "always"
    )
  } else {
    stop("mcp has no default priors for poisson(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    validate = function(y, data, response_columns) {
      checkmate::assert_integerish(y, lower = 0, .var.name = response_columns$y)
      invisible(TRUE)
    },
    is_discrete = TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    log_lik = function(y, dpars, data) {
      validate_count_mean(dpars$mu)
      stats::dpois(y, dpars$mu, log = TRUE)
    },
    rng = function(n, dpars, data, rate = FALSE) {
      validate_count_mean(dpars$mu)
      stats::rpois(n, dpars$mu)
    },
    cdf = function(q, dpars, data, rate = FALSE) {
      validate_count_mean(dpars$mu)
      stats::ppois(q, dpars$mu)
    }
  )
  jags = list(
    likelihood = function(context) {
      weights = context$aux("weights", "1")

      if (identical(weights, "1")) {
        paste0(context$y, " ~ dpois(", context$dpar("mu"), ")")
      } else {
        c(
          "# Poisson likelihood raised to the observation weight",
          paste0("likelihood_weight_[i_] = 1 + response_observed_[i_] * (", weights, " - 1)  # Ensures weight 1 if missing data"),
          paste0(context$y, " ~ dpois(", context$dpar("mu"), ")"),
          paste0("likelihood_zero_[i_] ~ dexp(pow(", context$dpar("mu"), ", (likelihood_weight_[i_] - 1) * ", context$y, ") * exp((1 - likelihood_weight_[i_]) * ", context$dpar("mu"), "))")
        )
      }
    }
  )
  garma = if (family$link == "log") {
    list(
      observed_r = function(y, data, boundary) pmax(y, boundary),
      observed_jags = function(context) paste0("max(", context$y, ", ", context$boundary, ")")
    )
  }

  new_mcpfamily(
    family,
    dpar_specs = new_dpar_spec("mu", family$link),
    default_prior = default_prior,
    response = response,
    r = r,
    backends = list(jags = jags),
    garma = garma
  )
}


mcpfamily_negbinomial = function(family) {
  if (family$link != "log" || family$link_shape != "log")
    stop("Negative-binomial models currently require log links for both mu and shape.")

  count_y = "log(pmax(.y, 0.1))"
  default_prior = tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
    "mu", "Intercept", paste0("dt(round(median(", count_y, "), 1), max(2.5, round(mad(", count_y, "), 1)), 3)"), "dt(0, 2.5, 3) T(0, )", "Robustly centered log-count intercept with a minimum scale of 2.5", "always",
    "mu", "dummy", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Regularizing categorical count contrast on the log scale", "always",
    "mu", "slope", "dt(0, 2.5 / predictor_scale(), 3)", "dt(0, 2.5 / predictor_scale(), 3) T(0, )", "Regularizing log-count coefficient scaled to a reference predictor change", "always",
    "shape", "Intercept", "dloginvgamma(0.4, 0.3)", NA_character_, "Weakly regularizing positive overdispersion shape", "constant",
    "shape", "Intercept", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Weakly regularizing modeled log-shape intercept", "modeled",
    "shape", "dummy", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Regularizing shape contrast on the log scale", "always",
    "shape", "slope", "dt(0, 2.5 / predictor_scale(), 3)", "dt(0, 2.5 / predictor_scale(), 3) T(0, )", "Regularizing log-shape coefficient scaled to a reference predictor change", "always"
  )
  response = list(
    validate = function(y, data, response_columns) {
      checkmate::assert_integerish(y, lower = 0, .var.name = response_columns$y)
      invisible(TRUE)
    },
    is_discrete = TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    log_lik = function(y, dpars, data) {
      validate_count_mean(dpars$mu)
      stats::dnbinom(y, mu = dpars$mu, size = dpars$shape, log = TRUE)
    },
    rng = function(n, dpars, data, rate = FALSE) {
      validate_count_mean(dpars$mu)
      stats::rnbinom(n, mu = dpars$mu, size = dpars$shape)
    },
    cdf = function(q, dpars, data, rate = FALSE) {
      validate_count_mean(dpars$mu)
      stats::pnbinom(q, mu = dpars$mu, size = dpars$shape)
    }
  )
  jags = list(likelihood = function(context) {
    mu = context$dpar("mu")
    shape = context$dpar("shape")
    prob = context$local("nb_prob")
    weights = context$aux("weights", "1")

    if (identical(weights, "1")) {
      c(
        paste0(prob, " = ", shape, " / (", shape, " + ", mu, ")"),
        paste0(context$y, " ~ dnegbin(", prob, ", ", shape, ")")
      )
    } else {
      log_rate = context$local("nb_log_rate")
      c(
        paste0(prob, " = ", shape, " / (", shape, " + ", mu, ")"),
        paste0(context$y, " ~ dnegbin(", prob, ", ", shape, ")"),
        "# Negative binomial likelihood raised to the observation weight",
        paste0("likelihood_weight_[i_] = 1 + response_observed_[i_] * (", weights, " - 1)  # Ensures weight 1 if missing data"),
        paste0(log_rate, " = (likelihood_weight_[i_] - 1) * (loggam(", context$y, " + ", shape, ") - loggam(", shape, ") + ", shape, " * log(", shape, " / (", shape, " + ", mu, ")) + ", context$y, " * log(", mu, " / (", shape, " + ", mu, ")))"),
        paste0("likelihood_zero_[i_] ~ dexp(exp(", log_rate, "))")
      )
    }
  })
  garma = list(
    observed_r = function(y, data, boundary) pmax(y, boundary),
    observed_jags = function(context) paste0("max(", context$y, ", ", context$boundary, ")")
  )

  new_mcpfamily(
    family,
    dpar_specs = dplyr::bind_rows(
      new_dpar_spec("mu", family$link),
      new_dpar_spec("shape", family$link_shape, implicit = TRUE)
    ),
    default_prior = default_prior,
    response = response,
    r = r,
    backends = list(jags = jags),
    garma = garma
  )
}


# Describe a distributional parameter
new_dpar_spec = function(dpar, link, implicit = FALSE, lower = NA_real_,
                          link_modeled = link) {
  tibble::tibble(
    dpar = dpar,
    link = link,
    link_constant = link,
    link_modeled = link_modeled,
    modeled = FALSE,
    implicit = implicit,
    lower = lower
  )
}


known_dpar_wrappers = function() {
  c("sigma", "shape")
}


# Normalize and expose a family's distributional-parameter metadata.
add_dpar_specs = function(family) {
  family$dpar_specs = normalize_dpar_specs(family$dpar_specs)
  assert_dpar_specs(family$dpar_specs)
  family$dpars = family$dpar_specs$dpar
  family$links = stats::setNames(family$dpar_specs$link, family$dpar_specs$dpar)
  family
}


# Resolve each active dpar link from whether its formula is explicit.
resolve_dpar_specs = function(family, predictors, model = NULL) {
  checkmate::assert_true(is.mcpfamily(family), .var.name = "family")
  checkmate::assert_data_frame(predictors)
  assert_data_cols(predictors, c("dpar", "explicit"))
  family = add_dpar_specs(family)

  modeled_dpars = unique(predictors$dpar[predictors$explicit])
  if (!is.null(model)) {
    term_labels = unlist(lapply(model, function(segment) {
      rhs = get_rhs(segment)
      attr(stats::terms(rhs), "term.labels")
    }))
    declared_dpars = family$dpar_specs$dpar[vapply(
      family$dpar_specs$dpar,
      function(dpar) any(stringr::str_detect(term_labels, paste0("^", dpar, "\\("))),
      logical(1)
    )]
    modeled_dpars = union(modeled_dpars, declared_dpars)
  }

  family$dpar_specs$modeled = vapply(
    family$dpar_specs$dpar,
    function(dpar) dpar %in% modeled_dpars,
    logical(1)
  )
  family$dpar_specs$link = ifelse(
    family$dpar_specs$modeled,
    family$dpar_specs$link_modeled,
    family$dpar_specs$link_constant
  )

  add_dpar_specs(family)
}


# Add conditional-link fields missing from older dpar specifications.
normalize_dpar_specs = function(x) {
  checkmate::assert_data_frame(x)
  if ("link_constant" %notin% names(x))
    x$link_constant = x$link
  if ("link_modeled" %notin% names(x))
    x$link_modeled = x$link
  if ("modeled" %notin% names(x))
    x$modeled = FALSE
  x
}


# Validate distributional-parameter metadata and supported links.
assert_dpar_specs = function(x) {
  x = normalize_dpar_specs(x)
  required = c(
    "dpar", "link", "link_constant", "link_modeled", "modeled",
    "implicit", "lower"
  )
  assert_data_cols(x, required)

  if (!is.character(x$dpar) || anyNA(x$dpar) || any(x$dpar == ""))
    stop("`family$dpar_specs$dpar` must contain non-empty parameter names.")
  link_cols = c("link", "link_constant", "link_modeled")
  if (any(!vapply(x[link_cols], is.character, logical(1))) || anyNA(x[link_cols]))
    stop("Links in `family$dpar_specs` must be names without missing values.")
  if (anyDuplicated(x$dpar))
    stop("Each distributional parameter must occur exactly once in `family$dpar_specs`.")
  if ("mu" %notin% x$dpar)
    stop("`family$dpar_specs` must contain the response-mean parameter 'mu'.")
  if (any(x$dpar %in% c("epred", "ar", "ma")))
    stop("'epred', 'ar', and 'ma' are reserved and cannot be family distributional parameters.")

  supported_links = c("identity", "log", "logit", "probit")
  links = unique(unlist(x[link_cols], use.names = FALSE))
  if (any(links %notin% supported_links))
    stop("Unsupported dpar link(s): ", and_collapse(links[links %notin% supported_links]))
  if (!is.logical(x$modeled) || anyNA(x$modeled))
    stop("`family$dpar_specs$modeled` must be logical without missing values.")
  if (!is.logical(x$implicit) || anyNA(x$implicit))
    stop("`family$dpar_specs$implicit` must be logical without missing values.")
  if (!is.numeric(x$lower))
    stop("`family$dpar_specs$lower` must be numeric.")

  TRUE
}


# Return the single specification row for one distributional parameter.
get_dpar_spec = function(family, dpar) {
  spec = family$dpar_specs[family$dpar_specs$dpar == dpar, , drop = FALSE]
  if (nrow(spec) != 1)
    stop_github("Expected exactly one dpar specification for '", dpar, "'.")
  spec
}


# Find response-auxiliary data columns used by the requested operations.
get_family_aux_columns = function(family, segments, operations = NULL) {
  auxiliary = family$response$auxiliary
  if (!is.null(operations)) {
    used = vapply(auxiliary, function(x) any(x$operations %in% operations), logical(1))
    auxiliary = auxiliary[used]
  }
  aux_names = names(auxiliary)
  if (length(aux_names) == 0)
    return(stats::setNames(character(), character()))

  stats::setNames(vapply(aux_names, function(name) {
    if (name %notin% names(segments)) {
      return(NA_character_)
    } else {
      cols = unique(stats::na.omit(segments[[name]]))
      if (length(cols) > 1)
        stop("There should be exactly zero or one column used for ", name, "().")
      if (length(cols) == 0) NA_character_ else cols
    }
  }, character(1)), aux_names)
}


# Extract the response-auxiliary columns declared by a fitted family.
get_family_response_data = function(family, segments, data) {
  aux_columns = get_family_aux_columns(family, segments)
  out = lapply(aux_columns, function(column) {
    if (is.na(column) || column %notin% names(data)) NULL else data[[column]]
  })
  out[!vapply(out, is.null, logical(1))]
}


# Return the JAGS function name for a link or inverse link.
get_link_str = function(link, inverse = FALSE) {
  if (!inverse) {
    return(ifelse(link == "identity", "", link))
  } else {
    switch(link, logit = "ilogit", probit = "phi", log = "exp", identity = "")
  }
}


# Return the executable R function for a link or inverse link.
get_link_function = function(link, inverse = FALSE) {
  link_object = stats::make.link(link)
  if (inverse) link_object$linkinv else link_object$linkfun
}


#' @aliases format.mcpfamily
#' @describeIn mcpfamily Format an `mcpfamily` object.
#' @param ... Must be empty. Reserved for future use.
#' @export
format.mcpfamily = function(x, ...) {
  rlang::check_dots_empty()
  dpar_links = if (!is.null(x$dpar_specs)) {
    paste0(x$dpar_specs$dpar, " = ", x$dpar_specs$link, collapse = "; ")
  } else {
    paste0("mu = ", x$link)
  }
  paste0("Family: ", x$family, "\nLinks: ", dpar_links)
}


#' @aliases print.mcpfamily
#' @describeIn mcpfamily Print an `mcpfamily` object.
#' @export
print.mcpfamily = function(x, ...) {
  rlang::check_dots_empty()
  cat(format(x), "\n", sep = "")
  invisible(x)
}


# Check that an object implements the complete mcp family contract.
#' @aliases is.mcpfamily
#' @describeIn mcpfamily Checks whether x is an `mcpfamily`.
#' @export
is.mcpfamily = function(x) {
  if (!inherits(x, "mcpfamily"))
    return(FALSE)

  required = c("default_prior", "dpar_specs", "dpars", "links", "response", "r", "backends")
  if (any(required %notin% names(x)))
    return(FALSE)

  checkmate::assert_true(is.family(x), .var.name = "x")
  checkmate::assert_data_frame(x$default_prior)
  assert_data_cols(x$default_prior, c("dpar", "par_type", "prior"))
  assert_dpar_specs(x$dpar_specs)
  checkmate::assert_character(x$dpars)
  checkmate::assert_character(x$links)
  checkmate::assert_list(x$response$auxiliary)
  auxiliary_names = names(x$response$auxiliary)
  if (length(auxiliary_names) > 0) {
    checkmate::assert_character(auxiliary_names, any.missing = FALSE, min.len = 1)
    if (anyDuplicated(auxiliary_names))
      stop("Family response auxiliary names must be unique.")
    reserved = intersect(auxiliary_names, c("par_x", "response", "series"))
    if (length(reserved) > 0)
      stop(
        "Family response auxiliary names cannot use reserved column roles: ",
        and_collapse(reserved), "."
      )
  }
  checkmate::assert_function(x$response$validate)
  checkmate::assert_function(x$response$observed)
  checkmate::assert_function(x$response$probability)
  checkmate::assert_string(x$response$point_size, null.ok = TRUE)
  for (auxiliary in x$response$auxiliary) {
    checkmate::assert_flag(auxiliary$required)
    checkmate::assert_character(auxiliary$operations)
  }
  checkmate::assert_function(x$r$epred)
  checkmate::assert_function(x$r$log_lik)
  checkmate::assert_function(x$r$rng)
  checkmate::assert_function(x$backends$jags$likelihood)
  if (!is.null(x$garma)) {
    checkmate::assert_function(x$garma$observed_r)
    checkmate::assert_function(x$garma$observed_jags)
    checkmate::assert_string(x$garma$generate_message, null.ok = TRUE)
  }
  checkmate::assert_string(x$linkfun_str)
  checkmate::assert_string(x$linkinv_str)
  checkmate::assert_function(x$linkfun)
  checkmate::assert_function(x$linkinv)

  TRUE
}


logit = stats::binomial(link = "logit")$linkfun

ilogit = stats::binomial(link = "logit")$linkinv

probit = stats::binomial(link = "probit")$linkfun

phi = stats::binomial(link = "probit")$linkinv
