# Family specifications ----------------------------------------------------

#' Bernoulli Family for mcp
#'
#' Binary responses (0 or 1). See [Binomial and Bernoulli](https://lindeloev.github.io/mcp/articles/binomial.html)
#' for worked examples.
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
#' @examples
#' # Fit a binary-response model with a probit link
#' data = data.frame(time = 1:6, y = c(0, 0, 0, 1, 1, 1))
#' fit = mcp(list(y ~ 1), data, family = bernoulli(link = "probit"), par_x = "time", sample = FALSE)
#' mcp_pars(fit)$parameter  # Show the parameters of the fitted Bernoulli model
bernoulli = function(link = "logit") {
  link = rlang::arg_match0(link, c("identity", "logit", "probit"))

  family = stats::binomial(link = link)
  family$family = "bernoulli"
  as_mcpfamily(family)
}


#' Negative Binomial for mcp
#'
#' Parameterized as `mu` (the conditional mean) and `shape` (the same quantity
#' as `size` in `rnbinom()`). Thus `Var(y) = mu + mu^2 / shape`, which approaches
#' the Poisson variance as `shape` approaches infinity.
#' [Read more](https://lindeloev.github.io/mcp/articles/poisson.html#negative-binomial-extension).
#'
#' @aliases negbinomial
#' @param link Link function for `mu`.
#' @param link_shape Link function for `shape` when it is modeled using `shape()`.
#' @details `shape` is required by the likelihood, so it includes an intercept in segment 1.
#'   Later segments that do not include `shape()` are joined. For example, `y ~ 1 + x + shape(1 + x)`
#'   models both the mean and shape. Regression coefficients for both dpars are
#'   on their link scales, e.g., `shape_Intercept_1` and `shape_x_1`.
#'
#'   Without `shape()`, `shape_1` is the shape itself with a `dinvgamma(0.4, 0.3)` prior,
#'   as in `brms`, much like the implicit `sigma_1` in `gaussian()`.
#' @export
#' @examples
#' # Fit an overdispersed count model with the default log links
#' data = data.frame(time = 1:6, count = c(1, 2, 8, 3, 12, 5))
#' fit = mcp(list(count ~ 1), data, family = negbinomial(), par_x = "time", sample = FALSE)
#' mcp_pars(fit)$parameter  # Show the mean and shape parameters
negbinomial = function(link = "log", link_shape = "log") {
  link = rlang::arg_match0(link, "log")
  link_shape = rlang::arg_match0(link_shape, "log")

  family = list(family = "negbinomial", link = link, link_shape = link_shape)
  class(family) = "family"
  as_mcpfamily(family)
}


# Convert a standard R family object (e.g., `binomial(link = "identity")`) into the
# `mcpfamily` used internally by mcp. Built-in families are (re)built from their family and
# link(s), so this is idempotent. Not exported until user-supplied families are supported.
as_mcpfamily = function(x) {
  checkmate::assert_true(is_family(x), .var.name = "x")

  switch(
    x$family,
    gaussian = mcpfamily_gaussian(x),
    binomial = mcpfamily_binomial(x),
    bernoulli = mcpfamily_bernoulli(x),
    poisson = mcpfamily_poisson(x),
    negbinomial = mcpfamily_negbinomial(x),
    if (is_mcpfamily(x)) assert_mcpfamily(x) else stop(
      "family = ", x$family, "() is not supported by mcp. Supported families are gaussian(), binomial(), ",
      "bernoulli(), poisson(), and negbinomial().", call. = FALSE
    )
  )
}


# Complete and validate an mcpfamily from the parts that each mcpfamily_*() supplies:
# - dpar_specs: one row per dpar (see new_dpar_spec()).
# - default_prior: one row per dpar and par_type. `condition` selects priors for a dpar
#   with ("modeled") or without ("constant") its own formula, e.g., sigma(1 + x).
# - response: auxiliary data columns (e.g., trials()) and how to validate and present the response.
# - r: R functions for predictions, log_lik(), etc. (see assert_mcpfamily()).
# - backends$jags: the likelihood in JAGS code and when to use the JAGS glm module.
# - garma: the observed response for AR/MA on the link scale, bounded by `threshold`. NULL if unsupported.
new_mcpfamily = function(family, dpar_specs, default_prior, response, r, backends, garma) {
  # Default values and funcs for the dependent variable (the response)
  response_defaults = list(
    auxiliary = list(weights = list(
      required = FALSE,  # Is this required in the formula? FALSE = need not be included.
      operations = "log_lik"  # For which operations is this auxiliary column needed?
    )),
    validate = function(y, data, data_columns) invisible(TRUE),  # Default to no validation
    observed = function(y, data, rate) y,
    probability = function(rate) FALSE,
    has_rate = FALSE,  # Does `rate = TRUE` apply, i.e., response per trial?
    is_discrete = FALSE  # Integer-valued responses? Used for predictive quantiles and weights
  )

  # Fill in missing top-level response elements
  missing_response = setdiff(names(response_defaults), names(response))
  response[missing_response] = response_defaults[missing_response]

  # Merge missing auxiliary elements (e.g., weights) from response_defaults
  missing_aux = setdiff(names(response_defaults$auxiliary), names(response$auxiliary))
  response$auxiliary[missing_aux] = response_defaults$auxiliary[missing_aux]

  family$dpar_specs = dpar_specs
  family$default_prior = standardize_family_default_priors(default_prior)
  family$response = response
  family$r = r
  family$backends = backends
  family$garma = garma

  family$linkfun = get_link_function(family$link)
  family$linkinv = get_link_function(family$link, inverse = TRUE)

  class(family) = c("mcpfamily", "family")
  assert_mcpfamily(family)
  family
}


# JAGS code raising a discrete likelihood to the observation weight using the Poisson zeros trick:
# P(Z = 0) from dpois(phi) contributes -phi to the log-likelihood, so phi = weight * -log(density).
# Missing responses (response_observed_ = 0) contribute nothing.
jags_weighted_likelihood = function(label, weights, neg_log_density) {
  c(
    paste0("# ", label, " likelihood raised to the observation weight"),
    paste0("likelihood_phi_[i_] = response_observed_[i_] * ", weights, " * max(0, ", neg_log_density, ")"),
    "likelihood_zero_[i_] ~ dpois(likelihood_phi_[i_])"
  )
}


validate_count_mean = function(mu) {
  if (any(is.infinite(mu)))
    stop("Modelled count mean overflowed to Inf. Consider narrower priors.")
  invisible(TRUE)
}


standardize_family_default_priors = function(defaults) {
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


# Default mu priors for log links, robustly centered on the log-response. `what` names the
# response in descriptions, e.g., "count" in "log-count intercept".
log_mu_priors = function(what) {
  log_y = "log(ifelse(.y <= 0, 0.1, .y)) - .offset"
  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
    "mu", "Intercept", paste0("dnorm(round(median(", log_y, "), 1), max(2.5, round(mad(", log_y, "), 1)))"), "dnorm(0, 2.5) T(0, )", paste0("Robustly centered log-", what, " intercept with a minimum scale of 2.5"), "always",
    "mu", "dummy", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5) T(0, )", paste0("Autoscaled log-", what, " contrast (rstanarm default)"), "always",
    "mu", "slope", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5 / predictor_scale()) T(0, )", paste0("Autoscaled log-", what, " coefficient (rstanarm default)"), "always"
  )
}


# Default mu priors for logit and probit links. The intercept prior is approximately (logit)
# or exactly (probit) uniform on the probability scale.
probability_mu_priors = function(link) {
  scale = if (link == "logit") "1.5" else "1"
  uniform = if (link == "logit") "approximately uniform" else "uniform"
  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
    "mu", "Intercept", paste0("dnorm(0, ", scale, ")"), paste0("dt(0, ", scale, ", 3) T(0, )"), paste0("Link-scale intercept, ", uniform, " on the probability scale"), "always",
    "mu", "dummy", paste0("dnorm(0, ", scale, " / predictor_scale())"), paste0("dt(0, ", scale, ", 3) T(0, )"), "Autoscaled link-scale contrast (rstanarm autoscaling)", "always",
    "mu", "slope", paste0("dnorm(0, ", scale, " / predictor_scale())"), paste0("dt(0, ", scale, " / predictor_scale(), 3) T(0, )"), "Autoscaled link-scale coefficient (rstanarm autoscaling)", "always"
  )
}


mcpfamily_gaussian = function(family) {
  if (family$link %notin% c("identity", "log"))
    stop("mcp has no default priors for gaussian(link = \"", family$link, "\") so it's likely not supported.")

  response_scale = "max(2.5, round(mad(.y), 1))"

  if (family$link == "identity") {
    mu_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dnorm(mean(.y - .offset), 2.5 * sd(.y - .offset))", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Mean intercept (rstanarm default)", "always",
      "mu", "dummy", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Autoscaled mean contrast (rstanarm default)", "always",
      "mu", "slope", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, " / predictor_scale(), 3) T(0, )"), "Autoscaled mean coefficient (rstanarm default)", "always"
    )
  } else {
    mu_prior = log_mu_priors("mean")
  }

  sigma_prior = tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
    "sigma", "Intercept", paste0("dt(0, ", response_scale, ", 3) T(0, )"), paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Positive residual SD calibrated on the response scale", "constant",
    "sigma", "Intercept", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Weakly regularizing modeled log-SD intercept", "modeled",
    "sigma", "dummy", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Regularizing log-SD contrast", "always",
    "sigma", "slope", "dt(0, 2.5 / predictor_scale(), 3)", "dt(0, 2.5 / predictor_scale(), 3) T(0, )", "Regularizing log-SD coefficient scaled by the predictor SD", "always"
  )

  default_prior = dplyr::bind_rows(mu_prior, sigma_prior)

  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    variance = function(dpars, data, rate = FALSE) dpars$sigma^2,
    log_lik = function(y, dpars, data) stats::dnorm(y, dpars$mu, dpars$sigma, log = TRUE),
    rng = function(n, dpars, data, rate = FALSE) {
      stats::rnorm(n, dpars$mu, dpars$sigma)
    },
    cdf = function(q, dpars, data) stats::pnorm(q, dpars$mu, dpars$sigma),
    quantile = function(p, dpars, data) stats::qnorm(p, dpars$mu, dpars$sigma)
  )
  jags = list(
    glm_module = "always",  # JAGS glm module: see run_jags()
    likelihood = function(context) {
      weights = context$aux("weights", "1")

      if (identical(weights, "1")) {
        # No weight
        paste0(context$y, " ~ dnorm(", context$dpar("mu"), ", 1 / ", context$dpar("sigma"), "^2)")
      } else {
        # Precision w / sigma^2 raises the kernel to weight w. The normalizing term needs
        # an extra -(w - 1) * log(sigma), added with the Poisson zeros trick (see
        # jags_weighted_likelihood()). The constant |w - 1| * 50 keeps phi >= 0 for
        # sigma within e^-50 to e^50 and cancels out in MCMC.
        c(
          "# Gaussian likelihood raised to the observation weight",
          paste0("likelihood_weight_[i_] = 1 + response_observed_[i_] * (", weights, " - 1)  # Ensures weight 1 if missing data"),
          paste0(context$y, " ~ dnorm(", context$dpar("mu"), ", likelihood_weight_[i_] / ", context$dpar("sigma"), "^2)"),
          paste0("likelihood_phi_[i_] = response_observed_[i_] * max(0, (likelihood_weight_[i_] - 1) * log(", context$dpar("sigma"), ") + abs(likelihood_weight_[i_] - 1) * 50)"),
          paste0("likelihood_zero_[i_] ~ dpois(likelihood_phi_[i_])")
        )
      }
    }
  )

  garma = if (family$link == "identity") {
    list(
      observed_r = function(y, data, threshold) y,
      observed_jags = function(context) context$y
    )
  }

  # Return
  new_mcpfamily(
    family,
    dpar_specs = dplyr::bind_rows(
      new_dpar_spec("mu", family$link),
      new_dpar_spec("sigma", "identity", implicit = TRUE, require_initial_term = TRUE, link_modeled = "log", lower = 0, floor = 1e-5)
    ),
    default_prior = default_prior,
    response = list(),
    r = r,
    backends = list(jags = jags),
    garma = garma
  )
}


mcpfamily_binomial = function(family) {
  if (family$link == "identity") {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dunif(0, 1)", "dt(0, 1, 3) T(0, )", "Uniform success-probability intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "dt(0, 1, 3) T(0, )", "Success-probability difference between levels", "always",
      "mu", "slope", "dnorm(0, 1 / predictor_scale())", "dt(0, 1 / predictor_scale(), 3) T(0, )", "Autoscaled success-probability coefficient", "always"
    )
  } else if (family$link %in% c("logit", "probit")) {
    default_prior = probability_mu_priors(family$link)
  } else {
    stop("mcp has no default priors for binomial(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    auxiliary = list(trials = list(
      required = TRUE,
      operations = c("epred", "log_lik", "rng", "garma")
    )),
    validate = function(y, data, data_columns) {
      checkmate::assert_integerish(y, lower = 0, .var.name = data_columns$response)
      checkmate::assert_integerish(data$trials, lower = 1, any.missing = FALSE, .var.name = data_columns$trials)
      invalid = which(!is.na(y) & !is.na(data$trials) & y > data$trials)
      if (length(invalid) > 0)
        stop(
          "For family = binomial(), responses in '", data_columns$response,
          "' cannot exceed trials in '", data_columns$trials,
          "'. Found invalid data in row(s): ", paste(invalid, collapse = ", "), "."
        )
      invisible(TRUE)
    },
    observed = function(y, data, rate) if (rate) y / data$trials else y,
    probability = function(rate) isTRUE(rate),
    has_rate = TRUE,
    is_discrete = TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) if (rate) dpars$mu else data$trials * dpars$mu,
    variance = function(dpars, data, rate = FALSE) if (rate) dpars$mu * (1 - dpars$mu) / data$trials else data$trials * dpars$mu * (1 - dpars$mu),
    log_lik = function(y, dpars, data) stats::dbinom(y, data$trials, dpars$mu, log = TRUE),
    rng = function(n, dpars, data, rate = FALSE) {
      y = stats::rbinom(n, data$trials, dpars$mu)
      if (rate) y / data$trials else y
    },
    cdf = function(q, dpars, data) stats::pbinom(q, data$trials, dpars$mu),
    quantile = function(p, dpars, data) stats::qbinom(p, data$trials, dpars$mu)
  )
  jags = list(
    glm_module = "group_effects",  # JAGS glm module: see run_jags()
    likelihood = function(context) {
      weights = context$aux("weights", "1")

      # Unweighted likelihood
      if (identical(weights, "1")) {
        paste0(context$y, " ~ dbin(", context$dpar("mu"), ", ", context$aux("trials"), ")")
      } else {
        trials = context$aux("trials")
        mu = context$dpar("mu")
        y = context$y
        jags_weighted_likelihood("Binomial", weights, paste0("loggam(", y, " + 1) + loggam(", trials, " - ", y, " + 1) - loggam(", trials, " + 1) - ", y, " * log(", mu, ") - (", trials, " - ", y, ") * log(1 - ", mu, ")"))
      }
    }
  )
  garma = if (family$link == "logit") {
    list(
      observed_r = function(y, data, threshold) pmin(pmax(y, threshold), data$trials - threshold) / data$trials,
      observed_jags = function(context) paste0(
        "min(max(", context$y, ", ", context$threshold, "), ",
        context$aux("trials"), " - ", context$threshold, ") / ", context$aux("trials")
      ),
      validate_threshold = function(threshold, data = list()) {
        min_trials = min(data$trials)
        max_threshold = min_trials / 2
        if (any(threshold <= 0 | threshold >= max_threshold))
          stop("`threshold` for family = binomial() must be strictly between 0 and half the smallest trial count (", max_threshold, ").")
      }
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
      "mu", "Intercept", "dunif(0, 1)", "dt(0, 1, 3) T(0, )", "Uniform P(y = TRUE) intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "dt(0, 1, 3) T(0, )", "P(y = TRUE) difference between levels", "always",
      "mu", "slope", "dnorm(0, 1 / predictor_scale())", "dt(0, 1 / predictor_scale(), 3) T(0, )", "Autoscaled P(y = TRUE) coefficient", "always"
    )
  } else if (family$link %in% c("logit", "probit")) {
    default_prior = probability_mu_priors(family$link)
  } else {
    stop("mcp has no default priors for bernoulli(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    validate = function(y, data, data_columns) {
      invalid = which(!is.na(y) & y %notin% c(0, 1))
      if (length(invalid) > 0)
        stop("Only responses 0 and 1 are allowed for family = bernoulli() in column '", data_columns$response, "'")
      invisible(TRUE)
    },
    probability = function(rate) TRUE,
    has_rate = TRUE,
    is_discrete = TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    variance = function(dpars, data, rate = FALSE) dpars$mu * (1 - dpars$mu),
    log_lik = function(y, dpars, data) stats::dbinom(y, 1, dpars$mu, log = TRUE),
    rng = function(n, dpars, data, rate = FALSE) stats::rbinom(n, 1, dpars$mu),
    cdf = function(q, dpars, data) stats::pbinom(q, 1, dpars$mu),
    quantile = function(p, dpars, data) stats::qbinom(p, 1, dpars$mu)
  )
  jags = list(
    glm_module = "group_effects",  # JAGS glm module: see run_jags()
    likelihood = function(context) {
      weights = context$aux("weights", "1")

      if (identical(weights, "1")) {
        paste0(context$y, " ~ dbern(", context$dpar("mu"), ")")
      } else {
        mu = context$dpar("mu")
        y = context$y
        jags_weighted_likelihood("Bernoulli", weights, paste0("- ", y, " * log(", mu, ") - (1 - ", y, ") * log(1 - ", mu, ")"))
      }
    }
  )

  garma = if (family$link == "logit") {
    list(
      observed_r = function(y, data, threshold) pmin(pmax(y, threshold), 1 - threshold),
      observed_jags = function(context) paste0(
        "min(max(", context$y, ", ", context$threshold, "), 1 - ", context$threshold, ")"
      ),
      validate_threshold = function(threshold, data = list()) {
        if (any(threshold <= 0 | threshold >= 0.5))
          stop("`threshold` for family = bernoulli() must be strictly between 0 and 0.5.")
      }
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


mcpfamily_poisson = function(family) {
  response_scale = "max(2.5, round(mad(.y), 1))"
  if (family$link == "identity") {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dnorm(mean(.y - .offset), 2.5 * sd(.y - .offset)) T(-min(.offset), )", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Count intercept on the response scale with a positive mean", "always",
      "mu", "dummy", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Autoscaled count contrast on the response scale", "always",
      "mu", "slope", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, " / predictor_scale(), 3) T(0, )"), "Autoscaled count coefficient on the response scale", "always"
    )
  } else if (family$link == "log") {
    default_prior = log_mu_priors("count")
  } else {
    stop("mcp has no default priors for poisson(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    validate = function(y, data, data_columns) {
      checkmate::assert_integerish(y, lower = 0, .var.name = data_columns$response)
      invisible(TRUE)
    },
    is_discrete = TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    variance = function(dpars, data, rate = FALSE) dpars$mu,
    log_lik = function(y, dpars, data) {
      validate_count_mean(dpars$mu)
      stats::dpois(y, dpars$mu, log = TRUE)
    },
    rng = function(n, dpars, data, rate = FALSE) {
      validate_count_mean(dpars$mu)
      stats::rpois(n, dpars$mu)
    },
    cdf = function(q, dpars, data) {
      validate_count_mean(dpars$mu)
      stats::ppois(q, dpars$mu)
    },
    quantile = function(p, dpars, data) {
      validate_count_mean(dpars$mu)
      stats::qpois(p, dpars$mu)
    }
  )
  jags = list(
    glm_module = "group_effects",  # JAGS glm module: see run_jags()
    likelihood = function(context) {
      weights = context$aux("weights", "1")

      # Unweighted likelihood
      if (identical(weights, "1")) {
        paste0(context$y, " ~ dpois(", context$dpar("mu"), ")")
      } else {
        mu = context$dpar("mu")
        y = context$y
        jags_weighted_likelihood("Poisson", weights, paste0(mu, " - ", y, " * log(", mu, ") + loggam(", y, " + 1)"))
      }
    }
  )
  garma = if (family$link == "log") {
    list(
      observed_r = function(y, data, threshold) pmax(y, threshold),
      observed_jags = function(context) paste0("max(", context$y, ", ", context$threshold, ")")
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

  shape_prior = tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
    "shape", "Intercept", "dinvgamma(0.4, 0.3)", NA_character_, "Weakly regularizing positive overdispersion shape (brms default)", "constant",
    "shape", "Intercept", "dnorm(0, 2.5)", "dnorm(0, 2.5) T(0, )", "Weakly regularizing modeled log-shape intercept", "modeled",
    "shape", "dummy", "dnorm(0, 2.5)", "dnorm(0, 2.5) T(0, )", "Regularizing shape contrast on the log scale", "always",
    "shape", "slope", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5 / predictor_scale()) T(0, )", "Regularizing log-shape coefficient scaled by the predictor SD", "always"
  )
  default_prior = dplyr::bind_rows(log_mu_priors("count"), shape_prior)
  response = list(
    validate = function(y, data, data_columns) {
      checkmate::assert_integerish(y, lower = 0, .var.name = data_columns$response)
      invisible(TRUE)
    },
    is_discrete = TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    variance = function(dpars, data, rate = FALSE) dpars$mu + dpars$mu^2 / dpars$shape,
    log_lik = function(y, dpars, data) {
      validate_count_mean(dpars$mu)
      stats::dnbinom(y, mu = dpars$mu, size = dpars$shape, log = TRUE)
    },
    rng = function(n, dpars, data, rate = FALSE) {
      validate_count_mean(dpars$mu)
      stats::rnbinom(n, mu = dpars$mu, size = dpars$shape)
    },
    cdf = function(q, dpars, data) {
      validate_count_mean(dpars$mu)
      stats::pnbinom(q, mu = dpars$mu, size = dpars$shape)
    },
    quantile = function(p, dpars, data) {
      validate_count_mean(dpars$mu)
      stats::qnbinom(p, mu = dpars$mu, size = dpars$shape)
    }
  )
  jags = list(
    glm_module = "never",  # JAGS glm module: see run_jags()
    likelihood = function(context) {
      mu = context$dpar("mu")
      shape = context$dpar("shape")
      prob = context$local("nb_prob")
      weights = context$aux("weights", "1")
      y = context$y
      prob_code = paste0(prob, " = ", shape, " / (", shape, " + ", mu, ")")

      # Unweighted likelihood
      if (identical(weights, "1"))
        return(c(prob_code, paste0(y, " ~ dnegbin(", prob, ", ", shape, ")")))

      c(prob_code, jags_weighted_likelihood("Negative binomial", weights, paste0("loggam(", shape, ") + loggam(", y, " + 1) - loggam(", y, " + ", shape, ") - ", shape, " * log(", prob, ") - ", y, " * log(1 - ", prob, ")")))
    }
  )
  garma = list(
    observed_r = function(y, data, threshold) pmax(y, threshold),
    observed_jags = function(context) paste0("max(", context$y, ", ", context$threshold, ")")
  )

  new_mcpfamily(
    family,
    dpar_specs = dplyr::bind_rows(
      new_dpar_spec("mu", family$link),
      new_dpar_spec("shape", "identity", implicit = TRUE, require_initial_term = TRUE, link_modeled = family$link_shape, lower = 0)
    ),
    default_prior = default_prior,
    response = response,
    r = r,
    backends = list(jags = jags),
    garma = garma
  )
}


# Describe a distributional parameter. `lower` is the exclusive lower bound of its
# support (e.g., 0 for sigma) and bounds constant priors. `floor` is an optional numerical
# floor on per-observation values (see dev/DECISIONS.md for sigma).
new_dpar_spec = function(dpar, link, implicit = FALSE, require_initial_term = FALSE, lower = NA_real_,
                          link_modeled = link, floor = NA_real_) {
  tibble::tibble(
    dpar = dpar,
    link = link,
    link_constant = link,
    link_modeled = link_modeled,
    modeled = FALSE,
    implicit = implicit,
    require_initial_term = require_initial_term,
    lower = lower,
    floor = floor
  )
}


# Distributional parameter wrapper names of the built-in families. A family's own dpars
# are always recognized; this list only gives informative errors when, e.g., `sigma()`
# is used with a family that lacks it.
known_dpar_wrappers = function() {
  c("mu", "sigma", "shape")
}


# Lower bounds of constant (unmodeled) dpars on their parameter scale, e.g., c(sigma = 0).
# Bounds that the link already enforces (e.g., log(0) = -Inf) are omitted.
get_constant_dpar_bounds = function(family) {
  dpar_specs = family$dpar_specs[!family$dpar_specs$modeled & !is.na(family$dpar_specs$lower), ]
  bounds = vapply(seq_len(nrow(dpar_specs)), function(i) get_link_function(dpar_specs$link[i])(dpar_specs$lower[i]), numeric(1))
  stats::setNames(bounds, dpar_specs$dpar)[is.finite(bounds)]
}


# Resolve each dpar link from whether the model includes a formula for it, e.g., sigma(1).
# - components: A list with `get_segment_components()` for each segment.
resolve_dpar_specs = function(family, components) {
  checkmate::assert_class(family, "mcpfamily")
  family$dpar_specs$modeled = family$dpar_specs$dpar %in% unlist(lapply(components, names))
  family$dpar_specs$link = ifelse(
    family$dpar_specs$modeled,
    family$dpar_specs$link_modeled,
    family$dpar_specs$link_constant
  )
  family
}


# Validate distributional-parameter metadata and supported links.
assert_dpar_specs = function(x) {
  checkmate::assert_data_frame(x)
  required = c(
    "dpar", "link", "link_constant", "link_modeled", "modeled",
    "implicit", "require_initial_term", "lower", "floor"
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

  supported_links = names(.mcplinks)
  links = unique(unlist(x[link_cols], use.names = FALSE))
  if (any(links %notin% supported_links))
    stop("Unsupported dpar link(s): ", and_collapse(links[links %notin% supported_links]))
  if (!is.logical(x$modeled) || anyNA(x$modeled))
    stop("`family$dpar_specs$modeled` must be logical without missing values.")
  if (!is.logical(x$implicit) || anyNA(x$implicit))
    stop("`family$dpar_specs$implicit` must be logical without missing values.")
  if (!is.logical(x$require_initial_term) || anyNA(x$require_initial_term))
    stop("`family$dpar_specs$require_initial_term` must be logical without missing values.")
  if (!is.numeric(x$lower))
    stop("`family$dpar_specs$lower` must be numeric.")
  if (!is.numeric(x$floor))
    stop("`family$dpar_specs$floor` must be numeric.")

  TRUE
}


# Return the single specification row for one distributional parameter.
get_dpar_spec = function(family, dpar) {
  dpar_spec = family$dpar_specs[family$dpar_specs$dpar == dpar, , drop = FALSE]
  if (nrow(dpar_spec) != 1)
    stop_github("Expected exactly one dpar specification for '", dpar, "'.")
  dpar_spec
}


# Find response-auxiliary data columns used by the requested operations.
get_family_aux_columns = function(family, data_columns, operations = NULL) {
  auxiliary = family$response$auxiliary
  if (!is.null(operations)) {
    used = vapply(auxiliary, function(x) any(x$operations %in% operations), logical(1))
    auxiliary = auxiliary[used]
  }
  aux_names = names(auxiliary)
  stats::setNames(vapply(aux_names, function(name) {
    if (is.null(data_columns[[name]])) NA_character_ else data_columns[[name]]
  }, character(1)), aux_names)
}


# Extract the response-auxiliary columns defined by a fitted family.
get_family_response_data = function(family, data_columns, data) {
  aux_columns = get_family_aux_columns(family, data_columns)
  out = lapply(aux_columns, function(column) {
    if (is.na(column) || column %notin% names(data)) NULL else data[[column]]
  })
  out[!vapply(out, is.null, logical(1))]
}


# Internal specification of a link function and its backend representations
mcplink = function(name, linkfun, linkinv, jags_link = name, jags_linkinv = "") {
  structure(
    list(
      name = name,
      linkfun = linkfun,
      linkinv = linkinv,
      jags_link = jags_link,
      jags_linkinv = jags_linkinv
    ),
    class = "mcplink"
  )
}

.mcplinks = list(
  identity = mcplink("identity", identity, identity, jags_link = "", jags_linkinv = ""),
  logit = mcplink("logit", stats::qlogis, stats::plogis, jags_link = "logit", jags_linkinv = "ilogit"),
  probit = mcplink("probit", stats::qnorm, stats::pnorm, jags_link = "probit", jags_linkinv = "phi"),
  log = mcplink("log", base::log, base::exp, jags_link = "log", jags_linkinv = "exp")
)

get_mcplink = function(link) {
  .mcplinks[[rlang::arg_match0(link, names(.mcplinks))]]
}


# Return the JAGS function name for a link or inverse link.
get_link_str = function(link, inverse = FALSE) {
  mcplink = get_mcplink(link)
  if (inverse) mcplink$jags_linkinv else mcplink$jags_link
}


# Return the executable R function for a link or inverse link.
get_link_function = function(link, inverse = FALSE) {
  mcplink = get_mcplink(link)
  if (inverse) mcplink$linkinv else mcplink$linkfun
}


#' @export
#' @noRd
format.mcpfamily = function(x, ...) {
  rlang::check_dots_empty()
  dpar_links = paste0(x$dpar_specs$dpar, " = ", x$dpar_specs$link, collapse = "; ")
  paste0("Family: ", x$family, "\nLinks: ", dpar_links)
}


#' @export
#' @noRd
print.mcpfamily = function(x, ...) {
  rlang::check_dots_empty()
  cat(format(x), "\n", sep = "")
  invisible(x)
}


is_mcpfamily = function(x) {
  inherits(x, "mcpfamily")
}


# Check that an mcpfamily implements the complete mcp family contract. Returns `x` invisibly.
assert_mcpfamily = function(x) {
  checkmate::assert_class(x, "mcpfamily")
  checkmate::assert_true(is_family(x), .var.name = "x")
  checkmate::assert_data_frame(x$default_prior)
  assert_data_cols(x$default_prior, c("dpar", "par_type", "prior"))
  assert_dpar_specs(x$dpar_specs)
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
  for (auxiliary in x$response$auxiliary) {
    checkmate::assert_flag(auxiliary$required)
    checkmate::assert_character(auxiliary$operations)
  }
  
  # R-side functions take `dpars` as rows x draws matrices on the response scale and `data` (response
  # auxiliaries, e.g., `trials`) with one value per row. Like stats::dbinom() etc., they must recycle `data` down the columns.
  checkmate::assert_function(x$r$epred)
  checkmate::assert_function(x$r$variance, null.ok = TRUE)
  checkmate::assert_function(x$r$log_lik)
  checkmate::assert_function(x$r$rng)
  checkmate::assert_function(x$r$cdf, null.ok = TRUE)
  checkmate::assert_function(x$r$quantile, null.ok = is.null(x$r$cdf))
  checkmate::assert_function(x$backends$jags$likelihood)
  if (!is.null(x$garma)) {
    checkmate::assert_function(x$garma$observed_r)
    checkmate::assert_function(x$garma$observed_jags)
    checkmate::assert_function(x$garma$validate_threshold, null.ok = TRUE)
  }
  checkmate::assert_function(x$linkfun)
  checkmate::assert_function(x$linkinv)

  invisible(x)
}
