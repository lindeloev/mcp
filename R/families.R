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
#' @param link_shape Link function for `shape` when it is modeled using `shape()`.
#' @details `shape` is required by the likelihood, so it includes an intercept in segment 1.
#'   Later segments that do not include `shape()` are joined. For example, `y ~ 1 + x + shape(1 + x)`
#'   models both the mean and shape. Regression coefficients for both dpars are
#'   on their link scales.
#'
#'   Without `shape()`, `shape_1` is the shape itself with an `inv_gamma(0.4, 0.3)` prior,
#'   as in `brms`, much like the implicit `sigma_1` in `gaussian()`.
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
#' `mcp`. Supported family and link combinations are:
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

  if (inherits(x, "mcpfamily")) {
    link_matches = is.null(x$links) || (
      identical(unname(x$link), unname(x$links[["mu"]])) &&
      (is.null(x$link_shape) || identical(unname(x$link_shape), x$dpar_specs$link_modeled[x$dpar_specs$dpar == "shape"]))
    )
    if (link_matches) {
      assert_mcpfamily(x)
      return(x)
    }
  }

  switch(
    x$family,
    gaussian = mcpfamily_gaussian(x),
    binomial = mcpfamily_binomial(x),
    bernoulli = mcpfamily_bernoulli(x),
    poisson = mcpfamily_poisson(x),
    negbinomial = mcpfamily_negbinomial(x),
    new_mcpfamily(x)
  )
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
    probability = function(rate) FALSE
  )

  # Fill in missing top-level response elements
  missing_response = setdiff(names(response_defaults), names(response))
  response[missing_response] = response_defaults[missing_response]

  # Merge missing auxiliary elements (e.g., weights) from response_defaults
  missing_aux = setdiff(names(response_defaults$auxiliary), names(response$auxiliary))
  response$auxiliary[missing_aux] = response_defaults$auxiliary[missing_aux]

  # Each family will define a validate(); this applies that.
  response$validate = make_weighted_validate(response$validate)

  # R-side log-likelihood. Weights are defined as applying directly here
  r$log_lik = make_weighted_log_lik(r$log_lik)

  family$dpar_specs = dpar_specs
  family$default_prior = standardize_family_default_priors(default_prior)
  family$response = response
  family$r = r
  family$backends = backends
  family$garma = garma
  family = add_dpar_specs(family)

  family$linkfun =get_link_function(family$link)
  family$linkinv = get_link_function(family$link, inverse = TRUE)

  class(family) = c("mcpfamily", "family")
  assert_mcpfamily(family)
  family
}


# Wrap family response validation to enforce positive, finite observation weights
make_weighted_validate = function(family_validate) {
  force(family_validate)
  function(y, data, response_columns) {
    # Check weights if present
    if (!is.null(data$weights)) {
      if (!is.numeric(data$weights) || anyNA(data$weights) || any(data$weights <= 0))
        stop("All weights must be numeric and greater than zero.")
    }

    # Apply the family-specific validation
    family_validate(y, data, response_columns)
  }
}


# Wrap family log-likelihood to apply observation weights and isolate closure state
make_weighted_log_lik = function(family_log_lik) {
  force(family_log_lik)
  function(y, dpars, data) {
    weights = if (is.null(data$weights)) 1 else data$weights
    weights * family_log_lik(y, dpars, data)
  }
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
  if (any(mu > 2146275819, na.rm = TRUE))
    stop("Modelled extremely large count mean (> 2146275819).")
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


mcpfamily_gaussian = function(family) {
  if (family$link %notin% c("identity", "log"))
    stop("mcp has no default priors for gaussian(link = \"", family$link, "\") so it's likely not supported.")

  response_scale = "max(2.5, round(mad(.y), 1))"

  if (family$link == "identity") {
    mu_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dnorm(mean(.y), 2.5 * sd(.y))", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Mean intercept (rstanarm default)", "always",
      "mu", "dummy", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Autoscaled mean contrast (rstanarm default)", "always",
      "mu", "slope", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, " / predictor_scale(), 3) T(0, )"), "Autoscaled mean coefficient (rstanarm default)", "always"
    )
  } else {
    log_y = "log(ifelse(.y <= 0, 0.1, .y))"
    mu_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", paste0("dnorm(round(median(", log_y, "), 1), max(2.5, round(mad(", log_y, "), 1)))"), "dnorm(0, 2.5) T(0, )", "Robustly centered log-mean intercept with a minimum scale of 2.5", "always",
      "mu", "dummy", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5) T(0, )", "Autoscaled log-mean contrast (rstanarm default)", "always",
      "mu", "slope", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5 / predictor_scale()) T(0, )", "Autoscaled log-mean coefficient (rstanarm default)", "always"
    )
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
        # The weighted normal contributes sigma^-1 from its normalizing term.
        # Raising the likelihood to weight w requires adding -(w - 1) * log(sigma)
        # to the log-likelihood. Computing dexp(sigma^(1 - w)) directly underflows/overflows
        # for large weights. Instead, we use the Poisson zeros trick on the log scale:
        # P(Z = 0) from dpois(phi) contributes -phi to the log-likelihood.
        # Setting phi = (w - 1) * log(sigma) + |w - 1| * 50 ensures phi >= 0 for any
        # reasonable sigma (e^-50 to e^50). The constant offset cancels out in MCMC.
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
      new_dpar_spec("sigma", "identity", implicit = TRUE, require_initial_predictor = TRUE, link_modeled = "log", lower = 0)
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
    # Intercept prior is approximately (logit) or exactly (probit) uniform on the probability scale
    scale = if (family$link == "logit") "1.5" else "1"
    uniform = if (family$link == "logit") "approximately uniform" else "uniform"
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", paste0("dnorm(0, ", scale, ")"), paste0("dt(0, ", scale, ", 3) T(0, )"), paste0("Link-scale intercept, ", uniform, " on the probability scale"), "always",
      "mu", "dummy", paste0("dnorm(0, ", scale, " / predictor_scale())"), paste0("dt(0, ", scale, ", 3) T(0, )"), "Autoscaled link-scale contrast (rstanarm autoscaling)", "always",
      "mu", "slope", paste0("dnorm(0, ", scale, " / predictor_scale())"), paste0("dt(0, ", scale, " / predictor_scale(), 3) T(0, )"), "Autoscaled link-scale coefficient (rstanarm autoscaling)", "always"
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
    # Intercept prior is approximately (logit) or exactly (probit) uniform on the probability scale
    scale = if (family$link == "logit") "1.5" else "1"
    uniform = if (family$link == "logit") "approximately uniform" else "uniform"
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", paste0("dnorm(0, ", scale, ")"), paste0("dt(0, ", scale, ", 3) T(0, )"), paste0("Link-scale intercept, ", uniform, " on the probability scale"), "always",
      "mu", "dummy", paste0("dnorm(0, ", scale, " / predictor_scale())"), paste0("dt(0, ", scale, ", 3) T(0, )"), "Autoscaled link-scale contrast (rstanarm autoscaling)", "always",
      "mu", "slope", paste0("dnorm(0, ", scale, " / predictor_scale())"), paste0("dt(0, ", scale, " / predictor_scale(), 3) T(0, )"), "Autoscaled link-scale coefficient (rstanarm autoscaling)", "always"
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
      "mu", "Intercept", "dnorm(mean(.y), 2.5 * sd(.y)) T(0, )", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Positive count intercept on the response scale", "always",
      "mu", "dummy", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Autoscaled count contrast on the response scale", "always",
      "mu", "slope", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, " / predictor_scale(), 3) T(0, )"), "Autoscaled count coefficient on the response scale", "always"
    )
  } else if (family$link == "log") {
    count_y = "log(ifelse(.y <= 0, 0.1, .y))"
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", paste0("dnorm(round(median(", count_y, "), 1), max(2.5, round(mad(", count_y, "), 1)))"), "dnorm(0, 2.5) T(0, )", "Robustly centered log-count intercept with a minimum scale of 2.5", "always",
      "mu", "dummy", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5) T(0, )", "Autoscaled log-count contrast (rstanarm default)", "always",
      "mu", "slope", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5 / predictor_scale()) T(0, )", "Autoscaled log-count coefficient (rstanarm default)", "always"
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

  count_y = "log(ifelse(.y <= 0, 0.1, .y))"
  default_prior = tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
    "mu", "Intercept", paste0("dnorm(round(median(", count_y, "), 1), max(2.5, round(mad(", count_y, "), 1)))"), "dnorm(0, 2.5) T(0, )", "Robustly centered log-count intercept with a minimum scale of 2.5", "always",
    "mu", "dummy", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5) T(0, )", "Autoscaled log-count contrast (rstanarm default)", "always",
    "mu", "slope", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5 / predictor_scale()) T(0, )", "Autoscaled log-count coefficient (rstanarm default)", "always",
    "shape", "Intercept", "dinvgamma(0.4, 0.3)", NA_character_, "Weakly regularizing positive overdispersion shape (brms default)", "constant",
    "shape", "Intercept", "dnorm(0, 2.5)", "dnorm(0, 2.5) T(0, )", "Weakly regularizing modeled log-shape intercept", "modeled",
    "shape", "dummy", "dnorm(0, 2.5)", "dnorm(0, 2.5) T(0, )", "Regularizing shape contrast on the log scale", "always",
    "shape", "slope", "dnorm(0, 2.5 / predictor_scale())", "dnorm(0, 2.5 / predictor_scale()) T(0, )", "Regularizing log-shape coefficient scaled by the predictor SD", "always"
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
      new_dpar_spec("shape", "identity", implicit = TRUE, require_initial_predictor = TRUE, link_modeled = family$link_shape, lower = 0)
    ),
    default_prior = default_prior,
    response = response,
    r = r,
    backends = list(jags = jags),
    garma = garma
  )
}


# Describe a distributional parameter. `lower` is the exclusive lower bound of its
# support (e.g., 0 for sigma). It bounds constant priors and floors observation-level values.
new_dpar_spec = function(dpar, link, implicit = FALSE, require_initial_predictor = FALSE, lower = NA_real_,
                          link_modeled = link) {
  tibble::tibble(
    dpar = dpar,
    link = link,
    link_constant = link,
    link_modeled = link_modeled,
    modeled = FALSE,
    implicit = implicit,
    require_initial_predictor = require_initial_predictor,
    lower = lower
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
  specs = family$dpar_specs[!family$dpar_specs$modeled & !is.na(family$dpar_specs$lower), ]
  bounds = vapply(seq_len(nrow(specs)), function(i) get_link_function(specs$link[i])(specs$lower[i]), numeric(1))
  stats::setNames(bounds, specs$dpar)[is.finite(bounds)]
}


# Standardize and expose a family's distributional-parameter metadata.
add_dpar_specs = function(family) {
  assert_dpar_specs(family$dpar_specs)
  family$dpars = family$dpar_specs$dpar
  family$links = stats::setNames(family$dpar_specs$link, family$dpar_specs$dpar)
  family
}


# Resolve each dpar link from whether the model includes a formula for it, e.g., sigma(1).
# - components: A list with `get_segment_components()` for each segment.
resolve_dpar_specs = function(family, components) {
  checkmate::assert_class(family, "mcpfamily")
  family = add_dpar_specs(family)

  family$dpar_specs$modeled = family$dpar_specs$dpar %in% unlist(lapply(components, names))
  family$dpar_specs$link = ifelse(
    family$dpar_specs$modeled,
    family$dpar_specs$link_modeled,
    family$dpar_specs$link_constant
  )

  add_dpar_specs(family)
}


# Validate distributional-parameter metadata and supported links.
assert_dpar_specs = function(x) {
  checkmate::assert_data_frame(x)
  required = c(
    "dpar", "link", "link_constant", "link_modeled", "modeled",
    "implicit", "require_initial_predictor", "lower"
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
  if (!is.logical(x$require_initial_predictor) || anyNA(x$require_initial_predictor))
    stop("`family$dpar_specs$require_initial_predictor` must be logical without missing values.")
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


# Extract the response-auxiliary columns defined by a fitted family.
get_family_response_data = function(family, segments, data) {
  aux_columns = get_family_aux_columns(family, segments)
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
  obj = get_mcplink(link)
  if (inverse) obj$jags_linkinv else obj$jags_link
}


# Return the executable R function for a link or inverse link.
get_link_function = function(link, inverse = FALSE) {
  obj = get_mcplink(link)
  if (inverse) obj$linkinv else obj$linkfun
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


#' @aliases is.mcpfamily
#' @describeIn mcpfamily Checks whether x is an `mcpfamily`.
#' @export
is.mcpfamily = function(x) {
  inherits(x, "mcpfamily")
}


# Check that an mcpfamily implements the complete mcp family contract. Returns `x` invisibly.
assert_mcpfamily = function(x) {
  checkmate::assert_class(x, "mcpfamily")
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
  for (auxiliary in x$response$auxiliary) {
    checkmate::assert_flag(auxiliary$required)
    checkmate::assert_character(auxiliary$operations)
  }
  checkmate::assert_function(x$r$epred)
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

# This renaming enables fit$simulate() with JAGS-like code
logit = stats::qlogis

ilogit = stats::plogis

probit = stats::qnorm

phi = stats::pnorm
