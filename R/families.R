# Family specifications ----------------------------------------------------

#' Bernoulli Family for mcp
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
bernoulli = function(link = "logit") {
  assert_value(link, allowed = c("identity", "logit", "probit"))

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
negbinomial = function(link = "log", link_shape = "log") {
  assert_value(link, allowed = "log")
  assert_value(link_shape, allowed = "log")

  family = list(family = "negbinomial", link = link, link_shape = link_shape)
  class(family) = "family"
  mcpfamily(family)
}


#' Create or Test Objects of Class "mcpfamily"
#'
#' @aliases mcpfamily
#' @param x A family object, e.g., `binomial(link = "identity")`.
#' @seealso \code{\link{family}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
mcpfamily = function(x) {
  assert_types(x, "family")

  family = switch(
    x$family,
    gaussian = mcpfamily_gaussian(x),
    binomial = mcpfamily_binomial(x),
    bernoulli = mcpfamily_bernoulli(x),
    poisson = mcpfamily_poisson(x),
    negbinomial = mcpfamily_negbinomial(x),
    new_mcpfamily(x)
  )

  family
}


new_mcpfamily = function(family, dpar_specs = family$dpar_specs,
                         default_prior = family$default_prior,
                         response = family$response, r = family$r,
                         backends = family$backends, garma = family$garma) {
  if (is.null(dpar_specs) || is.null(default_prior) || is.null(response) ||
      is.null(r) || is.null(backends$jags)) {
    stop(
      "family = ", family$family, "() does not have a complete internal mcp ",
      "family specification."
    )
  }

  response_defaults = list(
    auxiliary = list(),
    validate = function(y, data, columns) invisible(TRUE),
    observed = function(y, data, rate) y,
    probability = function(rate) FALSE,
    point_size = NULL
  )
  missing_response = setdiff(names(response_defaults), names(response))
  response[missing_response] = response_defaults[missing_response]

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
  assert_types(defaults, "tibble", "data.frame")
  assert_data_cols(defaults, c("dpar", "par_type", "prior"))
  defaults = tibble::as_tibble(defaults)
  if (!"description" %in% names(defaults))
    defaults$description = NA_character_
  if (!"condition" %in% names(defaults))
    defaults$condition = "always"
  defaults[, c("dpar", "par_type", "prior", "description", "condition")]
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
    mu_location = "round(log(max(median(.y), 0.1)), 1)"
    mu_scale = "2.5"
  }

  default_prior = tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "mu", "Intercept", paste0("dt(", mu_location, ", ", mu_scale, ", 3)"), "Robustly centered mean intercept with a brms-inspired minimum scale", "always",
    "mu", "dummy", paste0("dt(0, ", mu_scale, ", 3)"), "Regularizing mean contrast on the link scale", "always",
    "mu", "slope", paste0("dt(0, ", mu_scale, " / segment_width(.x), 3)"), "Regularizing mean slope scaled to the expected segment width", "always",
    "sigma", "Intercept", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Positive residual SD calibrated on the response scale", "always",
    "sigma", "dummy", paste0("dt(0, ", response_scale, ", 3)"), "Residual-SD contrast calibrated on the response scale", "always",
    "sigma", "slope", paste0("dt(0, ", response_scale, " / segment_width(.x), 3)"), "Residual-SD slope scaled to the expected segment width", "always"
  )

  response = list(
    auxiliary = list(weights = list(
      required = FALSE,
      operations = c("log_lik", "rng", "garma")
    )),
    validate = function(y, data, columns) {
      if (!is.null(data$weights)) {
        if (!is.numeric(data$weights) || anyNA(data$weights) || any(data$weights <= 0))
          stop("All weights must be numeric and greater than zero.")
      }
      invisible(TRUE)
    },
    point_size = "weights"
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    log_lik = function(y, dpars, data) {
      weights = if (is.null(data$weights)) 1 else data$weights
      stats::dnorm(y, dpars$mu, dpars$sigma / sqrt(weights), log = TRUE)
    },
    rng = function(n, dpars, data, rate = FALSE) {
      weights = if (is.null(data$weights)) 1 else data$weights
      stats::rnorm(n, dpars$mu, dpars$sigma / sqrt(weights))
    }
  )
  jags = list(likelihood = function(context) paste0(context$y, " ~ dnorm(", context$dpar("mu"), ", ", context$aux("weights", "1"), " / ", context$dpar("sigma"), "^2)  # SD as precision"))
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
      new_dpar_spec("sigma", "identity", implicit = TRUE, lower = 1e-9)
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
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", "dbeta(1, 1)", "Uniform probability intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "Probability difference between levels", "always",
      "mu", "slope", "dt(0, 1 / segment_width(.x), 3)", "Current slope prior scaled to the expected change-point segment width", "always"
    )
  } else if (family$link %in% c("logit", "probit")) {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", "dt(0, 2.5, 3)", "brms-inspired weak prior for the link-scale intercept", "always",
      "mu", "dummy", "dt(0, 2.5, 3)", "Regularizing categorical contrast on the link scale", "always",
      "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing link-scale slope scaled to the expected segment width", "always"
    )
  } else {
    stop("mcp has no default priors for binomial(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    auxiliary = list(trials = list(
      required = TRUE,
      operations = c("epred", "log_lik", "rng", "garma")
    )),
    validate = function(y, data, columns) {
      assert_integer(y, columns$y, lower = 0)
      assert_integer(data$trials, columns$trials, lower = 1)
      invalid = which(!is.na(y) & !is.na(data$trials) & y > data$trials)
      if (length(invalid) > 0)
        stop(
          "For family = binomial(), responses in '", columns$y,
          "' cannot exceed trials in '", columns$trials,
          "'. Found invalid data in row(s): ", paste(invalid, collapse = ", "), "."
        )
      invisible(TRUE)
    },
    observed = function(y, data, rate) if (rate) y / data$trials else y,
    probability = function(rate) isTRUE(rate)
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) if (rate) dpars$mu else data$trials * dpars$mu,
    log_lik = function(y, dpars, data) stats::dbinom(y, data$trials, dpars$mu, log = TRUE),
    rng = function(n, dpars, data, rate = FALSE) {
      y = stats::rbinom(n, data$trials, dpars$mu)
      if (rate) y / data$trials else y
    }
  )
  jags = list(likelihood = function(context) paste0(context$y, " ~ dbin(", context$dpar("mu"), ", ", context$aux("trials"), ")"))
  garma = if (family$link == "logit") {
    list(
      observed_r = function(y, data, boundary) pmin(pmax(y, boundary), data$trials - boundary) / data$trials,
      observed_jags = function(context) paste0("min(max(", context$y, ", ", context$boundary, "), ", context$aux("trials"), " - ", context$boundary, ") / ", context$aux("trials"))
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
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", "dbeta(1, 1)", "Uniform probability intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "Probability difference between levels", "always",
      "mu", "slope", "dt(0, 1 / segment_width(.x), 3)", "Current slope prior scaled to the expected change-point segment width", "always"
    )
  } else if (family$link %in% c("logit", "probit")) {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", "dt(0, 2.5, 3)", "brms-inspired weak prior for the link-scale intercept", "always",
      "mu", "dummy", "dt(0, 2.5, 3)", "Regularizing categorical contrast on the link scale", "always",
      "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing link-scale slope scaled to the expected segment width", "always"
    )
  } else {
    stop("mcp has no default priors for bernoulli(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    validate = function(y, data, columns) {
      invalid = which(!is.na(y) & y %notin% c(0, 1))
      if (length(invalid) > 0)
        stop("Only responses 0 and 1 are allowed for family = bernoulli() in column '", columns$y, "'")
      invisible(TRUE)
    },
    probability = function(rate) TRUE
  )
  r = list(
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    log_lik = function(y, dpars, data) stats::dbinom(y, 1, dpars$mu, log = TRUE),
    rng = function(n, dpars, data, rate = FALSE) stats::rbinom(n, 1, dpars$mu)
  )
  jags = list(likelihood = function(context) paste0(context$y, " ~ dbern(", context$dpar("mu"), ")"))

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
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", paste0("dt(", response_location, ", ", response_scale, ", 3) T(0, )"), "Positive count intercept calibrated on the response scale", "always",
      "mu", "dummy", paste0("dt(0, ", response_scale, ", 3)"), "Count contrast calibrated on the response scale", "always",
      "mu", "slope", paste0("dt(0, ", response_scale, " / segment_width(.x), 3)"), "Count slope scaled to the expected segment width", "always"
    )
  } else if (family$link == "log") {
    count_y = "log(pmax(.y, 0.1))"
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", paste0("dt(round(median(", count_y, "), 1), max(2.5, round(mad(", count_y, "), 1)), 3)"), "Robustly centered log-count intercept with a brms-inspired minimum scale", "always",
      "mu", "dummy", "dt(0, 2.5, 3)", "Regularizing categorical contrast on the log scale", "always",
      "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing log-count slope scaled to the expected segment width", "always"
    )
  } else {
    stop("mcp has no default priors for poisson(link = \"", family$link, "\") so it's likely not supported.")
  }

  response = list(
    validate = function(y, data, columns) {
      assert_integer(y, columns$y, lower = 0)
      invisible(TRUE)
    }
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
    }
  )
  jags = list(likelihood = function(context) paste0(context$y, " ~ dpois(", context$dpar("mu"), ")"))
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
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "mu", "Intercept", paste0("dt(round(median(", count_y, "), 1), max(2.5, round(mad(", count_y, "), 1)), 3)"), "Robustly centered log-count intercept with a brms-inspired minimum scale", "always",
    "mu", "dummy", "dt(0, 2.5, 3)", "Regularizing categorical count contrast on the log scale", "always",
    "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing log-count slope scaled to the expected segment width", "always",
    "shape", "Intercept", "dloginvgamma(0.4, 0.3)", "brms-compatible negative-binomial shape prior", "constant",
    "shape", "Intercept", "dt(0, 2.5, 3)", "brms-inspired prior for a modeled log-shape intercept", "modeled",
    "shape", "dummy", "dt(0, 2.5, 3)", "Regularizing shape contrast on the log scale", "always",
    "shape", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing log-shape slope scaled to the expected segment width", "always"
  )
  response = list(
    validate = function(y, data, columns) {
      assert_integer(y, columns$y, lower = 0)
      invisible(TRUE)
    }
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
    }
  )
  jags = list(likelihood = function(context) {
    mu = context$dpar("mu")
    shape = context$dpar("shape")
    prob = context$local("nb_prob")
    c(
      paste0(prob, " = ", shape, " / (", shape, " + ", mu, ")"),
      paste0(context$y, " ~ dnegbin(", prob, ", ", shape, ")")
    )
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


#' Describe a distributional parameter
#'
#' @keywords internal
#' @noRd
new_dpar_spec = function(dpar, link, implicit = FALSE, lower = NA_real_) {
  tibble::tibble(dpar = dpar, link = link, implicit = implicit, lower = lower)
}


known_dpar_wrappers = function() {
  c("sigma", "shape")
}


add_dpar_specs = function(family) {
  assert_dpar_specs(family$dpar_specs)
  family$dpars = family$dpar_specs$dpar
  family$links = stats::setNames(family$dpar_specs$link, family$dpar_specs$dpar)
  family
}


assert_dpar_specs = function(x) {
  assert_types(x, "data.frame", "tibble")
  required = c("dpar", "link", "implicit", "lower")
  assert_data_cols(x, required)

  if (!is.character(x$dpar) || anyNA(x$dpar) || any(x$dpar == ""))
    stop("`family$dpar_specs$dpar` must contain non-empty parameter names.")
  if (!is.character(x$link) || anyNA(x$link))
    stop("`family$dpar_specs$link` must contain link names without missing values.")
  if (anyDuplicated(x$dpar))
    stop("Each distributional parameter must occur exactly once in `family$dpar_specs`.")
  if ("mu" %notin% x$dpar)
    stop("`family$dpar_specs` must contain the response-mean parameter 'mu'.")
  if (any(x$dpar %in% c("epred", "ar", "ma")))
    stop("'epred', 'ar', and 'ma' are reserved and cannot be family distributional parameters.")

  supported_links = c("identity", "log", "logit", "probit")
  if (any(x$link %notin% supported_links))
    stop("Unsupported dpar link(s): ", and_collapse(unique(x$link[x$link %notin% supported_links])))
  if (!is.logical(x$implicit) || anyNA(x$implicit))
    stop("`family$dpar_specs$implicit` must be logical without missing values.")
  if (!is.numeric(x$lower))
    stop("`family$dpar_specs$lower` must be numeric.")

  TRUE
}


get_dpar_spec = function(family, dpar) {
  spec = family$dpar_specs[family$dpar_specs$dpar == dpar, , drop = FALSE]
  if (nrow(spec) != 1)
    stop_github("Expected exactly one dpar specification for '", dpar, "'.")
  spec
}


get_family_aux_columns = function(family, ST, operations = NULL) {
  auxiliary = family$response$auxiliary
  if (!is.null(operations)) {
    used = vapply(auxiliary, function(x) any(x$operations %in% operations), logical(1))
    auxiliary = auxiliary[used]
  }
  aux_names = names(auxiliary)
  if (length(aux_names) == 0)
    return(stats::setNames(character(), character()))

  stats::setNames(vapply(aux_names, function(name) {
    if (name %notin% names(ST))
      return(NA_character_)
    cols = unique(stats::na.omit(ST[[name]]))
    if (length(cols) > 1)
      stop("There should be exactly zero or one column used for ", name, "().")
    if (length(cols) == 0) NA_character_ else cols
  }, character(1)), aux_names)
}


get_family_response_data = function(family, ST, data) {
  columns = get_family_aux_columns(family, ST)
  out = lapply(columns, function(column) {
    if (is.na(column) || column %notin% names(data)) NULL else data[[column]]
  })
  out[!vapply(out, is.null, logical(1))]
}


get_link_str = function(link, inverse = FALSE) {
  if (!inverse)
    return(ifelse(link == "identity", "", link))

  switch(link, logit = "ilogit", probit = "phi", log = "exp", identity = "")
}


get_link_function = function(link, inverse = FALSE) {
  link_object = stats::make.link(link)
  if (inverse) link_object$linkinv else link_object$linkfun
}


#' @aliases is.mcpfamily
#' @describeIn mcpfamily Checks whether x is an `mcpfamily`.
#' @export
is.mcpfamily = function(x) {
  if (!inherits(x, "mcpfamily"))
    return(FALSE)

  required = c("default_prior", "dpar_specs", "dpars", "links", "response", "r", "backends")
  if (any(required %notin% names(x)))
    return(FALSE)

  assert_types(x, "family")
  assert_types(x$default_prior, "tibble", "data.frame")
  assert_data_cols(x$default_prior, c("dpar", "par_type", "prior"))
  assert_dpar_specs(x$dpar_specs)
  assert_types(x$dpars, "character")
  assert_types(x$links, "character")
  assert_types(x$response$auxiliary, "list")
  assert_types(x$response$validate, "function")
  assert_types(x$response$observed, "function")
  assert_types(x$response$probability, "function")
  assert_types(x$response$point_size, "null", "character", len = c(0, 1))
  for (auxiliary in x$response$auxiliary) {
    assert_logical(auxiliary$required, len = 1)
    assert_types(auxiliary$operations, "character")
  }
  assert_types(x$r$epred, "function")
  assert_types(x$r$log_lik, "function")
  assert_types(x$r$rng, "function")
  assert_types(x$backends$jags$likelihood, "function")
  if (!is.null(x$garma)) {
    assert_types(x$garma$observed_r, "function")
    assert_types(x$garma$observed_jags, "function")
    assert_types(x$garma$generate_message, "null", "character", len = c(0, 1))
  }
  assert_types(x$linkfun_str, "character", len = 1)
  assert_types(x$linkinv_str, "character", len = 1)
  assert_types(x$linkfun, "function")
  assert_types(x$linkinv, "function")

  TRUE
}


#' Logit function
#'
#' @aliases logit
#' @param mu A vector of probabilities (0.0 to 1.0)
#' @return A vector with same length as `mu`
#' @export
logit = stats::binomial(link = "logit")$linkfun

#' Inverse logit function
#'
#' @aliases ilogit
#' @param eta A vector of logits
#' @return A vector with same length as `eta`
#' @export
ilogit = stats::binomial(link = "logit")$linkinv

#' Probit function
#'
#' @aliases probit
#' @param mu A vector of probabilities (0.0 to 1.0)
#' @return A vector with same length as `mu`
#' @export
probit = stats::binomial(link = "probit")$linkfun

#' Inverse probit function
#'
#' @aliases phi
#' @param eta A vector of probits
#' @return A vector with same length as `mu`
#' @export
phi = stats::binomial(link = "probit")$linkinv
