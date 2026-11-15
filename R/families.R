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
  mcpfamily_bernoulli(link)
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
  mcpfamily_negbinomial(link, link_shape)
}


# Convert a stats family object (e.g., `binomial(link = "identity")`) into an `mcpfamily`.
# An `mcpfamily`, e.g., from bernoulli() or mcpfamily(), is returned as is.
as_mcpfamily = function(x) {
  checkmate::assert_true(is_family(x), .var.name = "x")
  if (is_mcpfamily(x))
    return(assert_mcpfamily(x))

  switch(
    x$family,
    gaussian = mcpfamily_gaussian(x$link),
    binomial = mcpfamily_binomial(x$link),
    poisson = mcpfamily_poisson(x$link),
    stop(
      "family = ", x$family, "() is not supported by mcp. Supported families are gaussian(), binomial(), ",
      "bernoulli(), poisson(), and negbinomial(). Use mcpfamily() to define other families.", call. = FALSE
    )
  )
}


#' Define a Family for mcp
#'
#' `r lifecycle::badge("experimental")`
#'
#' Define a response distribution for [mcp()], e.g., one that `mcp` does not ship.
#' The built-in families, such as `gaussian()` and [negbinomial()], are defined using `mcpfamily()` too.
#' See [Custom families](https://lindeloev.github.io/mcp/articles/custom_families.html) for a worked example.
#'
#' Only the elements needed to fit the model are required: `family`, `links`, and `jags$likelihood`, with priors supplied as defaults or in [mcp()].
#' Methods that are needed later (e.g., `r$rng` for `predict()`) raise an error on use.
#'
#' @param family String. Name the family for printing, e.g., `"student"`.
#' @param links Named character vector with the link function of each distributional parameter (dpar),
#'   e.g., `c(mu = "identity", sigma = "log")`. Supported links: `"identity"`, `"log"`, `"logit"`, and `"probit"`.
#'   The dpar `mu` is required.
#' @param jags Named list with the JAGS likelihood:
#'   * `likelihood`: A function taking a `context` list and returning JAGS code for observation `i_`
#'     (a string or a character vector of lines). See examples. `context` contains:
#'     * `y`: The response node, e.g., `"y[i_]"`.
#'     * `dpar(name)`: A dpar on the response scale, e.g., `context$dpar("sigma")` returns `"sigma_[i_]"`.
#'     * `aux(name, default)`: A response auxiliary column, e.g., `context$aux("weights", "1")`.
#'       Returns `default` if the model does not include it.
#'     * `local(name)`: A node for intermediate values, e.g., `context$local("prob")` returns `"prob_[i_]"`.
#'     * `zeros_trick(neg_log_density)`: For likelihoods that JAGS does not have. Takes the negative log-density of `y`
#'       as a string. Add a constant to keep it non-negative. Applies weights but does not support missing responses. See the [Gumbel example](https://lindeloev.github.io/mcp/articles/custom_families.html#likelihoods-that-jags-does-not-have).
#'   * `glm_module`: When to load the JAGS glm module: `"never"` (default), `"group_effects"`, or `"always"`.
#' @param r Named list of R functions. They take `dpars`, a list with a rows-by-draws matrix for each dpar
#'   on its response scale, and `data`, a list with the response auxiliary columns, e.g., `data$trials`.
#'   Like [stats::dnorm()], they must recycle `data` over the draws (columns). See examples and [Custom families article](https://lindeloev.github.io/mcp/articles/custom_families.html) for applied examples.
#'   * `log_lik(y, dpars, data)`: The log-density of each response, including normalizing constants.
#'     When supplied, `mcp` checks `jags$likelihood` against it after sampling. Required by [loo()], [waic()], and [log_lik()].
#'     For weighted models, the observation weights are applied by `mcp`.
#'   * `epred(dpars, data, rate = FALSE)`: The expected response. Used by `fitted()`, `plot()`, etc.
#'   * `rng(n, dpars, data, rate = FALSE)`: Random responses. Used by `predict()`, [pp_check()], `fit$simulate()`, etc.
#'   * `variance(dpars, data, rate = FALSE)`: The variance of the response. Used for exact `predict()` SDs.
#'     If missing, the SD of the `rng` draws is used.
#'   * `cdf(q, dpars, data)` and `quantile(p, dpars, data)`: Used for exact predictive intervals.
#'     If missing, the quantiles of the `rng` draws are used, which is less precise.
#' @param default_prior `NULL` or a data frame with columns `dpar`, `par_type` (`"Intercept"`, `"dummy"`, or
#'   `"slope"`), and `prior`, e.g., `"dnorm(0, 2.5 / predictor_scale())"`. Each parameter gets the row matching
#'   its `dpar` and `par_type` in [mcp_pars()]. Optional columns:
#'   * `group_sd_prior`: Prior for the SD of group-level effects.
#'   * `description`: Shown in `prior_summary(fit, verbose = TRUE)`.
#'   * `condition`: one of
#'     - `"always"`: (default)
#'     - `"included"`: for dpar terms explicitly included in formula, e.g., `sigma_Intercept_1` in `sigma(1)`.
#'     - `"not_included"` for dpar terms not included in any formula, e.g., `sigma_1` implicit when `sigma()` term is absent.
#'
#'   Priors use the syntax of `mcp(..., prior)` (see [mcp-priors]). `.y` and `.x` are the response and the change-point variable,
#'   e.g., `"dnorm(mean(.y), 2.5 * sd(.y))"`. `predictor_scale()` is the SD of each coefficient's predictor, for autoscaling.
#'   [prior_summary()] shows the result: try `prior_summary(demo_fit, verbose = TRUE)`, where `rule` is the default `prior`
#'   with these filled in, and `prior` is the resolved prior.
#' @param response Named list describing the response:
#'   * `auxiliary`: Named list of response auxiliaries, used like `y | trials(N)`, e.g.,
#'     `list(trials = list(required = TRUE, operations = c("epred", "log_lik", "rng")))`.
#'     `operations` lists where the column is needed (`"epred"`, `"log_lik"`, `"rng"`, and `"garma"`).
#'     `weights()` is always available. Use it in `jags$likelihood` to support weights.
#'   * `validate(y, data, data_columns)`: Function that raises an error for invalid responses.
#'   * `is_discrete`: Logical. Integer-valued responses? Defaults to `FALSE`.
#'   * `has_rate`, `observed(y, data, rate)`, `probability(rate)`: For responses per trial, as in `binomial()`.
#' @param lower Named numeric vector. The exclusive lower bound of dpars, e.g., `c(sigma = 0)`.
#'   Priors on constant dpars (e.g., `sigma_1`) are truncated to it.
#' @param floor Named numeric vector. Numerical floor on the dpar for each observation, e.g., `c(sigma = 1e-5)`.
#' @param garma `NULL` (AR/MA unsupported) or a named list with `observed_r(y, data, threshold)` and
#'   `observed_jags(context)`: The response, bounded by `threshold`, before the link function is applied, in R and JAGS.
#'   Optionally `validate_threshold(threshold, data)` which stop()s if invalid.
#' @param glm_family `NULL` or a stats family, e.g., `poisson()`, that approximates the model of `mu`.
#'   [stats::glm.fit()] uses it to find initial values for change points. If `NULL`, they are equally spaced.
#' @return An `mcpfamily` object.
#' @aliases mcpfamily
#' @export
#' @examples
#' # A Student-t family with mu, sigma, and nu
#' student = mcpfamily(
#'   family = "student",
#'   links = c(mu = "identity", sigma = "log", nu = "log"),
#'   jags = list(likelihood = function(context) {
#'     paste0(context$y, " ~ dt(", context$dpar("mu"), ", 1 / ", context$dpar("sigma"), "^2, ", context$dpar("nu"), ")")
#'   }),
#'   r = list(
#'     log_lik = function(y, dpars, data) stats::dt((y - dpars$mu) / dpars$sigma, dpars$nu, log = TRUE) - log(dpars$sigma),
#'     rng = function(n, dpars, data, rate = FALSE) dpars$mu + dpars$sigma * stats::rt(n, dpars$nu)
#'   ),
#'   lower = c(sigma = 0, nu = 1)
#' )
#'
#' # Priors are required for parameters without a default prior
#' data = data.frame(x = 1:20, y = c(rnorm(10), rnorm(10, 3)))
#' prior = list(Intercept_1 = "dnorm(0, 10)", Intercept_2 = "dnorm(0, 10)", sigma_1 = "dnorm(0, 5) T(0, )", nu_1 = "dgamma(2, 0.1) T(1, )")
#' fit = mcp(list(y ~ 1, ~ 1), data, family = student, prior = prior, par_x = "x", sample = FALSE)
#' fit$family
mcpfamily = function(family, links, jags, r = list(), default_prior = NULL, response = list(), lower = NULL, floor = NULL, garma = NULL, glm_family = NULL) {
  checkmate::assert_string(family)
  checkmate::assert_character(links, names = "unique", any.missing = FALSE)
  checkmate::assert_list(jags, names = "unique")
  checkmate::assert_list(r, names = "unique")
  checkmate::assert_numeric(lower, names = "unique", null.ok = TRUE)
  checkmate::assert_numeric(floor, names = "unique", null.ok = TRUE)
  if ("mu" %notin% names(links))
    stop("`links` must include the link function for 'mu', e.g., `links = c(mu = \"identity\")`.", call. = FALSE)

  # mu first. Other dpars are implicit: always part of the model, and constant (identity link) unless included in a formula.
  dpars = c("mu", setdiff(names(links), "mu"))
  dpar_specs = dplyr::bind_rows(lapply(dpars, function(dpar) {
    is_mu = dpar == "mu"
    new_dpar_spec(
      dpar,
      link_included = links[[dpar]],
      implicit = !is_mu,
      require_initial_term = !is_mu,
      lower = if (dpar %in% names(lower)) lower[[dpar]] else NA_real_,
      floor = if (dpar %in% names(floor)) floor[[dpar]] else NA_real_
    )
  }))

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

  if (is.null(jags$glm_module))
    jags$glm_module = "never"

  out = list(
    family = family,
    link = links[["mu"]],
    dpar_specs = dpar_specs,
    default_prior = standardize_family_default_priors(default_prior),
    response = response,
    r = r,
    backends = list(jags = jags),
    garma = garma,
    glm_family = glm_family,
    linkfun = get_link_function(links[["mu"]]),
    linkinv = get_link_function(links[["mu"]], inverse = TRUE)
  )
  class(out) = c("mcpfamily", "family")
  assert_mcpfamily(out)
  out
}


validate_count_mean = function(mu) {
  if (any(is.infinite(mu)))
    stop("Modelled count mean overflowed to Inf. Consider narrower priors.")
  invisible(TRUE)
}


standardize_family_default_priors = function(defaults) {
  if (is.null(defaults))
    defaults = data.frame(dpar = character(), par_type = character(), prior = character())
  checkmate::assert_data_frame(defaults)
  assert_data_cols(defaults, c("dpar", "par_type", "prior"))
  defaults = tibble::as_tibble(defaults)
  if (!"group_sd_prior" %in% names(defaults))
    defaults$group_sd_prior = NA_character_
  if (!"description" %in% names(defaults))
    defaults$description = NA_character_
  if (!"condition" %in% names(defaults))
    defaults$condition = "always"
  checkmate::assert_subset(defaults$condition, c("always", "included", "not_included"), .var.name = "default_prior$condition")
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


mcpfamily_gaussian = function(link) {
  if (link %notin% c("identity", "log"))
    stop("mcp has no default priors for gaussian(link = \"", link, "\") so it's likely not supported.")

  response_scale = "max(2.5, round(mad(.y), 1))"

  if (link == "identity") {
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
    "sigma", "Intercept", paste0("dt(0, ", response_scale, ", 3) T(0, )"), paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Positive residual SD calibrated on the response scale", "not_included",
    "sigma", "Intercept", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "Weakly regularizing log-SD intercept", "included",
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
        # get_jags_family_context()). The constant |w - 1| * 50 keeps phi >= 0 for
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

  garma = if (link == "identity") {
    list(
      observed_r = function(y, data, threshold) y,
      observed_jags = function(context) context$y
    )
  }

  mcpfamily(
    family = "gaussian",
    links = c(mu = link, sigma = "log"),
    jags = jags,
    r = r,
    default_prior = default_prior,
    response = list(),
    lower = c(sigma = 0),
    floor = c(sigma = 1e-5),
    garma = garma,
    glm_family = stats::gaussian(link)
  )
}


mcpfamily_binomial = function(link) {
  if (link == "identity") {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dunif(0, 1)", "dt(0, 1, 3) T(0, )", "Uniform success-probability intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "dt(0, 1, 3) T(0, )", "Success-probability difference between levels", "always",
      "mu", "slope", "dnorm(0, 1 / predictor_scale())", "dt(0, 1 / predictor_scale(), 3) T(0, )", "Autoscaled success-probability coefficient", "always"
    )
  } else if (link %in% c("logit", "probit")) {
    default_prior = probability_mu_priors(link)
  } else {
    stop("mcp has no default priors for binomial(link = \"", link, "\") so it's likely not supported.")
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
        c("# Binomial likelihood raised to the observation weight", context$zeros_trick(paste0("max(0, loggam(", y, " + 1) + loggam(", trials, " - ", y, " + 1) - loggam(", trials, " + 1) - ", y, " * log(", mu, ") - (", trials, " - ", y, ") * log(1 - ", mu, "))")))
      }
    }
  )
  garma = if (link == "logit") {
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

  mcpfamily(
    family = "binomial",
    links = c(mu = link),
    jags = jags,
    r = r,
    default_prior = default_prior,
    response = response,
    garma = garma,
    glm_family = stats::binomial(link)
  )
}


mcpfamily_bernoulli = function(link) {
  if (link == "identity") {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dunif(0, 1)", "dt(0, 1, 3) T(0, )", "Uniform P(y = TRUE) intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "dt(0, 1, 3) T(0, )", "P(y = TRUE) difference between levels", "always",
      "mu", "slope", "dnorm(0, 1 / predictor_scale())", "dt(0, 1 / predictor_scale(), 3) T(0, )", "Autoscaled P(y = TRUE) coefficient", "always"
    )
  } else if (link %in% c("logit", "probit")) {
    default_prior = probability_mu_priors(link)
  } else {
    stop("mcp has no default priors for bernoulli(link = \"", link, "\") so it's likely not supported.")
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
        c("# Bernoulli likelihood raised to the observation weight", context$zeros_trick(paste0("max(0, - ", y, " * log(", mu, ") - (1 - ", y, ") * log(1 - ", mu, "))")))
      }
    }
  )

  garma = if (link == "logit") {
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

  mcpfamily(
    family = "bernoulli",
    links = c(mu = link),
    jags = jags,
    r = r,
    default_prior = default_prior,
    response = response,
    garma = garma,
    glm_family = stats::binomial(link)
  )
}


mcpfamily_poisson = function(link) {
  response_scale = "max(2.5, round(mad(.y), 1))"
  if (link == "identity") {
    default_prior = tibble::tribble(
      ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
      "mu", "Intercept", "dnorm(mean(.y - .offset), 2.5 * sd(.y - .offset)) T(-min(.offset), )", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Count intercept on the response scale with a positive mean", "always",
      "mu", "dummy", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Autoscaled count contrast on the response scale", "always",
      "mu", "slope", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", paste0("dt(0, ", response_scale, " / predictor_scale(), 3) T(0, )"), "Autoscaled count coefficient on the response scale", "always"
    )
  } else if (link == "log") {
    default_prior = log_mu_priors("count")
  } else {
    stop("mcp has no default priors for poisson(link = \"", link, "\") so it's likely not supported.")
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
        c("# Poisson likelihood raised to the observation weight", context$zeros_trick(paste0("max(0, ", mu, " - ", y, " * log(", mu, ") + loggam(", y, " + 1))")))
      }
    }
  )
  garma = if (link == "log") {
    list(
      observed_r = function(y, data, threshold) pmax(y, threshold),
      observed_jags = function(context) paste0("max(", context$y, ", ", context$threshold, ")")
    )
  }

  mcpfamily(
    family = "poisson",
    links = c(mu = link),
    jags = jags,
    r = r,
    default_prior = default_prior,
    response = response,
    garma = garma,
    glm_family = stats::poisson(link)
  )
}


mcpfamily_negbinomial = function(link, link_shape) {
  if (link != "log" || link_shape != "log")
    stop("Negative-binomial models currently require log links for both mu and shape.")

  shape_prior = tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~description, ~condition,
    "shape", "Intercept", "dinvgamma(0.4, 0.3)", NA_character_, "Weakly regularizing positive overdispersion shape (brms default)", "not_included",
    "shape", "Intercept", "dnorm(0, 2.5)", "dnorm(0, 2.5) T(0, )", "Weakly regularizing log-shape intercept", "included",
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

      c(prob_code, "# Negative binomial likelihood raised to the observation weight", context$zeros_trick(paste0("max(0, loggam(", shape, ") + loggam(", y, " + 1) - loggam(", y, " + ", shape, ") - ", shape, " * log(", prob, ") - ", y, " * log(1 - ", prob, "))")))
    }
  )
  garma = list(
    observed_r = function(y, data, threshold) pmax(y, threshold),
    observed_jags = function(context) paste0("max(", context$y, ", ", context$threshold, ")")
  )

  mcpfamily(
    family = "negbinomial",
    links = c(mu = link, shape = link_shape),
    jags = jags,
    r = r,
    default_prior = default_prior,
    response = response,
    lower = c(shape = 0),
    garma = garma,
    glm_family = stats::poisson(link)  # Quasi-likelihood for the mean
  )
}


# Describe a distributional parameter. `link_included` applies when the dpar is included in a formula.
# Implicit dpars not included in any formula are constant with an identity link (see resolve_dpar_specs()).
# `lower` is the exclusive lower bound of its support (e.g., 0 for sigma) and bounds constant priors.
# `floor` is an optional numerical floor on per-observation values (see dev/DECISIONS.md for sigma).
new_dpar_spec = function(dpar, link_included, implicit = FALSE, require_initial_term = FALSE, lower = NA_real_, floor = NA_real_) {
  tibble::tibble(
    dpar = dpar,
    link = if (implicit) "identity" else link_included,
    link_included = link_included,
    included = FALSE,
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


# Lower bounds of constant dpars (not included in any formula) on their parameter scale, e.g., c(sigma = 0).
# Bounds that the link already enforces (e.g., log(0) = -Inf) are omitted.
get_constant_dpar_bounds = function(family) {
  dpar_specs = family$dpar_specs[!family$dpar_specs$included & !is.na(family$dpar_specs$lower), ]
  bounds = vapply(seq_len(nrow(dpar_specs)), function(i) get_link_function(dpar_specs$link[i])(dpar_specs$lower[i]), numeric(1))
  stats::setNames(bounds, dpar_specs$dpar)[is.finite(bounds)]
}


# Resolve each dpar link from whether the model includes a formula for it, e.g., sigma(1).
# - components: A list with `get_segment_components()` for each segment.
resolve_dpar_specs = function(family, components) {
  checkmate::assert_class(family, "mcpfamily")
  family$dpar_specs$included = family$dpar_specs$dpar %in% unlist(lapply(components, names))
  family$dpar_specs$link = ifelse(family$dpar_specs$included | !family$dpar_specs$implicit, family$dpar_specs$link_included, "identity")
  family
}


# Validate distributional-parameter metadata and supported links.
assert_dpar_specs = function(x) {
  checkmate::assert_data_frame(x)
  required = c(
    "dpar", "link", "link_included", "included",
    "implicit", "require_initial_term", "lower", "floor"
  )
  assert_data_cols(x, required)

  if (!is.character(x$dpar) || anyNA(x$dpar) || any(x$dpar == ""))
    stop("`family$dpar_specs$dpar` must contain non-empty parameter names.")
  link_cols = c("link", "link_included")
  if (any(!vapply(x[link_cols], is.character, logical(1))) || anyNA(x[link_cols]))
    stop("Links in `family$dpar_specs` must be names without missing values.")
  if (anyDuplicated(x$dpar))
    stop("Each distributional parameter must occur exactly once in `family$dpar_specs`.")
  if ("mu" %notin% x$dpar)
    stop("`family$dpar_specs` must contain the response-mean parameter 'mu'.")
  if (any(x$dpar %in% c("epred", "ar", "ma", "cp", "same")))
    stop("'epred', 'ar', 'ma', 'cp', and 'same' are reserved and cannot be family distributional parameters.")
  if (any(make.names(x$dpar) != x$dpar))
    stop("Distributional parameter names must be syntactic, e.g., 'nu' or 'zi'.")

  supported_links = names(.mcplinks)
  links = unique(unlist(x[link_cols], use.names = FALSE))
  if (any(links %notin% supported_links))
    stop("Unsupported dpar link(s): ", and_collapse(links[links %notin% supported_links]))
  if (!is.logical(x$included) || anyNA(x$included))
    stop("`family$dpar_specs$included` must be logical without missing values.")
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


# Get an R-side function of the family, e.g., `family_r(family, "rng")`. Families need not define those
# that are not needed to fit the model, so this errors if it is missing.
family_r = function(family, name) {
  fun = family$r[[name]]
  if (is.null(fun))
    stop("family = ", family$family, "() has no `r$", name, "` function, which is needed here. ",
         "Add it in `mcpfamily(..., r = list(", name, " = ...))` and refit.", call. = FALSE)
  fun
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
  # R-side functions are optional. Methods get them using family_r().
  for (name in c("log_lik", "epred", "variance", "rng", "cdf"))
    checkmate::assert_function(x$r[[name]], null.ok = TRUE, .var.name = paste0("r$", name))
  
  checkmate::assert_function(x$r$quantile, null.ok = is.null(x$r$cdf), .var.name = "r$quantile")
  checkmate::assert_function(x$backends$jags$likelihood, .var.name = "jags$likelihood")
  checkmate::assert_choice(x$backends$jags$glm_module, c("never", "group_effects", "always"), .var.name = "jags$glm_module")
  checkmate::assert_class(x$glm_family, "family", null.ok = TRUE, .var.name = "glm_family")
  
  if (!is.null(x$garma)) {
    checkmate::assert_function(x$garma$observed_r)
    checkmate::assert_function(x$garma$observed_jags)
    checkmate::assert_function(x$garma$validate_threshold, null.ok = TRUE)
  }
  checkmate::assert_function(x$linkfun)
  checkmate::assert_function(x$linkinv)

  invisible(x)
}
