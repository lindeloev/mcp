# A Student-t family like in vignettes/custom_families.Rmd. Default priors for mu, sigma, and nu intercepts with nu() not included/included.
student = mcpfamily(
  family = "student",
  links = c(mu = "identity", sigma = "log", nu = "log"),
  jags = list(likelihood = function(context) {
    paste0(context$y, " ~ dt(", context$dpar("mu"), ", 1 / ", context$dpar("sigma"), "^2, ", context$dpar("nu"), ")")
  }),
  r = list(
    log_lik = function(y, dpars, data) stats::dt((y - dpars$mu) / dpars$sigma, dpars$nu, log = TRUE) - log(dpars$sigma),
    epred = function(dpars, data, rate = FALSE) dpars$mu,
    rng = function(n, dpars, data, rate = FALSE) dpars$mu + dpars$sigma * stats::rt(n, dpars$nu)
  ),
  default_prior = tibble::tribble(
    ~dpar, ~par_type, ~prior, ~group_sd_prior, ~condition,
    "mu", "Intercept", "dnorm(mean(.y), 2.5 * sd(.y))", "dt(0, sd(.y), 3) T(0, )", "always",
    "mu", "dummy", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", "dt(0, sd(.y), 3) T(0, )", "always",
    "mu", "slope", "dnorm(0, 2.5 * sd(.y) / predictor_scale())", "dt(0, sd(.y) / predictor_scale(), 3) T(0, )", "always",
    "sigma", "Intercept", "dt(0, sd(.y), 3) T(0, )", NA, "not_included",
    "sigma", "Intercept", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "included",
    "sigma", "dummy", "dt(0, 2.5, 3)", "dt(0, 2.5, 3) T(0, )", "always",
    "sigma", "slope", "dt(0, 2.5 / predictor_scale(), 3)", "dt(0, 2.5 / predictor_scale(), 3) T(0, )", "always",
    "nu", "Intercept", "dgamma(2, 0.1) T(1, )", NA, "not_included",
    "nu", "Intercept", "dnorm(2.5, 1)", NA, "included"
  ),
  lower = c(sigma = 0, nu = 1)
)


bad_student = list(
  list(y | weights(weights_ok) ~ 1),  # The likelihood does not use weights
  list(y ~ 1 + ar(1)),  # No `garma`
  list(y ~ 1 + nu(1 + x)),  # No default prior for nu_x_1
  list(y ~ 1 + shape(1))
)

test_bad(
  bad_student,
  family = student
)


good_student = list(
  list(y ~ 1 + x),
  list(
    y ~ 1 + x,
    ~ 1 + x
  ),
  list(y ~ 1 + sigma(1 + x)),
  list(
    y ~ 1 + nu(1),
    ~ 0 + x
  ),
  list(y ~ 1 + (1 | id))
)

test_good(
  good_student,
  family = student
)


test_that("custom families work with all mcpfit methods", {
  test_runs(list(y ~ 1 + x, ~ 1 + x), family = student, test_s3 = TRUE)
})


# The same distribution using the zeros trick, which supports weights. The constant 100 keeps it non-negative.
student_zeros = mcpfamily(
  family = "student_zeros",
  links = c(mu = "identity", sigma = "log", nu = "log"),
  jags = list(likelihood = function(context) {
    nu = context$dpar("nu")
    sigma = context$dpar("sigma")
    context$zeros_trick(paste0(
      "loggam(", nu, " / 2) - loggam((", nu, " + 1) / 2) + 0.5 * log(", nu, " * 3.141593) + log(", sigma, ") + (", nu, " + 1) / 2 * ",
      "log(1 + ((", context$y, " - ", context$dpar("mu"), ") / ", sigma, ")^2 / ", nu, ") + 100"
    ))
  }),
  r = student$r,
  default_prior = student$default_prior,
  lower = c(sigma = 0, nu = 1)
)

good_student_zeros = list(
  list(y ~ 1 + x),
  list(
    y ~ 1 + x,
    ~ 1 + sigma(1)
  ),
  list(y | weights(weights_ok) ~ 1 + x)
)

test_good(
  good_student_zeros,
  family = student_zeros
)


test_that("zeros-trick likelihoods keep the response as data", {
  fit = mcp(list(y ~ 1 + x), data_gauss, family = student_zeros, par_x = "x", sample = FALSE)
  expect_match(get_fit_jags_code(fit), "likelihood_phi_[i_] = response_observed_[i_] * (loggam(", fixed = TRUE)
  prior_jags = as_jags(fit, prior = TRUE)
  expect_identical(prior_jags$data$y, data_gauss$y)
  expect_identical(prior_jags$data$response_observed_, rep(0L, nrow(data_gauss)))

  data_na = transform(data_gauss, y = replace(y, 2, NA))
  expect_error(
    suppressMessages(mcp(list(y ~ 1 + x), data_na, family = student_zeros, par_x = "x")),
    "family = student_zeros() uses the response as data, e.g., to apply weights, so it cannot impute missing responses",
    fixed = TRUE
  )
})


test_that("mcpfamily() validates its arguments", {
  likelihood = function(context) paste0(context$y, " ~ dnorm(", context$dpar("mu"), ", 1)")
  r = list(log_lik = function(y, dpars, data) stats::dnorm(y, dpars$mu, 1, log = TRUE))

  expect_s3_class(mcpfamily("unit", c(mu = "identity"), list(likelihood = likelihood), r), c("mcpfamily", "family"))
  expect_error(mcpfamily("unit", c(sigma = "log"), list(likelihood = likelihood), r), "must include the link function for 'mu'")
  expect_error(mcpfamily("unit", c(mu = "cloglog"), list(likelihood = likelihood), r), "must be one of")
  expect_error(mcpfamily("unit", c(mu = "identity", cp = "log"), list(likelihood = likelihood), r), "reserved")
  expect_error(mcpfamily("unit", c(mu = "identity", `1nu` = "log"), list(likelihood = likelihood), r), "must be syntactic")
  expect_s3_class(mcpfamily("unit", c(mu = "identity"), list(likelihood = likelihood)), "mcpfamily")
  expect_error(mcpfamily("unit", c(mu = "identity"), list(likelihood = likelihood), list(log_lik = 1)), "r\\$log_lik")
  expect_error(mcpfamily("unit", c(mu = "identity"), list(), r), "jags\\$likelihood")
  expect_error(mcpfamily("unit", c(mu = "identity"), list(likelihood = likelihood, glm_module = "sometimes"), r), "jags\\$glm_module")

  # Built-in families are mcpfamily() objects and returned as is
  expect_identical(as_mcpfamily(student), student)
  expect_identical(as_mcpfamily(negbinomial())$dpar_specs, negbinomial()$dpar_specs)
})


test_that("families without default priors require specified priors", {
  minimal = mcpfamily(
    family = "minimal",
    links = c(mu = "identity", sigma = "log", nu = "log"),
    jags = student$backends$jags,
    r = list(log_lik = student$r$log_lik),
    lower = c(sigma = 0, nu = 1)
  )
  model = list(y ~ 1 + x)
  expect_error(
    mcp(model, data_gauss, family = minimal, par_x = "x", sample = FALSE),
    "family = minimal() has no default prior for `Intercept_1` and `x_1` and `nu_1` and `sigma_1`",
    fixed = TRUE
  )

  prior = list(Intercept_1 = "dnorm(0, 10)", x_1 = 0.5, sigma_1 = "dnorm(0, 5)", nu_1 = "dgamma(2, 0.1)")
  fit = mcp(model, data_gauss, family = minimal, prior = prior, par_x = "x", sample = FALSE)
  table = prior_summary(fit, verbose = TRUE)
  expect_equal(table$prior, c("dnorm(0, 10)", "0.5", "dnorm(0, 5) T(0, )", "dgamma(2, 0.1) T(1, )"))  # Truncated to `lower`
  expect_true(all(table$source == "user"))

  # Methods error when they need a missing R-side function
  expect_error(
    fit$simulate(fit, data_gauss, Intercept_1 = 0, x_1 = 1, sigma_1 = 1, nu_1 = 3, .type = "predict"),
    "family = minimal() has no `r$rng` function",
    fixed = TRUE
  )
})


test_that("custom families can fit without R-side functions", {
  minimal = mcpfamily(
    family = "minimal",
    links = c(mu = "identity", sigma = "log", nu = "log"),
    jags = student$backends$jags,
    default_prior = student$default_prior,
    lower = c(sigma = 0, nu = 1)
  )
  expect_no_warning({
    fit = mcp(list(y ~ 1 + x), data_gauss, family = minimal, par_x = "x", chains = 1, warmup = 300, iter = 50, diagnostics = FALSE, quiet = TRUE)
  })
  expect_s3_class(coda::as.mcmc(fit), "mcmc.list")
  expect_silent(check_jags_likelihood(fit, NULL))
  for (method in list(log_lik, loo, waic))
    expect_error(method(fit), "family = minimal() has no `r$log_lik` function", fixed = TRUE)
})


test_that("mcp() warns when the JAGS likelihood and R-side log_lik() disagree", {
  wrong = student
  wrong$r$log_lik = function(y, dpars, data) stats::dnorm(y, dpars$mu, dpars$sigma, log = TRUE)
  expect_warning(
    mcp(list(y ~ 1 + x), data_gauss, family = wrong, par_x = "x", chains = 1, warmup = 300, iter = 50, diagnostics = FALSE, quiet = TRUE),
    "The JAGS likelihood and the R-side `log_lik()` of family = student() disagree",
    fixed = TRUE
  )
})


test_that("as_jags() returns the JAGS inputs that mcp() sampled", {
  data = mcp_example_data("demo")
  model = list(response ~ 1, ~ 0 + time)
  fit = mcp(model, data, chains = 2, warmup = 300, iter = 50, inits = list(Intercept_1 = 10), diagnostics = FALSE, quiet = TRUE)
  jags = as_jags(fit)
  expect_s3_class(jags, "mcpjags")
  expect_named(jags, c("code", "data", "inits", "variable.names", "modules", "levels"))
  expect_identical(jags$code, get_fit_jags_code(fit))
  expect_identical(jags$variable.names, mcp_pars(fit)$parameter)
  expect_identical(jags$modules, "glm")
  expect_length(jags$inits, 2)
  expect_equal(jags$inits[[1]]$Intercept_1, 10)
  expect_output(print(jags), "model {", fixed = TRUE)

  # Sampling these inputs reproduces the draws
  draws = run_jags(jags$code, jags$data, jags$variable.names, n.chains = 2, n.iter = 50, n.adapt = 300, inits = jags$inits, quiet = TRUE, use_glm = TRUE)
  expect_equal(as.matrix(draws)[, jags$variable.names], as.matrix(as_draws_matrix(fit))[, jags$variable.names], ignore_attr = TRUE)

  # Prior sampling and grouping levels
  prior_jags = as_jags(fit, prior = TRUE)
  expect_true(all(is.na(prior_jags$data$response)))
  expect_false(identical(prior_jags$inits[[1]]$.RNG.seed, jags$inits[[1]]$.RNG.seed))
  group_fit = mcp(list(y ~ 1 + (1 | id)), data_gauss, par_x = "x", sample = FALSE)
  expect_identical(as_jags(group_fit)$levels, list(id = c("a", "b", "c", "d", "e")))
})


test_that("jags_code is deprecated", {
  fit = mcp(list(y ~ 1 + x), data_gauss, par_x = "x", sample = FALSE)
  lifecycle::expect_deprecated(fit$jags_code)
  lifecycle::expect_deprecated(fit[["jags_code"]])
  lifecycle::expect_deprecated(mcp(list(y ~ 1 + x), data_gauss, par_x = "x", jags_code = get_fit_jags_code(fit), sample = FALSE))
})
