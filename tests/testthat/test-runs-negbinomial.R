bad_negbinomial = list(
  list(y | trials(N) ~ 1),
  list(y ~ 1 + sigma(1)),
  list(y_bad_numeric ~ 1)
)

test_bad(
  bad_negbinomial,
  data = data_binomial,
  family = negbinomial()
)


good_negbinomial = list(
  list(y ~ 1 + x),
  list(
    y ~ 1 + x,
    ~ 1 + x
  ),
  list(y ~ 1 + x + shape(1 + x)),
  list(y ~ 1 + shape(1 + (1 | id))),
  list(
    y ~ 1 + x + shape(1),
    ~ 0 + x + shape(1)
  ),
  list(y ~ 1 + ar(1) + ma(1)),
  list(y | weights(weights_ok) ~ 1)  # With weights
)

test_good(
  good_negbinomial,
  data = data_binomial,
  family = negbinomial()
)


test_that("negative-binomial links are explicit and currently log-only", {
  expect_error(negbinomial(link = "identity"), '`link` must be one of "log"')
  expect_error(negbinomial(link_shape = "identity"), '`link_shape` must be one of "log"')
})


test_that("Negative-binomial JAGS weights implement a likelihood power and sample on ordinary large counts", {
  # Test reproduction failure case: count of 200, shape 10, weight 2
  df_nb = data.frame(x = 1:5, y = rep(200, 5), w = 2)
  fit = mcp(list(y | weights(w) ~ 1), data = df_nb, family = negbinomial(), par_x = "x", diagnostics = FALSE, quiet = TRUE)  # disabling diagnostics is justified because this tests and edge case - not the fit
  expect_true(is.list(fit$model))
  capture.output({ summary_fit = summary(fit) })
  expect_equal(nrow(summary_fit), 2)
  expect_true(abs(summary_fit$mean[summary_fit$variable == "Intercept_1"] - log(200)) < 0.1)

  expect_match(fit$jags_code, "likelihood_phi_[i_] = response_observed_[i_] * w[i_] * max(0, loggam(shape_[i_]) + loggam(y[i_] + 1) - loggam(y[i_] + shape_[i_]) - shape_[i_] * log(nb_prob_[i_]) - y[i_] * log(1 - nb_prob_[i_]))", fixed = TRUE)
  expect_match(fit$jags_code, "likelihood_zero_[i_] ~ dpois(likelihood_phi_[i_])", fixed = TRUE)

  # Check that package R-side log_lik applies the weights directly
  expect_equal(
    fit$family$r$log_lik(df_nb$y, list(mu = rep(200, 5), shape = 10), list(weights = df_nb$w)),
    df_nb$w * stats::dnbinom(df_nb$y, mu = 200, size = 10, log = TRUE)
  )
})


test_that("Weighted Negative-Binomial targets exact posterior without clamping", {
  df_nb_large = data.frame(x = 1:5, y = rep(20, 5), w = 100)
  fit_nb = mcp(
    list(y | weights(w) ~ 1),
    data = df_nb_large,
    family = negbinomial(),
    prior = list(Intercept_1 = "dnorm(3, 1)", shape_1 = "dnorm(10, 2)"),
    par_x = "x",
    diagnostics = FALSE,
    quiet = TRUE
  )
  draws = as.matrix(coda::as.mcmc(fit_nb))[, "Intercept_1"]
  expect_equal(mean(draws), log(20), tolerance = 0.05)
})


test_that("negative-binomial R evaluation matches stats", {
  data = data.frame(x = seq(-1, 1, length.out = 10000), y = rep(0, 10000))
  fit = mcp(
    list(y ~ 1 + x),
    data,
    family = negbinomial(),
    sample = FALSE
  )

  mu = 4
  shape = 2
  args = list(
    fit,
    data,
    Intercept_1 = log(mu),
    x_1 = 0,
    shape_1 = log(shape)
  )

  predictors = add_rhs_predictors(data, fit)
  loglik = rlang::exec(
    simulate_vectorized,
    fit,
    !!!predictors,
    Intercept_1 = rep(log(mu), nrow(data)),
    x_1 = rep(0, nrow(data)),
    shape_1 = rep(log(shape), nrow(data)),
    .type = "loglik"
  )
  expect_equal(
    as.numeric(loglik),
    stats::dnbinom(data$y, mu = mu, size = shape, log = TRUE)
  )

  set.seed(123)
  prediction = rlang::exec(fit$simulate, !!!args, .type = "predict")
  expect_equal(mean(prediction), mu, tolerance = 0.1)
  expect_equal(stats::var(prediction), mu + mu^2 / shape, tolerance = 0.5)
})


test_that("negative-binomial exposes mean-shape metadata and defaults", {
  data = data.frame(x = 1:5, y = c(0, 1, 2, 4, 8))
  fit = mcp(list(y ~ 1 + x), data, family = negbinomial(), sample = FALSE)

  expect_equal(fit$family$dpars, c("mu", "shape"))
  expect_equal(fit$family$links, c(mu = "log", shape = "log"))
  expect_equal(fit$prior$shape_1, "dloginvgamma(0.4, 0.3)")
})


test_that("modeled shape uses link-scale coefficient priors", {
  data = data.frame(x = 1:5, y = c(0, 1, 2, 4, 8))
  fit = mcp(
    list(y ~ 1 + x + shape(1 + x)),
    data,
    family = negbinomial(),
    sample = FALSE
  )

  expect_equal(fit$prior$shape_1, "dnorm(0, 2.5)")
  expect_equal(
    fit$prior$shape_x_1,
    "dnorm(0, 0.625)"
  )
})


test_that("negative-binomial support does not alter Poisson metadata or priors", {
  data = data.frame(x = 1:5, y = c(0, 1, 2, 4, 8))
  fit = mcp(list(y ~ 1 + x), data, family = poisson(), sample = FALSE)

  expect_equal(fit$family$dpars, "mu")
  expect_equal(fit$family$links, c(mu = "log"))
  expect_equal(fit$prior$Intercept_1, "dnorm(0.7, 2.5)")
  expect_equal(fit$prior$x_1, "dnorm(0, 0.625)")

  negbin_fit = mcp(list(y ~ 1 + x), data, family = negbinomial(), sample = FALSE)
  expect_equal(negbin_fit$prior$Intercept_1, fit$prior$Intercept_1)
  expect_equal(negbin_fit$prior$x_1, fit$prior$x_1)

  identity_fit = mcp(
    list(y ~ 1 + x), data,
    family = poisson(link = "identity"), sample = FALSE
  )
  expect_equal(identity_fit$prior$Intercept_1, "dt(2, 3, 3) T(0, )")
  expect_equal(identity_fit$prior$x_1, "dt(0, 0.75, 3)")
})

