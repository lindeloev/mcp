bad_poisson = list(
  # Misspecification of y and trials
  list(y | trials(N) ~ 1),  # bad response format
  list(y ~ 1 + x,
       y | trials(N) ~ 1 ~ 1),  # misspecification in later segment

  # Bad data
  list(y_bad_numeric ~ 1),

  # Does not work with sigma
  list(y ~ 1 + sigma(1)),

  # Coefficient-free model
  list(y ~ 0)
)

test_bad(bad_poisson,
         data = data_binomial,
         family = poisson())


good_poisson = list(
  list(y ~ 1),  # one segment
  list(y ~ 1 + x,  # specified multiple times
       y  ~ 1 ~ 1 + x,
       1 ~ 0),
  list(y ~ 1,  # With varying
       1 + (1|id) ~ 1),
  list(y ~ 1 + ar(1) + ma(1)),
  list(y | weights(weights_ok) ~ 1)  # With weights
)

test_good(good_poisson,
          data = data_binomial,
          family = poisson())


test_that("Poisson JAGS weights implement a likelihood power and sample on ordinary large counts", {
  # Test reproduction failure case: count of 200 with weight of 2
  df_pois = data.frame(x = 1:5, y = rep(200, 5), w = 2)
  fit = mcp(list(y | weights(w) ~ 1), data = df_pois, family = poisson(), par_x = "x", quiet = TRUE)
  expect_true(is.list(fit$model))
  capture.output({ summary_fit = summary(fit) })
  expect_equal(nrow(summary_fit), 1)
  expect_true(abs(summary_fit$mean[1] - log(200)) < 0.1)

  expect_match(fit$jags_code, "likelihood_phi_[i_] = response_observed_[i_] * w[i_] * max(0, mu_[i_] - y[i_] * log(mu_[i_]) + loggam(y[i_] + 1))", fixed = TRUE)
  expect_match(fit$jags_code, "likelihood_zero_[i_] ~ dpois(likelihood_phi_[i_])", fixed = TRUE)

  # Check that package R-side log_lik applies the weights directly
  expect_equal(
    fit$family$r$log_lik(df_pois$y, list(mu = rep(200, 5)), list(weights = df_pois$w)),
    df_pois$w * stats::dpois(df_pois$y, lambda = 200, log = TRUE)
  )
})


test_that("Weighted Poisson targets exact conjugate posterior without clamping", {
  data_conj = data.frame(x = 1:3, y = 10, w = 1000)
  fit = mcp(
    list(y | weights(w) ~ 1),
    data = data_conj,
    family = poisson("identity"),
    prior = list(Intercept_1 = "dgamma(1, 0.1)"),
    par_x = "x",
    quiet = TRUE
  )
  draws = as.matrix(coda::as.mcmc(fit))[, "Intercept_1"]
  expect_equal(mean(draws), 10, tolerance = 0.05)
  # Correct posterior is Gamma(shape = 30001, rate = 3000.1) with SD = 0.0577.
  # Before the fix, the clamp at -700 caused SD to be ~ 1.82.
  expect_equal(stats::sd(draws), 0.0577, tolerance = 0.005)

  # Check that weights with missing responses errors cleanly
  data_na = data.frame(x = 1:3, y = c(10, NA, 10), w = 10)
  expect_error(
    mcp(list(y | weights(w) ~ 1), data = data_na, family = poisson(), par_x = "x"),
    "Weights with missing responses"
  )
})


test_that("Coefficient-free models error with an informative message", {
  expect_error(
    mcp(list(y ~ 0), data.frame(x = 1:6, y = 1:6), family = poisson(), par_x = "x", sample = FALSE),
    "The model does not contain any parameters to estimate."
  )
})

