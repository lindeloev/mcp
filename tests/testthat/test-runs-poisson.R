bad_poisson = list(
  # Misspecification of y and trials
  list(y | trials(N) ~ 1),  # bad response format
  list(y ~ 1 + x,
       y | trials(N) ~ 1 ~ 1),  # misspecification in later segment

  # Bad data
  list(y_bad_numeric ~ 1),

  # Does not work with sigma
  list(y ~ 1 + sigma(1))
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

  expect_match(fit$jags_code, "likelihood_weight_[i_] = 1 + response_observed_[i_] * (w[i_] - 1)", fixed = TRUE)
  expect_match(fit$jags_code, "likelihood_zero_[i_] ~ dexp(exp(max(-700, (likelihood_weight_[i_] - 1) * (y[i_] * log(mu_[i_]) - mu_[i_] - loggam(y[i_] + 1)))))", fixed = TRUE)

  # Check that package R-side log_lik applies the weights directly
  expect_equal(
    fit$family$r$log_lik(df_pois$y, list(mu = rep(200, 5)), list(weights = df_pois$w)),
    df_pois$w * stats::dpois(df_pois$y, lambda = 200, log = TRUE)
  )
})

