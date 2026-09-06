#################
# TEST BINOMIAL #
#################

bad_binomial = list(
  # Misspecification of y and trials
  list(y ~ 1),  # no trials
  list(y | N ~ 1),  # wrong format
  list(trials(N) | y ~ 1),  # Wrong order
  list(y | trials() ~ 1),  # trials missing
  list(trials(N) ~ 1),  # no y
  list(y | trials(N) ~ 1 + x,
       y | N ~ 1 ~ 1),  # misspecification in later segment

  # Bad data
  list(y_bad_numeric | trials(N) ~ 1),
  list(y | trials(N_bad_numeric) ~ 1),
  list(y | trials(N_bad_factor) ~ 1),
  list(y | trials(N_bad_char) ~ 1),

  # Does not work with sigma
  list(y | trials(N) ~ 1 + sigma(1))
)

test_bad(bad_binomial,
         data = data_binomial,
         family = binomial())


good_binomial = list(
  list(y | trials(N) ~ 1),  # one segment
  list(y | trials(N) ~ 1 + x,  # specified multiple times
       y | trials(N) ~ 1 ~ 1 + x,
       ~ 0),
  list(y | trials(N) ~ 1,  # With varying
       1 + (1|id) ~ 1),
  list(y | trials(N) ~ 1 + ar(1) + ma(1)),
  list(y | trials(N) ~ 1,
       1 ~ N),  # N can be both trials and slope
  list(y | trials(N) + weights(weights_ok) ~ 1)  # With weights
)

test_good(good_binomial,
          data = data_binomial,
          family = binomial())

test_that("binomial and Bernoulli links share weakly regularizing defaults", {
  families = list(
    mcpfamily(binomial("logit")),
    mcpfamily(binomial("probit")),
    bernoulli("logit"),
    bernoulli("probit")
  )
  expected = c(
    "dt(0, 1.5, 3)",
    "dt(0, 1.5, 3)",
    "dt(0, 1.5 / predictor_scale(), 3)"
  )

  for (family in families)
    expect_equal(family$default_prior$prior, expected)
})

test_that("binomial responses cannot exceed trials", {
  invalid_data = data_binomial
  invalid_data$y[1] = invalid_data$N[1] + 1

  expect_error(
    mcp(
      list(y | trials(N) ~ 1),
      data = invalid_data,
      family = binomial(),
      par_x = "x",
      sample = FALSE
    ),
    "responses in 'y' cannot exceed trials in 'N'. Found invalid data in row(s): 1.",
    fixed = TRUE
  )
})




test_that("binomial prediction defaults agree with posterior methods and warn on migration", {
  withr::local_options(rlib_warning_verbosity = "verbose")
  data = data.frame(x = 1:3, N = c(4, 10, 6), y = 0)
  fit = mcp(list(y | trials(N) ~ 1), data, family = binomial(), par_x = "x", sample = FALSE)
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(
    matrix(0, nrow = 20, ncol = 1, dimnames = list(NULL, "Intercept_1"))
  ))

  expect_warning((counts = fitted(fit)), "Set `rate = TRUE`")
  expect_equal(counts$fitted, data$N / 2)
  expect_warning(predict(fit), "Set `rate = TRUE`")
  expect_no_warning(fitted(fit, rate = FALSE))
  expect_no_warning(fitted(fit, rate = TRUE))
  expect_no_warning(predict(fit, rate = FALSE))
  expect_no_warning(predict(fit, rate = TRUE))
  expect_no_warning(fitted(fit, dpar = "mu"))
  expect_no_warning(fitted(fit, scale = "linear"))
  expect_no_warning(fitted(fit, newdata = transform(data, N = 1)))
  expect_no_warning(predict(fit, newdata = transform(data, N = 1)))
  expect_warning(fitted(fit, dpar = NULL), "Set `rate = TRUE`")
  expect_equal(residuals(fit)$residuals, data$y - counts$fitted)

  for (rate in c(FALSE, TRUE)) {
    expected = fitted(fit, rate = rate, summary = FALSE, draws_format = "matrix")
    expect_no_warning((actual = posterior_epred.mcpfit(fit, rate = rate)))
    expect_equal(actual, expected)
    expect_equal(unname(actual[1, ]), if (rate) rep(0.5, 3) else data$N / 2)

    set.seed(42)
    expected = predict(fit, rate = rate, summary = FALSE, draws_format = "matrix")
    expect_no_warning((actual = posterior_predict.mcpfit(fit, rate = rate, seed = 42)))
    expect_equal(actual, expected)
  }
})


##################
# TEST BERNOULLI #
##################
# This is rather short since most is tested via binomial
bad_bernoulli = list(
  # Misspecification of y and trials
  list(y_bern | trials(N) ~ 1),  # trials
  list(y_bern ~ 1 + x,
       y_bern | trials(N) ~ 1 ~ 1),  # misspecification in later segment

  # Bad data
  list(y_bad_numeric ~ 1),
  list(y ~ 1),  # binomial response

  # Does not work with sigma
  list(y_bern ~ 1 + sigma(1)),

  # Bernoulli does not take trials
  list(y | trials(N) + weights(weights_ok) ~ 1)
)

test_bad(bad_bernoulli,
         data = data_binomial,
         family = bernoulli())


good_bernoulli = list(
  list(y_bern ~ 1),  # one segment
  list(y_bern ~ 1 + x,  # specified multiple times
       y_bern ~ 1 ~ 1 + x,
       1 ~ 0),
  list(y_bern ~ 1,  # With varying
       1 + (1|id) ~ 1),
  list(y_bern | weights(weights_ok) ~ 1),  # With weights
  list(y_bern ~ 1 + ar(1) + ma(1))  # With AR and MA
)

test_good(good_bernoulli,
          data = data_binomial,
          family = bernoulli())


test_that("Binomial JAGS weights implement a likelihood power and sample on ordinary large counts", {
  # Test reproduction failure case: 1000 successes out of 2000 trials with weight of 2
  df_bin = data.frame(x = 1:5, y = rep(1000, 5), N = 2000, w = 2)
  fit = mcp(list(y | trials(N) + weights(w) ~ 1), data = df_bin, family = binomial(), par_x = "x")
  expect_true(is.list(fit$model))
  expect_equal(nrow(summary(fit)), 1)
  expect_true(abs(summary(fit)$mean[1] - 0) < 0.1)

  expect_match(fit$jags_code, "likelihood_weight_[i_] = 1 + response_observed_[i_] * (w[i_] - 1)", fixed = TRUE)
  expect_match(fit$jags_code, "likelihood_zero_[i_] ~ dexp(exp(max(-700, (likelihood_weight_[i_] - 1) * (loggam(N[i_] + 1) - loggam(y[i_] + 1) - loggam(N[i_] - y[i_] + 1) + y[i_] * log(mu_[i_]) + (N[i_] - y[i_]) * log(1 - mu_[i_])))))", fixed = TRUE)

  # Check that package R-side log_lik applies the weights directly
  expect_equal(
    fit$family$r$log_lik(df_bin$y, list(mu = rep(0.5, 5)), list(trials = df_bin$N, weights = df_bin$w)),
    df_bin$w * stats::dbinom(df_bin$y, size = df_bin$N, prob = 0.5, log = TRUE)
  )
})
