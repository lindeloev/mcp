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

  expect_match(fit$jags_code, "likelihood_weight_[i_] = 1 + response_observed_[i_] * (w[i_] - 1)", fixed = TRUE)
  expect_match(fit$jags_code, "likelihood_zero_[i_] ~ dexp(exp(max(-700, (likelihood_weight_[i_] - 1) * (loggam(y[i_] + shape_[i_]) - loggam(shape_[i_]) - loggam(y[i_] + 1) + shape_[i_] * log(nb_prob_[i_]) + y[i_] * log(1 - nb_prob_[i_])))))", fixed = TRUE)

  # Check that package R-side log_lik applies the weights directly
  expect_equal(
    fit$family$r$log_lik(df_nb$y, list(mu = rep(200, 5), shape = 10), list(weights = df_nb$w)),
    df_nb$w * stats::dnbinom(df_nb$y, mu = 200, size = 10, log = TRUE)
  )
})



