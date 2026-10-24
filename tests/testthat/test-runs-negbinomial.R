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
  fit = mcp(list(y | weights(w) ~ 1), data = df_nb, family = negbinomial(), par_x = "x", diagnostics = FALSE, quiet = TRUE, seed = 42)  # disabling diagnostics is justified because this tests and edge case - not the fit
  expect_true(is.list(fit$model))
  capture.output({ summary_fit = summary(fit) })
  expect_equal(nrow(summary_fit), 2)
  expect_true(abs(summary_fit$mean[summary_fit$parameter == "Intercept_1"] - log(200)) < 0.1)

  expect_match(fit$jags_code, "likelihood_phi_[i_] = response_observed_[i_] * w[i_] * max(0, loggam(shape_[i_]) + loggam(y[i_] + 1) - loggam(y[i_] + shape_[i_]) - shape_[i_] * log(nb_prob_[i_]) - y[i_] * log(1 - nb_prob_[i_]))", fixed = TRUE)
  expect_match(fit$jags_code, "likelihood_zero_[i_] ~ dpois(likelihood_phi_[i_])", fixed = TRUE)

  # Check that R-side log_lik() applies the weights directly
  draws = as.data.frame(fit)
  expect_equal(
    unname(log_lik(fit, summary = FALSE)),
    2 * stats::dnbinom(matrix(df_nb$y, nrow(draws), 5, byrow = TRUE), mu = exp(draws$Intercept_1), size = draws$shape_1, log = TRUE)
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
    shape_1 = shape
  )

  loglik = evaluate_model(
    fit, get_model_data(fit, data), list(Intercept_1 = log(mu), x_1 = 0, shape_1 = shape),
    n_draws = 1, y = data$y, type = "loglik"
  )$value
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

  expect_equal(fit$family$dpar_specs$dpar, c("mu", "shape"))
  expect_equal(fit$family$dpar_specs$link, c("log", "identity"))
  expect_equal(fit$prior$shape_1, "dinvgamma(0.4, 0.3)")
  expect_match(fit$jags_code, "shape_1 = 1 / shape_1_inverse", fixed = TRUE)

  # User priors on the constant shape are bounded to its positive support
  fit_user = mcp(list(y ~ 1 + x), data, family = negbinomial(), prior = list(shape_1 = "dnorm(10, 2)"), sample = FALSE)
  expect_equal(fit_user$prior$shape_1, "dnorm(10, 2) T(0, )")
  fit_user = mcp(list(y ~ 1 + x), data, family = negbinomial(), prior = list(shape_1 = "dinvgamma(1, 2)"), sample = FALSE)
  expect_equal(fit_user$prior$shape_1, "dinvgamma(1, 2)")
})


test_that("modeled shape uses link-scale coefficient priors", {
  data = data.frame(x = 1:5, y = c(0, 1, 2, 4, 8))
  fit = mcp(
    list(y ~ 1 + x + shape(1 + x)),
    data,
    family = negbinomial(),
    sample = FALSE
  )

  expect_equal(fit$family$dpar_specs$link, c("log", "log"))
  expect_equal(fit$prior$shape_Intercept_1, "dnorm(-shape_x_1 * 1, 2.5)")
  expect_prior(
    fit$prior$shape_x_1,
    paste0("dnorm(0, ", format_prior_number(2.5 / stats::sd(data$x)), ")")
  )
})


test_that("negative-binomial support does not alter Poisson metadata or priors", {
  data = data.frame(x = 1:5, y = c(0, 1, 2, 4, 8))
  fit = mcp(list(y ~ 1 + x), data, family = poisson(), sample = FALSE)

  expect_equal(fit$family$dpar_specs$dpar, "mu")
  expect_equal(fit$family$dpar_specs$link, "log")
  expect_equal(fit$prior$Intercept_1, "dnorm(0.7 - x_1 * 1, 2.5)")
  expect_prior(fit$prior$x_1, paste0("dnorm(0, ", format_prior_number(2.5 / stats::sd(data$x)), ")"))

  negbin_fit = mcp(list(y ~ 1 + x), data, family = negbinomial(), sample = FALSE)
  expect_equal(negbin_fit$prior$Intercept_1, fit$prior$Intercept_1)
  expect_equal(negbin_fit$prior$x_1, fit$prior$x_1)

  identity_fit = mcp(
    list(y ~ 1 + x), data,
    family = poisson(link = "identity"), sample = FALSE
  )
  expect_prior(identity_fit$prior$Intercept_1, paste0("dnorm(", format_prior_number(mean(data$y)), " - x_1 * 1, ", format_prior_number(2.5 * stats::sd(data$y)), ") T(-x_1 * 1, )"))
  expect_prior(identity_fit$prior$x_1, paste0("dnorm(0, ", format_prior_number(2.5 * stats::sd(data$y) / stats::sd(data$x)), ")"))
})

