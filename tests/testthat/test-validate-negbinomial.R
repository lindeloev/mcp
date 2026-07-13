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

  expect_equal(fit$family$dpars, c("mu", "shape", "ar", "ma"))
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

  expect_equal(fit$prior$shape_1, "dt(0, 2.5, 3)")
  expect_equal(
    fit$prior$shape_x_1,
    "dt(0, (N_CP + 1) * 2.5 / (MAXX - MINX), 3)"
  )
})


test_that("negative-binomial support does not alter Poisson metadata or priors", {
  data = data.frame(x = 1:5, y = c(0, 1, 2, 4, 8))
  fit = mcp(list(y ~ 1 + x), data, family = poisson(), sample = FALSE)

  expect_equal(fit$family$dpars, c("mu", "ar", "ma"))
  expect_equal(fit$family$links, c(mu = "log"))
  expect_equal(fit$prior$Intercept_1, "dt(0, 10, 3)")
  expect_equal(fit$prior$x_1, "dt(0, 10, 3)")
})


test_that("negative-binomial coefficients agree with MASS glm.nb", {
  testthat::skip_if_not_installed("MASS")
  set.seed(42)
  data = data.frame(x = seq(-1, 1, length.out = 300))
  data$y = stats::rnbinom(
    nrow(data),
    mu = exp(1 + 0.6 * data$x),
    size = 2.5
  )

  fit_mcp = quiet_mcp(
    list(y ~ 1 + x),
    data,
    family = negbinomial(),
    adapt = 500,
    iter = 2000,
    chains = 2
  )
  fit_mass = MASS::glm.nb(y ~ x, data = data)

  samples = as.matrix(fit_mcp$mcmc_post)
  expect_equal(mean(samples[, "Intercept_1"]), unname(stats::coef(fit_mass)[1]), tolerance = 0.1)
  expect_equal(mean(samples[, "x_1"]), unname(stats::coef(fit_mass)[2]), tolerance = 0.1)
  expect_equal(stats::median(exp(samples[, "shape_1"])), fit_mass$theta, tolerance = 0.5)
})
