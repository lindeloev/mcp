if (Sys.getenv("MCP_TEST_LEVEL") != "release") {
  testthat::skip("Time-consuming validation tests against reference implementations are only run when MCP_TEST_LEVEL='release'.")
}

test_that("negative-binomial coefficients agree with MASS glm.nb", {
  testthat::skip_if_not_installed("MASS")
  set.seed(42)
  data = data.frame(x = seq(-1, 1, length.out = 300))
  data$y = stats::rnbinom(
    nrow(data),
    mu = exp(1 + 0.6 * data$x),
    size = 2.5
  )

  fit_mcp = mcp(
    list(y ~ 1 + x),
    data,
    family = negbinomial(),
    warmup = 500,
    iter = 2000,
    chains = 2,
    diagnostics = FALSE,
    quiet = TRUE
  )
  fit_mass = MASS::glm.nb(y ~ x, data = data)

  samples = posterior::as_draws_matrix(fit_mcp)
  expect_equal(mean(samples[, "Intercept_1"]), unname(stats::coef(fit_mass)[1]), tolerance = 0.1)
  expect_equal(mean(samples[, "x_1"]), unname(stats::coef(fit_mass)[2]), tolerance = 0.1)
  expect_equal(stats::median(exp(samples[, "shape_1"])), fit_mass$theta, tolerance = 0.5)
})

