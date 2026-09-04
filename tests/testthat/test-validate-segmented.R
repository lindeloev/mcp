if (Sys.getenv("MCP_TEST_LEVEL") != "release") {
  testthat::skip("Time-consuming validation tests against reference implementations are only run when MCP_TEST_LEVEL='release'.")
}

# Compare mcp estimated change points and segment parameters against segmented::segmented()

test_that("Gaussian 2-changepoint model agrees with segmented()", {
  testthat::skip_if_not_installed("segmented")

  set.seed(42)
  x = 1:100
  y = 2 + 0.5 * x - 1.2 * pmax(0, x - 30) + 1.5 * pmax(0, x - 70) + rnorm(100, sd = 1.2)
  df = data.frame(x = x, y = y)

  # Fit with segmented
  fit_lm = lm(y ~ x, data = df)
  fit_seg = segmented::segmented(fit_lm, seg.Z = ~x, npsi = 2)

  # Fit with mcp (2 joined change points)
  model = list(
    y ~ 1 + x,
    ~ 0 + x,
    ~ 0 + x
  )
  fit_mcp = mcp(model, df, family = gaussian(), warmup = 500, iter = 2000, seed = 42, diagnostics = FALSE, quiet = TRUE)

  # Extract mcp estimates
  fix_mcp = fixef(fit_mcp)
  capture.output({ summary_mcp = summary(fit_mcp) })
  cp1_mcp = summary_mcp$mean[summary_mcp$variable == "cp_1"]
  cp2_mcp = summary_mcp$mean[summary_mcp$variable == "cp_2"]
  intercept_mcp = fix_mcp$mean[fix_mcp$variable == "Intercept_1"]
  x1_mcp = fix_mcp$mean[fix_mcp$variable == "x_1"]
  x2_mcp = fix_mcp$mean[fix_mcp$variable == "x_2"]
  x3_mcp = fix_mcp$mean[fix_mcp$variable == "x_3"]
  sigma_mcp = summary_mcp$mean[summary_mcp$variable == "sigma_1"]

  # Extract segmented estimates
  psi_seg = fit_seg$psi[, "Est."]
  slopes_seg = segmented::slope(fit_seg)$x[, "Est."]
  intercept_seg = unname(coef(fit_seg)["(Intercept)"])
  sigma_seg = summary(fit_seg)$sigma

  # Compare change points (within 0.5 units on x)
  testthat::expect_equal(cp1_mcp, unname(psi_seg["psi1.x"]), tolerance = 0.5)
  testthat::expect_equal(cp2_mcp, unname(psi_seg["psi2.x"]), tolerance = 0.5)

  # Compare parameters (within 0.05, except intercept posterior mean vs. plug-in MLE)
  testthat::expect_equal(intercept_mcp, intercept_seg, tolerance = 0.15)
  testthat::expect_equal(x1_mcp, unname(slopes_seg["slope1"]), tolerance = 0.05)
  testthat::expect_equal(x2_mcp, unname(slopes_seg["slope2"]), tolerance = 0.05)
  testthat::expect_equal(x3_mcp, unname(slopes_seg["slope3"]), tolerance = 0.05)
  testthat::expect_equal(sigma_mcp, sigma_seg, tolerance = 0.05)
})


test_that("Poisson 2-changepoint model agrees with segmented()", {
  testthat::skip_if_not_installed("segmented")

  set.seed(123)
  x = 1:100
  log_mu = 1 + 0.04 * x - 0.07 * pmax(0, x - 35) + 0.08 * pmax(0, x - 75)
  y = rpois(100, lambda = exp(log_mu))
  df = data.frame(x = x, y = y)

  # Fit with segmented
  fit_glm = glm(y ~ x, data = df, family = poisson())
  fit_seg = segmented::segmented(fit_glm, seg.Z = ~x, npsi = 2)

  # Fit with mcp
  model = list(
    y ~ 1 + x,
    ~ 0 + x,
    ~ 0 + x
  )
  fit_mcp = mcp(model, df, family = poisson(), warmup = 500, iter = 2000, seed = 123, diagnostics = FALSE, quiet = TRUE)

  # Extract mcp estimates
  fix_mcp = fixef(fit_mcp)
  capture.output({ summary_mcp = summary(fit_mcp) })
  cp1_mcp = summary_mcp$mean[summary_mcp$variable == "cp_1"]
  cp2_mcp = summary_mcp$mean[summary_mcp$variable == "cp_2"]
  intercept_mcp = fix_mcp$mean[fix_mcp$variable == "Intercept_1"]
  x1_mcp = fix_mcp$mean[fix_mcp$variable == "x_1"]
  x2_mcp = fix_mcp$mean[fix_mcp$variable == "x_2"]
  x3_mcp = fix_mcp$mean[fix_mcp$variable == "x_3"]

  # Extract segmented estimates
  psi_seg = fit_seg$psi[, "Est."]
  slopes_seg = segmented::slope(fit_seg)$x[, "Est."]
  intercept_seg = unname(coef(fit_seg)["(Intercept)"])

  # Compare change points (within 2 units due to Poisson sampling variation)
  testthat::expect_equal(cp1_mcp, unname(psi_seg["psi1.x"]), tolerance = 2.0)
  testthat::expect_equal(cp2_mcp, unname(psi_seg["psi2.x"]), tolerance = 2.0)

  # Compare parameters (within 0.05)
  testthat::expect_equal(intercept_mcp, intercept_seg, tolerance = 0.05)
  testthat::expect_equal(x1_mcp, unname(slopes_seg["slope1"]), tolerance = 0.05)
  testthat::expect_equal(x2_mcp, unname(slopes_seg["slope2"]), tolerance = 0.05)
  testthat::expect_equal(x3_mcp, unname(slopes_seg["slope3"]), tolerance = 0.05)
})


test_that("Binomial 2-changepoint model agrees with segmented()", {
  testthat::skip_if_not_installed("segmented")

  set.seed(42)
  N = round(runif(100, 10, 20))
  x = 1:100
  logit_p = -1 + 0.08 * x - 0.12 * pmax(0, x - 35) + 0.10 * pmax(0, x - 75)
  y = rbinom(100, N, stats::plogis(logit_p))
  df = data.frame(x = x, y = y, N = N)

  # Fit with segmented
  fit_glm = glm(cbind(y, N - y) ~ x, data = df, family = binomial())
  fit_seg = segmented::segmented(fit_glm, seg.Z = ~x, npsi = 2)

  # Fit with mcp
  model = list(
    y | trials(N) ~ 1 + x,
    ~ 0 + x,
    ~ 0 + x
  )
  fit_mcp = mcp(model, df, family = binomial(), warmup = 500, iter = 2000, seed = 42, diagnostics = FALSE, quiet = TRUE)

  # Extract mcp estimates
  fix_mcp = fixef(fit_mcp)
  capture.output({ summary_mcp = summary(fit_mcp) })
  cp1_mcp = summary_mcp$mean[summary_mcp$variable == "cp_1"]
  cp2_mcp = summary_mcp$mean[summary_mcp$variable == "cp_2"]
  intercept_mcp = fix_mcp$mean[fix_mcp$variable == "Intercept_1"]
  x1_mcp = fix_mcp$mean[fix_mcp$variable == "x_1"]
  x2_mcp = fix_mcp$mean[fix_mcp$variable == "x_2"]
  x3_mcp = fix_mcp$mean[fix_mcp$variable == "x_3"]

  # Extract segmented estimates
  psi_seg = fit_seg$psi[, "Est."]
  slopes_seg = segmented::slope(fit_seg)$x[, "Est."]
  intercept_seg = unname(coef(fit_seg)["(Intercept)"])

  # Compare change points (within 5.0 units given discrete trials binomial noise)
  testthat::expect_equal(cp1_mcp, unname(psi_seg["psi1.x"]), tolerance = 5.0)
  testthat::expect_equal(cp2_mcp, unname(psi_seg["psi2.x"]), tolerance = 8.0)

  # Compare parameters (within 0.15 given Bayesian regularizing priors vs MLE)
  testthat::expect_equal(intercept_mcp, intercept_seg, tolerance = 0.15)
  testthat::expect_equal(x1_mcp, unname(slopes_seg["slope1"]), tolerance = 0.15)
})
