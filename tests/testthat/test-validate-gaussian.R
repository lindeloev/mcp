if (Sys.getenv("MCP_TEST_LEVEL") != "release") {
  testthat::skip("Time-consuming validation tests against reference implementations are only run when MCP_TEST_LEVEL='release'.")
}

# Fit a simple Gaussian model on simulated data
set.seed(42)
model = list(
  y ~ 1 + x + group
)

df = tibble::tibble(
  x = seq(0, 20, length.out = 200),
  group = rep(c("A", "B"), 100),
  y = 3 + 0.5 * x + ifelse(group == "B", -2, 0) + rnorm(200, sd = 1.5)
)

fit_mcp = mcp(model, df, family = gaussian(), warmup = 500, iter = 2000, seed = 42, diagnostics = FALSE, quiet = TRUE)

# Tests
test_that("Gaussian inference against lm()", {
  fit_lm = lm(y ~ x + group, data = df)

  # Regression coefficients
  params_mcp = fixef(fit_mcp)$mean
  params_lm = as.numeric(coef(fit_lm))
  testthat::expect_lt(max(abs(params_mcp - params_lm)), 0.1)

  capture.output({ summary_mcp = summary(fit_mcp) })
  sigma_mcp = summary_mcp$mean[summary_mcp$variable == "sigma_1"]
  testthat::expect_lt(abs(sigma_mcp - summary(fit_lm)$sigma), 0.1)

  # Log-likelihood
  loglik_mcp = mean(rowSums(log_lik(fit_mcp)))
  loglik_lm = as.numeric(logLik(fit_lm))
  expect_equal(loglik_lm, loglik_mcp, tolerance = 0.02)
})


test_that("Gaussian fixed change-point inference against lm()", {
  set.seed(42)
  df_cp = tibble::tibble(
    x = seq(0, 20, length.out = 200),
    y = 2 + 0.5 * x + 1.2 * pmax(0, x - 10) + rnorm(200, sd = 1)
  )

  fit_mcp_cp = mcp(
    list(y ~ 1 + x, ~ 0 + x),
    df_cp,
    family = gaussian(),
    prior = list(cp_1 = 10),
    warmup = 500,
    iter = 2000,
    seed = 42,
    diagnostics = FALSE,
    quiet = TRUE
  )

  fit_lm_cp = lm(y ~ x + I(pmax(0, x - 10)), data = df_cp)
  coef_lm = coef(fit_lm_cp)

  # Regression coefficients are in fixef(); the change point and residual SD
  # are available from summary().
  mcp_fixef = fixef(fit_mcp_cp)
  capture.output({ mcp_summary = summary(fit_mcp_cp) })
  estimates_mcp = c(
    mcp_fixef$mean[mcp_fixef$variable == "Intercept_1"],
    mcp_fixef$mean[mcp_fixef$variable == "x_1"],
    mcp_fixef$mean[mcp_fixef$variable == "x_2"],
    mcp_summary$mean[mcp_summary$variable == "sigma_1"]
  )

  # lm model parameterization: Intercept, x_1, x_2 (which is x_1 + delta_x), sigma
  estimates_lm = c(
    unname(coef_lm[1]),
    unname(coef_lm[2]),
    unname(coef_lm[2] + coef_lm[3]),
    summary(fit_lm_cp)$sigma
  )

  testthat::expect_lt(max(abs(estimates_mcp - estimates_lm)), 0.1)
})


test_that("Gaussian simulation against lm()", {
  newdata = df %>% dplyr::select(-"y") %>% tidyr::expand_grid(rep = 1:100)
  newdata$y = fit_mcp$simulate(fit_mcp, newdata, Intercept_1 = 3, x_1 = 0.5, groupB_1 = -2, sigma_1 = 1.5)

  fit_lm_sim = lm(y ~ x + group, data = newdata)
  params_lm_sim = c(as.numeric(coef(fit_lm_sim)), summary(fit_lm_sim)$sigma)

  testthat::expect_equal(params_lm_sim, c(3, 0.5, -2, 1.5), tolerance = 0.05)
})


test_that("Gaussian modeled sigma against nlme::gls()", {
  testthat::skip_if_not_installed("nlme")

  set.seed(42)
  N = 300
  x = seq(0, 20, length.out = N)
  seg = factor(ifelse(x < 10, "1", "2"))
  y = 2 + 0.5 * x + 1.2 * pmax(0, x - 10) + rnorm(N, sd = ifelse(seg == "1", 0.8, 2.0))
  df_sigma = data.frame(x = x, y = y, seg = seg)

  fit_gls = nlme::gls(y ~ x + I(pmax(0, x - 10)), data = df_sigma, weights = nlme::varIdent(form = ~ 1 | seg))
  coef_gls = coef(fit_gls)
  gls_sigma1 = fit_gls$sigma
  gls_sigma2 = gls_sigma1 * coef(fit_gls$modelStruct$varStruct, unconstrained = FALSE)[["2"]]

  fit_mcp_sigma = mcp(
    list(y ~ 1 + x, ~ 0 + x + sigma(1)),
    df_sigma,
    prior = list(cp_1 = 10),
    warmup = 500,
    iter = 2000,
    seed = 42,
    diagnostics = FALSE,
    quiet = TRUE
  )

  capture.output({ sum_mcp = summary(fit_mcp_sigma) })
  mcp_sigma1 = exp(sum_mcp$mean[sum_mcp$variable == "sigma_1"])
  mcp_sigma2 = exp(sum_mcp$mean[sum_mcp$variable == "sigma_2"])

  expect_equal(sum_mcp$mean[sum_mcp$variable == "Intercept_1"], unname(coef_gls[1]), tolerance = 0.1)
  expect_equal(sum_mcp$mean[sum_mcp$variable == "x_1"], unname(coef_gls[2]), tolerance = 0.05)
  expect_equal(sum_mcp$mean[sum_mcp$variable == "x_2"], unname(coef_gls[2] + coef_gls[3]), tolerance = 0.05)
  expect_equal(mcp_sigma1, gls_sigma1, tolerance = 0.05)
  expect_equal(mcp_sigma2, gls_sigma2, tolerance = 0.05)
})


test_that("Gaussian likelihood weights against lm(weights = w)", {
  set.seed(42)
  N = 250
  x = seq(0, 10, length.out = N)
  w = runif(N, 0.5, 3.0)
  y = 1.5 + 0.8 * x + rnorm(N, sd = 1.2 / sqrt(w))
  df_w = data.frame(x = x, y = y, w = w)

  fit_lm_w = lm(y ~ x, data = df_w, weights = w)
  coef_lm_w = coef(fit_lm_w)
  sigma_lm_weighted = summary(fit_lm_w)$sigma * sqrt((N - 2) / sum(w))

  fit_mcp_w = mcp(
    list(y | weights(w) ~ 1 + x),
    data = df_w,
    family = gaussian(),
    warmup = 500,
    iter = 2000,
    seed = 42,
    diagnostics = FALSE,
    quiet = TRUE
  )

  fix_mcp_w = fixef(fit_mcp_w)
  capture.output({ sum_mcp_w = summary(fit_mcp_w) })
  sigma_mcp_w = sum_mcp_w$mean[sum_mcp_w$variable == "sigma_1"]

  expect_equal(fix_mcp_w$mean[fix_mcp_w$variable == "Intercept_1"], unname(coef_lm_w[1]), tolerance = 0.05)
  expect_equal(fix_mcp_w$mean[fix_mcp_w$variable == "x_1"], unname(coef_lm_w[2]), tolerance = 0.02)
  expect_equal(sigma_mcp_w, sigma_lm_weighted, tolerance = 0.05)

  loglik_expected = sum(w * dnorm(df_w$y, predict(fit_lm_w), sigma_mcp_w, log = TRUE))
  loglik_mcp = mean(rowSums(log_lik(fit_mcp_w)))
  expect_equal(loglik_mcp, loglik_expected, tolerance = 0.02)
})
