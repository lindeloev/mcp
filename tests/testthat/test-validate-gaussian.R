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

fit_mcp = mcp(model, df, family = gaussian(), adapt = 500, iter = 2000, seed = 42, diagnostics = FALSE, quiet = TRUE)

# Tests
test_that("Gaussian inference against lm()", {
  fit_lm = lm(y ~ x + group, data = df)

  # Parameter estimates: (Intercept, x, groupB, sigma)
  params_mcp = fixef(fit_mcp)$mean
  params_lm = c(as.numeric(coef(fit_lm)), summary(fit_lm)$sigma)
  testthat::expect_lt(max(abs(params_mcp - params_lm)), 0.1)

  # Log-likelihood
  fit_mcp = add_loglik(fit_mcp)
  loglik_mcp = mean(rowSums(fit_mcp$loglik))
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
    adapt = 500,
    iter = 2000,
    seed = 42,
    diagnostics = FALSE,
    quiet = TRUE
  )

  fit_lm_cp = lm(y ~ x + I(pmax(0, x - 10)), data = df_cp)
  coef_lm = coef(fit_lm_cp)

  # mcp parameter order for this model: cp_1, Intercept_1, x_1, x_2, sigma_1
  mcp_fixef = fixef(fit_mcp_cp)
  estimates_mcp = c(
    mcp_fixef$mean[mcp_fixef$name == "Intercept_1"],
    mcp_fixef$mean[mcp_fixef$name == "x_1"],
    mcp_fixef$mean[mcp_fixef$name == "x_2"],
    mcp_fixef$mean[mcp_fixef$name == "sigma_1"]
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
