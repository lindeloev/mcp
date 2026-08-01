if (Sys.getenv("MCP_TEST_LEVEL") != "release") {
  testthat::skip("Time-consuming validation tests against reference implementations are only run when MCP_TEST_LEVEL='release'.")
}

# Fit a simple Poisson model on simulated data
set.seed(42)
model = list(
  y ~ 1 + x + group
)

df = tibble::tibble(
  x = seq(0, 20, length.out = 200),
  group = rep(c("A", "B"), 100),
  y = rpois(200, lambda = exp(1 + 0.15 * x + ifelse(group == "B", -0.5, 0)))
)

fit_mcp = quiet_mcp(model, df, family = poisson(), adapt = 500, iter = 2000, seed = 42)

# Tests
test_that("Poisson inference against glm()", {
  fit_glm = glm(y ~ x + group, data = df, family = poisson())

  # Parameter estimates. Both are in (Intercept, x, groupB) order.
  params_mcp = fixef(fit_mcp)$mean
  params_glm = as.numeric(coef(fit_glm))
  testthat::expect_lt(max(abs(params_mcp - params_glm)), 0.05)

  # Log-likelihood
  fit_mcp = add_loglik(fit_mcp)
  loglik_mcp = mean(rowSums(fit_mcp$loglik))
  loglik_glm = as.numeric(logLik(fit_glm))
  expect_equal(loglik_glm, loglik_mcp, tolerance = 0.01)
})


test_that("Poisson fixed change-point inference against glm()", {
  set.seed(42)
  df_cp = tibble::tibble(
    x = seq(0, 10, length.out = 200),
    y = rpois(200, lambda = exp(1 + 0.1 * x + 0.2 * pmax(0, x - 5)))
  )

  fit_mcp_cp = quiet_mcp(
    list(y ~ 1 + x, ~ 0 + x),
    df_cp,
    family = poisson(),
    prior = list(cp_1 = 5),
    adapt = 500,
    iter = 2000,
    seed = 42
  )

  fit_glm_cp = glm(y ~ x + I(pmax(0, x - 5)), data = df_cp, family = poisson())
  coef_glm = coef(fit_glm_cp)

  mcp_fixef = fixef(fit_mcp_cp)
  estimates_mcp = c(
    mcp_fixef$mean[mcp_fixef$name == "Intercept_1"],
    mcp_fixef$mean[mcp_fixef$name == "x_1"],
    mcp_fixef$mean[mcp_fixef$name == "x_2"]
  )

  estimates_glm = c(
    unname(coef_glm[1]),
    unname(coef_glm[2]),
    unname(coef_glm[2] + coef_glm[3])
  )

  testthat::expect_lt(max(abs(estimates_mcp - estimates_glm)), 0.05)
})


test_that("Poisson simulation against glm()", {
  newdata = df %>% dplyr::select(-"y") %>% tidyr::expand_grid(rep = 1:100)
  newdata$y = fit_mcp$simulate(fit_mcp, newdata, Intercept_1 = 1, x_1 = 0.15, groupB_1 = -0.5)

  fit_glm_sim = glm(y ~ x + group, data = newdata, family = poisson())

  testthat::expect_equal(as.numeric(coef(fit_glm_sim)), c(1, 0.15, -0.5), tolerance = 0.02)
})
