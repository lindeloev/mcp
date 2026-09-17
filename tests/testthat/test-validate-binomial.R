if (Sys.getenv("MCP_TEST_LEVEL") != "release") {
  testthat::skip("Time-consuming validation tests against reference implementations are only run when MCP_TEST_LEVEL='release'.")
}

# Fit a simple binomial model on simulated data
set.seed(42)
model = list(
  y | trials(N) ~ 1 + x + group
)

df = tibble::tibble(
  N = round(runif(200, 1, 10)),
  x = seq(0, 20, length.out = 200),
  group = rep(c("A", "B"), 100),
  y = rbinom(200, N, ilogit(2 - 0.1 * x + ifelse(group == "B", -1, 0)))
)

fit_mcp = mcp(model, df, family = binomial(), warmup = 100, iter = 1000, seed = 42, diagnostics = FALSE, quiet = TRUE)

# Tests
test_that("Binomial inference against glm()", {
  fit_glm = glm(cbind(y, N - y) ~ x + group, data = df, family = binomial())

  # Parameter estimates. Both are in (Intercept, x, groupB) order.
  params_mcp = fixef(fit_mcp)$mean
  params_glm = as.numeric(fit_glm$coefficients)
  testthat::expect_lt(max(abs(params_mcp - params_glm)), 0.05)

  # Log-likelihood
  loglik_mcp = mean(rowSums(log_lik(fit_mcp)))
  loglik_glm = as.numeric(logLik(fit_glm))
  expect_equal(loglik_glm, loglik_mcp, tolerance = 0.01)
})


test_that("Binomial simulation against glm()", {
  newdata = df %>% dplyr::select(-"y") %>% tidyr::expand_grid(rep = c(1:100))
  newdata$y = fit_mcp$simulate(fit_mcp, newdata, Intercept_1 = 2, x_1 = -0.1, groupB_1 = -1)

  fit_glm_sim = glm(cbind(y, N - y) ~ x + group, data = newdata, family = binomial())

  testthat::expect_equal(as.numeric(fit_glm_sim$coefficients), c(2, -0.1, -1), tolerance = 0.02)
})


test_that("Binomial likelihood weights against glm(weights = w)", {
  set.seed(42)
  N_obs = 200
  x = seq(0, 10, length.out = N_obs)
  w = runif(N_obs, 0.5, 2.5)
  trials = round(runif(N_obs, 5, 20))
  eta = 0.5 - 0.2 * x
  p = 1 / (1 + exp(-eta))
  y = rbinom(N_obs, size = trials, prob = p)
  df_w = data.frame(x = x, y = y, trials = trials, w = w)

  fit_glm_w = glm(cbind(y, trials - y) ~ x, data = df_w, family = binomial(), weights = w)
  coef_glm_w = coef(fit_glm_w)

  fit_mcp_w = mcp(
    list(y | trials(trials) + weights(w) ~ 1 + x),
    data = df_w,
    family = binomial(),
    warmup = 500,
    iter = 2000,
    seed = 42,
    diagnostics = FALSE,
    quiet = TRUE
  )

  fix_mcp_w = fixef(fit_mcp_w)
  expect_equal(fix_mcp_w$mean[fix_mcp_w$variable == "Intercept_1"], unname(coef_glm_w[1]), tolerance = 0.05)
  expect_equal(fix_mcp_w$mean[fix_mcp_w$variable == "x_1"], unname(coef_glm_w[2]), tolerance = 0.02)

  # Log-likelihood
  loglik_mcp = mean(rowSums(log_lik(fit_mcp_w)))
  loglik_glm = as.numeric(logLik(fit_glm_w))
  expect_equal(loglik_glm - loglik_mcp, 1.0, tolerance = 0.2)
})
