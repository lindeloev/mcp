if (Sys.getenv("MCP_TEST_LEVEL") != "release") {
  testthat::skip("Time-consuming validation tests against reference implementations are only run when MCP_TEST_LEVEL='release'.")
}

# Fit a simple AR model on simulated data
model = list(
  y ~ 1 + ar(2) + group
)

N = 300
df = data.frame(
  x = 1:N,
  y = 2 + arima.sim(list(ar = c(0.7, -0.4)), N),
  group = rep(c("A", "B"), 150)
)

fit_mcp = mcp(model, df, par_x = "x", warmup = 100, iter = 1000, chains = 2, diagnostics = FALSE, quiet = TRUE)


# Test stuff
test_that("AR inference against arima()", {
  fit_arima = arima(df$y, order = c(2,0,0))

  # Parameter estimates. arima() returns c(ar1, ar2, intercept).
  params_arima = as.numeric(fit_arima$coef)
  capture.output({ summary_mcp = summary(fit_mcp) })
  params_mcp = summary_mcp$mean[match(c("ar1_1", "ar2_1", "Intercept_1"), summary_mcp$variable)]
  testthat::expect_equal(params_arima, params_mcp, tolerance = 0.03)

  # Log-likelihood
  loglik_mcp = mean(rowSums(log_lik(fit_mcp)))
  loglik_arima = as.numeric(logLik(fit_arima))
  expect_equal(loglik_arima, loglik_mcp, tolerance = 0.05)
})


test_that("AR simulation against arima.sim()", {
  newdata = tidyr::expand_grid(df, rep = 1:100)  # should work with y in data too
  expect_message(
    {
      y_simulated = fit_mcp$simulate(
        fit_mcp, newdata,
        Intercept_1 = 9, sigma_1 = 2, ar1_1 = 0.7, ar2_1 = -0.3, groupB_1 = 0
      )
    },
    "Generating residuals for AR\\(N\\) model"
  )
  y_arima = arima(y_simulated, order = c(2, 0, 0))

  expect_equal(as.numeric(y_arima$coef), c(0.7, -0.3, 9), tolerance = 0.01)
})


test_that("MA(1) inference against arima()", {
  set.seed(42)
  N_ma = 400
  eps = rnorm(N_ma, sd = 1.0)
  y_ma = numeric(N_ma)
  theta = 0.5
  mu = 3.0
  y_ma[1] = mu + eps[1]
  for (t in 2:N_ma) {
    y_ma[t] = mu + eps[t] + theta * eps[t - 1]
  }
  df_ma = data.frame(x = 1:N_ma, y = y_ma)

  fit_arima_ma = arima(df_ma$y, order = c(0, 0, 1))

  fit_mcp_ma = mcp(
    list(y ~ 1 + ma(1)),
    data = df_ma,
    par_x = "x",
    warmup = 300,
    iter = 1500,
    chains = 2,
    seed = 42,
    diagnostics = FALSE,
    quiet = TRUE
  )

  capture.output({ summary_mcp_ma = summary(fit_mcp_ma) })
  mcp_ma1 = summary_mcp_ma$mean[summary_mcp_ma$variable == "ma1_1"]
  mcp_int = summary_mcp_ma$mean[summary_mcp_ma$variable == "Intercept_1"]
  mcp_sigma = summary_mcp_ma$mean[summary_mcp_ma$variable == "sigma_1"]

  expect_equal(mcp_ma1, unname(fit_arima_ma$coef["ma1"]), tolerance = 0.05)
  expect_equal(mcp_int, unname(fit_arima_ma$coef["intercept"]), tolerance = 0.05)
  expect_equal(mcp_sigma, sqrt(fit_arima_ma$sigma2), tolerance = 0.05)
})


test_that("ARMA(1, 1) inference against arima()", {
  set.seed(42)
  N_arma = 500
  y_arma = 2.5 + as.numeric(arima.sim(list(ar = 0.5, ma = 0.4), n = N_arma, sd = 0.9))
  df_arma = data.frame(x = 1:N_arma, y = y_arma)

  fit_arima_arma = arima(df_arma$y, order = c(1, 0, 1))

  fit_mcp_arma = mcp(
    list(y ~ 1 + ar(1) + ma(1)),
    data = df_arma,
    par_x = "x",
    warmup = 300,
    iter = 1500,
    chains = 2,
    seed = 42,
    diagnostics = FALSE,
    quiet = TRUE
  )

  capture.output({ summary_mcp_arma = summary(fit_mcp_arma) })
  mcp_ar1 = summary_mcp_arma$mean[summary_mcp_arma$variable == "ar1_1"]
  mcp_ma1 = summary_mcp_arma$mean[summary_mcp_arma$variable == "ma1_1"]
  mcp_int = summary_mcp_arma$mean[summary_mcp_arma$variable == "Intercept_1"]
  mcp_sigma = summary_mcp_arma$mean[summary_mcp_arma$variable == "sigma_1"]

  expect_equal(mcp_ar1, unname(fit_arima_arma$coef["ar1"]), tolerance = 0.05)
  expect_equal(mcp_ma1, unname(fit_arima_arma$coef["ma1"]), tolerance = 0.05)
  expect_equal(mcp_int, unname(fit_arima_arma$coef["intercept"]), tolerance = 0.05)
  expect_equal(mcp_sigma, sqrt(fit_arima_arma$sigma2), tolerance = 0.05)
})


test_that("MA simulation against arima()", {
  df_sim = data.frame(x = 1:500, y = 0)
  fit_sim = mcp(list(y ~ 1 + ma(1)), data = df_sim, par_x = "x", sample = FALSE)
  suppressMessages({
    y_simulated = fit_sim$simulate(fit_sim, df_sim, Intercept_1 = 3, sigma_1 = 1, ma1_1 = 0.5)
  })
  fit_arima_sim = arima(y_simulated, order = c(0, 0, 1))

  expect_equal(unname(fit_arima_sim$coef["ma1"]), 0.5, tolerance = 0.05)
  expect_equal(unname(fit_arima_sim$coef["intercept"]), 3.0, tolerance = 0.05)
})
