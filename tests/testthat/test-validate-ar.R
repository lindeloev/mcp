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
  params_mcp = summary_mcp$mean[match(c("ar1_1", "ar2_1", "Intercept_1"), summary_mcp$name)]
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
