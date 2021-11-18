
################
# AR INFERENCE #
################
# AR
model = list(
  y ~ 1 + ar(2) + group
)

N = 300
df = data.frame(
  x = 1:N,
  y = arima.sim(list(ar = c(0.7, -0.4)), N),
  group = rep(c("A", "B"), 150)
)

fit_arima = arima(df$y, order = c(2,0,0))
fit_mcp = mcp(model, df, par_x = "x", adapt = 100, iter = 1000, chains = 2)

# Parameter estimates
testthat::expect_equal(as.numeric(fit_arima$coef), fixef(fit_mcp)$mean[c(1,2, 4)], tolerance = 0.03)

# Log-likelihood
fit_mcp = add_loglik(fit_mcp)
loglik_mcp = mean(rowSums(fit_mcp$loglik))
loglik_arima = as.numeric(logLik(fit_arima))
expect_equal(loglik_arima, loglik_mcp, tolerance = 0.01)


#################
# AR SIMULATION #
#################

newdata = dplyr::select(df, -y) %>% tidyr::expand_grid(rep = 1:100)
y_simulated = fit_mcp$simulate(
  fit_mcp, newdata,
  Intercept_1 = 9, sigma_1 = 2, ar1_1 = 0.7, ar2_1 = -0.3, groupB_1 = 0
)
y_arima = arima(y_simulated, order = c(2, 0, 0))

expect_equal(as.numeric(y_arima$coef), c(0.7, -0.3, 9), tolerance = 0.01)
