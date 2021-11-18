#############
# INFERENCE #
#############
model = list(
  y | trials(N) ~ 1 + x + group
)

df = tibble::tibble(
  N = round(runif(200, 1, 10)),
  x = seq(0, 20, length.out = 200),
  group = rep(c("A", "B"), 100),
  y = rbinom(200, N, ilogit(2 - 0.1 * x + ifelse(group == "B", -1, 0)))
)

fit_mcp = mcp(model, df, family = binomial(), adapt = 100, iter = 1000)
fit_glm = glm(cbind(y, N - y) ~ x + group, data = df, family = binomial())

# Parameter estimates
params_mcp = fixef(fit_mcp)$mean
params_glm = as.numeric(fit_glm$coefficients[c(3, 1, 2)])
testthat::expect_equal(params_mcp, params_glm, tolerance = 0.02)

# Log-likelihood
fit_mcp = add_loglik(fit_mcp)
loglik_mcp = mean(rowSums(fit_mcp$loglik))
loglik_glm = as.numeric(logLik(fit_glm))
expect_equal(loglik_glm, loglik_mcp, tolerance = 0.01)


##############
# SIMULATION #
##############
newdata = dplyr::select(df, -y) %>% tidyr::expand_grid(rep = c(1:100))
newdata$y = fit_mcp$simulate(fit_mcp, newdata, Intercept_1 = 2, x_1 = -0.1, groupB_1 = -1)

fit_glm_sim = glm(cbind(y, N - y) ~ x + group, data = newdata, family = binomial())


testthat::expect_equal(as.numeric(fit_glm_sim$coefficients), c(2, -0.1, -1), tolerance = 0.02)
