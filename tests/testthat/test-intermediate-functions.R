##########################################################
# Direct tests of interpolate_newdata() and fit$simulate() #
##########################################################
# These intermediate functions are otherwise only exercised indirectly
# (via plot()/fitted()/predict() and test_arma_simulation() in helper-runs.R).

test_that("interpolate_newdata() combines categorical and continuous predictors", {
  data = data.frame(
    x = c(1, 2, 3, 1, 2, 3),
    z = c(2, 4, 6, 8, 10, 12),
    group = rep(c("A", "B"), each = 3),
    y = 0
  )
  fit = mcp(
    list(y ~ 1 + z + group),
    data,
    par_x = "x",
    sample = FALSE
  )

  # The has_continuous check probes the full (un-grouped) data, which has
  # duplicate x values across groups; suppress the resulting approx() notice.
  newdata = suppressWarnings(interpolate_newdata(fit, x_values = c(1, 2, 3)))

  # Categorical levels are combined factorially with x, continuous
  # predictors are interpolated separately within each group.
  expect_equal(nrow(newdata), 6)
  expect_setequal(unique(newdata$group), c("A", "B"))
  expect_equal(
    newdata$z[newdata$group == "A" & newdata$x == 2],
    4
  )
  expect_equal(
    newdata$z[newdata$group == "B" & newdata$x == 2],
    10
  )
})

test_that("interpolate_newdata() carries the observed response for AR/MA models", {
  data = data.frame(x = 1:5, y = c(1, 2, 3, 4, 5))
  fit = mcp(
    list(y ~ ar(1)),
    data,
    par_x = "x",
    sample = FALSE
  )

  newdata = interpolate_newdata(fit)

  expect_equal(nrow(newdata), nrow(data))
  expect_equal(newdata[[fit$pars$y]], data$y)
})

test_that("fit$simulate() reproduces the deterministic linear predictor at .type = 'fitted'", {
  data = data.frame(x = 1:5, y = 0)
  fit = mcp(
    list(y ~ 1 + x),
    data,
    par_x = "x",
    sample = FALSE
  )

  simulated = fit$simulate(
    fit,
    data,
    Intercept_1 = 2,
    x_1 = 3,
    sigma_1 = 1,
    .type = "fitted"
  )

  expect_equal(as.numeric(simulated), 2 + 3 * data$x)
})

test_that("fit$simulate() errors with a helpful message on the pre-v0.4 calling convention", {
  data = data.frame(x = 1:5, y = 0)
  fit = mcp(
    list(y ~ 1 + x),
    data,
    par_x = "x",
    sample = FALSE
  )

  expect_error(
    fit$simulate(1:5, data),
    "breaking changes"
  )
})
