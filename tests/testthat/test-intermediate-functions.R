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

  newdata = interpolate_newdata(fit, x_values = c(1, 2, 3))

  # Categorical levels are combined factorially with x, continuous
  # predictors are held at their observed means.
  expect_equal(nrow(newdata), 6)
  expect_setequal(unique(newdata$group), c("A", "B"))
  expect_equal(unique(newdata$z), mean(data$z))

  custom = interpolate_newdata(fit, x_values = 1:3, at = list(z = 4))
  expect_equal(unique(custom$z), 4)
  expect_error(interpolate_newdata(fit, at = list(x = 2)), "Invalid: 'x'.", fixed = TRUE)
  expect_error(interpolate_newdata(fit, at = list(z = 1:2)), "must be a single number")

  population = mcp_pars(fit, scope = "population")$name
  draws = matrix(0, 1, length(population), dimnames = list(NULL, population))
  draws[, "z_1"] = 2
  draws[, "sigma_1"] = 1
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))
  plotted = plot(fit, at = list(z = 4), lines = 1)
  expect_match(plotted$labels$caption, "z = 4", fixed = TRUE)
  expect_equal(unique(ggplot2::ggplot_build(plotted)$data[[2]]$y), 8)
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
  expect_equal(newdata[[mcp_columns(fit)$response]], data$y)
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
