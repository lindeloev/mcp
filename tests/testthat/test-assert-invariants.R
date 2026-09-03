test_that("model-defining data cannot be missing", {
  binomial_data = data.frame(x = 1:4, y = c(0, 1, 1, 0), trials = c(2, NA, 2, 2))
  expect_error(
    mcp(list(y | trials(trials) ~ 1), binomial_data, family = binomial(), par_x = "x", sample = FALSE),
    "trials"
  )

  group_data = data.frame(x = 1:4, y = 0, id = c("a", NA, "b", "b"))
  expect_error(
    mcp(list(y ~ 1, 1 + (1 | id) ~ 1), group_data, par_x = "x", sample = FALSE),
    "id"
  )

  x_data = data.frame(x = c(1, NA, 3, 4), y = 0)
  expect_error(
    mcp(list(y ~ 1), x_data, par_x = "x", sample = FALSE),
    "x"
  )
})


test_that("simulation requires ordered population change points", {
  data = data.frame(x = 1:5, y = 0)
  fit = mcp(
    list(y ~ 1, 1 ~ 1, 1 ~ 1), data,
    par_x = "x", sample = FALSE
  )

  expect_error(
    fit$simulate(
      fit, data,
      cp_1 = 4, cp_2 = 2,
      Intercept_1 = 0, Intercept_2 = 0, Intercept_3 = 0, sigma_1 = 1,
      .type = "fitted"
    ),
    "Population-level change points must remain strictly ordered.",
    fixed = TRUE
  )
})


test_that("simulation validates required response auxiliaries", {
  data = data.frame(x = 1:4, y = c(0, 1, 1, 0), trials = 2)
  fit = mcp(
    list(y | trials(trials) ~ 1), data,
    family = binomial(), par_x = "x", sample = FALSE
  )

  newdata = data.frame(x = 1:4, trials = c(2, NA, 2, 2))
  expect_error(
    fit$simulate(fit, newdata, Intercept_1 = 0, .type = "fitted"),
    "trials"
  )
})


test_that("data columns colliding with reserved output namespace are rejected early", {
  reserved = c("sd", "fitted", "predict", "residuals", "loglik", ".draw", ".chain", ".iteration", "data_row")
  for (name in reserved) {
    bad_data = data.frame(x = 1:4, y = 1:4)
    bad_data[[name]] = c(1, 2, 1, 2)
    expect_error(
      mcp(list(stats::as.formula(paste("y ~ 1 +", name)), ~ 1), data = bad_data, par_x = "x", sample = FALSE),
      "reserved output namespace"
    )
  }
})


test_that("newdata columns colliding with reserved output namespace are rejected early", {
  data = data.frame(x = 1:4, y = 1:4)
  fit = mcp(list(y ~ 1, ~ 1), data = data, par_x = "x", sample = FALSE)

  expect_error(
    fitted(fit, newdata = data.frame(x = 1:2, sd = 1:2)),
    "reserved output namespace"
  )
  expect_error(
    predict(fit, newdata = data.frame(x = 1:2, data_row = 1:2)),
    "reserved output namespace"
  )
})

