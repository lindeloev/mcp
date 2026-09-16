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


test_that("data columns colliding with reserved offset namespace are rejected early", {
  for (name in c("offset", "offset_1", "offset_2")) {
    bad_data = data.frame(x = 1:4, y = 1:4)
    bad_data[[name]] = 1:4
    expect_error(
      mcp(list(y ~ 1, ~ 1), data = bad_data, par_x = "x", sample = FALSE),
      "reserved offset namespace"
    )
  }
})


test_that("response column named offset is rejected early", {
  bad_data = data.frame(x = 1:4, offset = 1:4, exposure = 10)
  expect_error(
    mcp(list(offset ~ 1 + offset(log(exposure))), data = bad_data, par_x = "x", sample = FALSE),
    "reserved offset namespace"
  )
})


test_that("transformed predictors and offsets producing NA or non-finite values are rejected", {
  df = data.frame(x = 1:3, y = 1:3, z = c(-1, 4, -1))
  expect_error(
    suppressWarnings(mcp(list(y ~ 0 + sqrt(z)), data = df, par_x = "x", sample = FALSE)),
    "Predictor transformation resulted in NA or non-finite values: sqrt(z)",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(mcp(list(y ~ 1 + offset(sqrt(z))), data = df, par_x = "x", sample = FALSE)),
    "Predictor transformation resulted in NA or non-finite values: offset(sqrt(z))",
    fixed = TRUE
  )

  df_valid = data.frame(x = 1:3, y = 1:3, z = c(4, 9, 16))
  fit = mcp(list(y ~ 1 + sqrt(z)), data = df_valid, par_x = "x", sample = "prior")
  newdata = data.frame(x = 1:2, z = c(-1, 4))
  expect_error(
    suppressWarnings(fitted(fit, newdata = newdata, prior = TRUE)),
    "Predictor transformation resulted in NA or non-finite values: sqrt(z)",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(predict(fit, newdata = newdata, prior = TRUE)),
    "Predictor transformation resulted in NA or non-finite values: sqrt(z)",
    fixed = TRUE
  )
})



