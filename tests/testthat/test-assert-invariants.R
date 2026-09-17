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
  fit = mcp(list(y ~ 1 + sqrt(z)), data = df_valid, par_x = "x", sample = "prior", quiet = TRUE)
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


test_that("mcp rejects non-syntactic data column names", {
  data = data.frame(`x value` = 1:3, y = 0, check.names = FALSE)

  expect_error(
    mcp(list(y ~ 1), data, par_x = "x value", sample = FALSE),
    "`data` has non-syntactic column name(s): `x value`",
    fixed = TRUE
  )
})


test_that("mcp rejects data names that collide with generated JAGS nodes", {
  data = data.frame(mu_ = 1:5, x = 1:5)
  expect_error(
    mcp(list(mu_ ~ 1), data, par_x = "x", sample = FALSE),
    "Data column name(s) collide with mcp's generated JAGS namespace: 'mu_'",
    fixed = TRUE
  )
})


test_that("mcp rejects generated parameter names that collide with change points", {
  data = data.frame(x = 1:5, cp = c(0, 1, 0, 1, 0), y = 1:5)
  expect_error(
    mcp(list(y ~ cp, ~ 1), data, par_x = "x", sample = FALSE),
    "Generated parameter name(s) collide in the JAGS namespace: 'cp_1'",
    fixed = TRUE
  )
})


test_that("parameter-name collisions give a useful error", {
  data = data.frame(
    y = 1:6,
    x = 1:6,
    a = c(0, 0, 0, 1, 1, 1),
    b = c(0, 1, 2, 0, 1, 2),
    ab = c(0, 1, 0, 1, 0, 1)
  )

  expect_error(
    mcp(list(y ~ a:b + ab), data, par_x = "x", sample = FALSE),
    "`ab_1`: `ab` (mu, segment 1) and `a:b` (mu, segment 1)",
    fixed = TRUE
  )
})


test_that("zero and negative weights are rejected with informative error", {
  data_zero = data.frame(x = 1:5, y = 1:5, w = c(1, 1, 0, 1, 1))
  expect_error(
    mcp(list(y | weights(w) ~ 1 + x), data = data_zero, sample = FALSE),
    "All weights must be numeric and greater than zero.",
    fixed = TRUE
  )

  data_neg = data.frame(x = 1:5, y = 1:5, w = c(1, 1, -0.5, 1, 1))
  expect_error(
    mcp(list(y | weights(w) ~ 1 + x), data = data_neg, sample = FALSE),
    "All weights must be numeric and greater than zero.",
    fixed = TRUE
  )
})


test_that("probs and quantiles must be strictly between 0 and 1", {
  expect_error(fitted(demo_fit, probs = 0), "strictly between 0 and 1")
  expect_error(fitted(demo_fit, probs = 1), "strictly between 0 and 1")
  expect_error(fitted(demo_fit, probs = -0.1), "strictly between 0 and 1")
  expect_error(fitted(demo_fit, probs = 1.1), "strictly between 0 and 1")
  expect_error(predict(demo_fit, probs = 0), "strictly between 0 and 1")
  expect_error(predict(demo_fit, probs = 1), "strictly between 0 and 1")
  expect_error(residuals(demo_fit, probs = 0), "strictly between 0 and 1")
  expect_error(residuals(demo_fit, probs = 1), "strictly between 0 and 1")

  expect_error(plot(demo_fit, q_fit = 0), "strictly between 0 and 1")
  expect_error(plot(demo_fit, q_fit = 1), "strictly between 0 and 1")
  expect_error(plot(demo_fit, q_predict = 0), "strictly between 0 and 1")
  expect_error(plot(demo_fit, q_predict = 1), "strictly between 0 and 1")
})





