cp_group_data = data.frame(
  x = rep(1:10, 3),
  y = seq_len(30),
  id = factor(rep(letters[1:3], each = 10)),
  site = factor(rep(LETTERS[1:3], 10))
)


test_that("group-level change points use one grouping factor", {
  expect_error(
    mcp(
      list(y ~ 1, 1 + (1 | id) ~ 1, 1 + (1 | site) ~ 1),
      cp_group_data,
      par_x = "x",
      sample = FALSE
    ),
    "must use the same grouping factor"
  )
})


test_that("group-level change-point priors and JAGS code enforce ordering", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    sample = FALSE
  )

  expect_match(fit$prior$cp_2_id, "cp_1_id\\[id_\\]")
  expect_match(fit$jags_code, "cp_1_id = cp_1_id_uncentered - mean")
  expect_match(fit$jags_code, "Order realized group-level change points")
  expect_match(fit$jags_code, "cp_order_\\[id_, 1\\] ~ dbern\\(step")
})


test_that("simulation generates centered change-point deviations from SDs", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    sample = FALSE
  )

  expect_true("cp_1_sd" %in% names(formals(fit$simulate)))
  expect_false("cp_1_id" %in% names(formals(fit$simulate)))
  set.seed(42)
  simulated = fit$simulate(
    fit,
    cp_group_data,
    cp_1 = 5,
    cp_1_sd = 1,
    Intercept_1 = 0,
    Intercept_2 = 0,
    sigma_1 = 1,
    .type = "fitted"
  )
  deviations = attr(simulated, "simulated")$cp_1_id
  by_group = vapply(split(deviations, cp_group_data$id), unique, numeric(1))
  expect_equal(mean(by_group), 0)
})


test_that("simulation rejects unordered group-level change points", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    sample = FALSE
  )

  set.seed(1)
  expect_error(
    fit$simulate(
      fit,
      cp_group_data,
      cp_1 = 3,
      cp_2 = 6,
      cp_1_sd = 10,
      cp_2_sd = 10,
      Intercept_1 = 0,
      Intercept_2 = 0,
      Intercept_3 = 0,
      sigma_1 = 1,
      .type = "fitted"
    ),
    "must remain ordered"
  )
})
