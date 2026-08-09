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
  expect_false(grepl("uncentered", fit$jags_code, fixed = TRUE))
  expect_match(fit$jags_code, "Order realized group-level change points")
  expect_match(fit$jags_code, "cp_order_\\[id_, 1\\] ~ dbern\\(step")
})


test_that("simulation rejects unordered group-level change points", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    sample = FALSE
  )

  expect_error(
    fit$simulate(
      fit,
      cp_group_data,
      cp_1 = 3,
      cp_2 = 6,
      cp_1_id = 4,
      cp_2_id = -4,
      Intercept_1 = 0,
      Intercept_2 = 0,
      Intercept_3 = 0,
      sigma_1 = 1,
      .type = "fitted"
    ),
    "must remain ordered"
  )
})
