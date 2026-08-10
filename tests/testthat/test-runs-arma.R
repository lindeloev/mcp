#############
# TEST ARMA #
#############
# We can assume that it will fail for the same mis-specifications on the formula
# ar(order, [formula]), since the formula runs through the exact same code as
# sigma and mu.
bad_arma = list(
  list(y ~ ar(0)),  # currently not implemented
  list(y ~ ar(-1)),  # must be positive
  list(y ~ ar(1.5)),  # Cannot be decimal
  list(y ~ ar(1) + ar(2)),  # Only one per segment
  list(y ~ ar("1")),  # Should not take strings
  list(y ~ ar(1 + x)),  # must have order
  list(y ~ ar(x)),  # must have order
  list(y ~ ma(0)),
  list(y ~ ma(-1)),
  list(y ~ ma(1.5)),
  list(y ~ ma(1) + ma(2))
)

test_bad(bad_arma)


good_arma = list(
  list(y ~ ar(1)),  # simple
  list(y ~ ma(1)),
  list(y ~ ar(1) + ma(1)),
  list(y ~ ar(3)),  # higher order
  list(y ~ ar(1, 1 + x + I(x^2) + exp(x))),  # complicated regression
  list(y ~ ar(1),
       ~ ar(2, 0 + x)),  # change in ar
  list(y ~ 1,
       ~ 0 + ar(2)),  # onset of AR
  list(y ~ 1,
       1 + (1|id) ~ 1 + I(x^2) + ar(2, 1 + x)),  # varying change point
  list(y ~ ar(1) + sigma(1 + x),
       ~ ar(2, 1 + I(x^2)) + sigma(1)),  # With sigma
  list(y ~ ar(1),
       ~ ar(2, 1)),
  list(y | weights(weights_ok) ~ 1 + ar(1),  # With weights
       ~ 0 + ar(2, 1 + x)),
  list(y ~ 0 + ar(1),
       ~ 0 + ar(2))  # mu is ~0 across all segments; only ar() gives structure
)

test_good(good_arma)


test_that("series resets generated AR/MA lags", {
  data = data.frame(
    id = rep(c("a", "b"), each = 2),
    x = c(1, 2, 1, 3),
    y = 1:4
  )
  fit = suppressMessages(mcp(
    list(y ~ 1 + ar(2) + ma(1)), data,
    par_x = "x", series = "id", sample = FALSE, quiet = TRUE
  ))

  expect_match(
    fit$jags_code,
    "equals(series_id_[i_], series_id_[i_ - 2]) * ar2_[i_]",
    fixed = TRUE
  )
  expect_match(
    fit$jags_code,
    "equals(series_id_[i_], series_id_[i_ - 1]) * ma1_[i_]",
    fixed = TRUE
  )
  expect_false(grepl("series_id_", mcp(
    list(y ~ 1 + ar(1)), data, par_x = "x", sample = FALSE, quiet = TRUE
  )$jags_code, fixed = TRUE))
})


test_that("series input is contiguous and survives interpolation", {
  data = data.frame(id = c("a", "b", "a"), x = 1:3, y = 1:3)

  expect_error(
    mcp(list(y ~ ar(1)), data, par_x = "x", series = "id", sample = FALSE),
    "Rows belonging to each `series` must be contiguous.",
    fixed = TRUE
  )
  ordered = data[order(data$id), ]
  fit = mcp(
    list(y ~ ar(1)), ordered,
    par_x = "x", series = "id", sample = FALSE, quiet = TRUE
  )
  interpolated = interpolate_newdata(fit, by = "id", x_values = 11:13)
  expect_identical(interpolated$id, ordered$id)
  expect_equal(interpolated$x, 11:13)
})


test_that("AR/MA model checks warn about conditional validation", {
  data = data.frame(x = 1:4, y = 1:4)
  fit = mcp(list(y ~ 1 + ar(1)), data = data, par_x = "x", sample = FALSE, quiet = TRUE)

  expect_warning(
    try(pp_check(fit, ndraws = 1), silent = TRUE),
    "one-step-ahead conditional predictions",
    fixed = TRUE
  )
  expect_warning(
    try(loo(fit), silent = TRUE),
    "Observationwise PSIS-LOO/WAIC is problematic",
    fixed = TRUE
  )
  expect_warning(
    try(waic(fit), silent = TRUE),
    "Observationwise PSIS-LOO/WAIC is problematic",
    fixed = TRUE
  )
  expect_no_warning(warn_arma_check(fit, arma = FALSE, "ppc"))
})


test_that("ar() and ma() require a family GARMA implementation", {
  data = data.frame(x = 1:4, y = 1:4)
  family = gaussian(link = "log")

  expect_null(mcpfamily(family)$garma)
  expect_s3_class(
    mcp(list(y ~ 1 + x), data, family = family, par_x = "x", sample = FALSE),
    "mcpfit"
  )

  expect_error(
    mcp(
      list(y ~ 1 + ar(1) + ma(1)),
      data,
      family = family,
      par_x = "x",
      sample = FALSE
    ),
    "family = gaussian(link = \"log\") does not define the GARMA behavior required by ar() or ma().",
    fixed = TRUE
  )
})
