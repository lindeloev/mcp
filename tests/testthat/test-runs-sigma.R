#################
# TEST GAUSSIAN STANDARD DEVIATION #
#################
bad_sigma = list(
  list(y ~ 1 + sigma(q))  # variable does not exist
)

test_bad(bad_sigma)


test_that("fixed residual SDs must meet the family floor", {
  data = data.frame(x = 1:4, y = 0)

  expect_error(
    mcp(
      list(y ~ 1), data,
      par_x = "x", prior = list(sigma_1 = 0), sample = FALSE
    ),
    "Fixed residual standard deviation parameter(s) must be at least 0.001: sigma_1.",
    fixed = TRUE
  )
  expect_error(
    mcp(
      list(y ~ 1), data,
      par_x = "x", prior = list(sigma_1 = 0.0001), sample = FALSE
    ),
    "Fixed residual standard deviation parameter(s) must be at least 0.001: sigma_1.",
    fixed = TRUE
  )
  expect_error(
    mcp(
      list(y ~ 1), data,
      par_x = "x", prior = list(sigma_1 = -1), sample = FALSE
    ),
    "Fixed residual standard deviation parameter(s) must be at least 0.001: sigma_1.",
    fixed = TRUE
  )
  expect_silent(
    mcp(
      list(y ~ 1), data,
      par_x = "x", prior = list(sigma_1 = 0.001), sample = FALSE
    )
  )
  expect_silent(
    mcp(
      list(y ~ 1 + sigma(1)), data,
      par_x = "x", prior = list(sigma_1 = 0), sample = FALSE
    )
  )
})


test_that("distribution priors for fixed residual SDs use the likelihood lower bound", {
  data = data.frame(x = 1:4, y = 0)
  fit = mcp(
    list(y ~ 1), data,
    par_x = "x", prior = list(sigma_1 = "dnorm(0, 1)"), sample = FALSE
  )

  expect_equal(fit$prior$sigma_1, "dnorm(0, 1) T(0.001, )")
  expect_match(fit$jags_code, "sigma_1 ~ dnorm(0, 1/(1)^2) T(0.001,)", fixed = TRUE)

  # Already bounded uniform prior is preserved without T()
  fit_unif_ok = mcp(
    list(y ~ 1), data,
    par_x = "x", prior = list(sigma_1 = "dunif(1, 2)"), sample = FALSE
  )
  expect_equal(fit_unif_ok$prior$sigma_1, "dunif(1, 2)")
  expect_match(fit_unif_ok$jags_code, "sigma_1 ~ dunif(1, 2)", fixed = TRUE)
  expect_false(grepl("sigma_1 ~ dunif.*T\\(", fit_unif_ok$jags_code))

  # Uniform prior reaching below floor has lower bound adjusted
  fit_unif_floor = mcp(
    list(y ~ 1), data,
    par_x = "x", prior = list(sigma_1 = "dunif(0, 2)"), sample = FALSE
  )
  expect_equal(fit_unif_floor$prior$sigma_1, "dunif(0.001, 2)")
  expect_match(fit_unif_floor$jags_code, "sigma_1 ~ dunif(0.001, 2)", fixed = TRUE)

  # Uniform prior with upper bound not exceeding floor errors
  expect_error(
    mcp(list(y ~ 1), data, par_x = "x", prior = list(sigma_1 = "dunif(-2, 0)"), sample = FALSE),
    "must allow values of at least 0.001"
  )

  # Deterministic expressions are not treated as distributions
  fit_expr = mcp(
    list(y ~ 1), data,
    par_x = "x", prior = list(sigma_1 = "sqrt(4)"), sample = FALSE
  )
  expect_equal(fit_expr$prior$sigma_1, "sqrt(4)")
  expect_match(fit_expr$jags_code, "sigma_1 = sqrt(4)", fixed = TRUE)
  expect_false(grepl("sigma_1 = sqrt\\(4\\) T\\(", fit_expr$jags_code))
})


good_sigma = list(
  list(y ~ 1 + sigma(1)),
  list(y ~ 1 + sigma(1 + (1 | id))),
  list(y ~ 1 + sigma(1 + (ok_id_factor || id))),
  list(y ~ 1 + sigma(x + I(x^2))),
  list(y ~ 1 + sigma(1 + sin(x))),
  list(y ~ 1,
       1 + (1|id) ~ 1 + I(x^2) + sigma(1 + x)),  # Test with varying change point and more mcp stuff
  list(y | weights(weights_ok) ~ 1 + sigma(1 + x),  # With weights
       ~ 0 + sigma(1 + x))
)

test_good(good_sigma)
