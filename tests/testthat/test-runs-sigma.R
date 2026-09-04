#################
# TEST GAUSSIAN STANDARD DEVIATION #
#################
bad_sigma = list(
  list(y ~ 1 + sigma(q))  # variable does not exist
)

test_bad(bad_sigma)


test_that("fixed residual SDs must be positive", {
  data = data.frame(x = 1:4, y = 0)

  expect_error(
    mcp(
      list(y ~ 1), data,
      par_x = "x", prior = list(sigma_1 = 0), sample = FALSE
    ),
    "Fixed residual standard deviation parameter(s) must be positive: sigma_1.",
    fixed = TRUE
  )
  expect_error(
    mcp(
      list(y ~ 1), data,
      par_x = "x", prior = list(sigma_1 = -1), sample = FALSE
    ),
    "Fixed residual standard deviation parameter(s) must be positive: sigma_1.",
    fixed = TRUE
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
