test_that("observed AR residuals are isolated by series", {
  resid_abs = c(1, 2, 3, 1, 2, 3)
  result = mcp:::simulate_ar(
    sigma_ = rep(1, 6),
    ar_list = list(ar1_ = rep(0.5, 6)),
    resid_abs = resid_abs,
    series_id = rep(1:2, each = 3)
  )

  expect_equal(result$resid_ar, c(0, 0.5, 1, 0, 0.5, 1))
})


test_that("series isolation does not depend on contiguous rows", {
  result = mcp:::simulate_ar(
    sigma_ = rep(1, 6),
    ar_list = list(ar1_ = rep(0.5, 6)),
    resid_abs = c(1, 1, 2, 2, 3, 3),
    series_id = rep(1:2, 3)
  )

  expect_equal(result$resid_ar, c(0, 0, 0.5, 0.5, 1, 1))
})


test_that("fit$simulate resets fresh residual histories by series", {
  data = data.frame(
    id = rep(c("a", "b"), each = 2),
    x = c(1, 2, 1, 3),
    y = 0
  )
  fit = mcp(
    list(y ~ 1 + ar(1)), data,
    par_x = "x", series = "id", sample = FALSE, quiet = TRUE
  )

  set.seed(42)
  innovations = stats::rnorm(4)
  expected = innovations
  expected[c(2, 4)] = innovations[c(2, 4)] + 0.5 * innovations[c(1, 3)]

  set.seed(42)
  actual = suppressMessages(fit$simulate(
    fit, data, Intercept_1 = 0, sigma_1 = 1, ar1_1 = 0.5
  ))
  expect_equal(as.numeric(actual), expected)
})


test_that("AR models use all available early lags", {
  result = mcp:::simulate_ar(
    sigma_ = rep(1, 4),
    ar_list = list(ar1_ = rep(0.5, 4), ar2_ = rep(0.25, 4)),
    resid_abs = 1:4
  )

  expect_equal(result$resid_ar, c(0, 0.5, 1.25, 2))
})


test_that("early lags use coefficients from the current observation", {
  result = mcp:::simulate_ar(
    sigma_ = rep(1, 3),
    ar_list = list(ar1_ = c(0.1, 0.5, 0.9), ar2_ = c(0.01, 0.05, 0.09)),
    resid_abs = c(2, 4, 8)
  )

  expect_equal(result$resid_ar, c(0, 1, 3.78))

  jags_code = mcp:::get_ar_jagscode(2, "x")
  expect_match(jags_code, "ar1_\\[2\\] \\* resid_abs_\\[2 - 1\\]")
})


test_that("generated GARMA code uses all available early AR and MA lags", {
  jags_code = mcp:::get_garma_jagscode(1, 2, "x")

  expect_match(jags_code, "ar1_\\[2\\] \\* resid_abs_\\[2 - 1\\]")
  expect_match(jags_code, "ma1_\\[2\\] \\* resid_ma_\\[2 - 1\\]")
  expect_false(grepl("ma2_\\[2\\]", jags_code))
  expect_match(jags_code, "ma2_\\[i_\\] \\* resid_ma_\\[i_ - 2\\]")
})


test_that("generated AR residuals use the same partial-lag recursion", {
  ar_list = list(
    ar1_ = c(0.1, 0.2, 0.3, 0.4),
    ar2_ = c(0.01, 0.02, 0.03, 0.04)
  )
  set.seed(42)
  innovations = stats::rnorm(4)
  expected_residuals = numeric(4)
  expected_residuals[1] = innovations[1]
  expected_residuals[2] = innovations[2] + ar_list$ar1_[2] * expected_residuals[1]
  expected_residuals[3] = innovations[3] + ar_list$ar1_[3] * expected_residuals[2] + ar_list$ar2_[3] * expected_residuals[1]
  expected_residuals[4] = innovations[4] + ar_list$ar1_[4] * expected_residuals[3] + ar_list$ar2_[4] * expected_residuals[2]

  set.seed(42)
  result = suppressMessages(mcp:::simulate_ar(rep(1, 4), ar_list))

  expect_equal(result$resid_sigma, innovations)
  expect_equal(result$resid_ar, expected_residuals - innovations)
})


test_that("pp_eval keeps posterior draws in separate AR series", {
  data = data.frame(x = 1:3, y = 1:3)
  fit = mcp(list(y ~ 1 + ar(1)), data, par_x = "x", sample = "none")
  draws = matrix(
    c(0, 1, 0.5,
      0, 1, 0.5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("Intercept_1", "sigma_1", "ar1_1"))
  )
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))  # direct slot assignment, not the $ generic

  result = fitted(fit, summary = FALSE, probs = FALSE)

  expect_equal(result$.draw, rep(1:2, each = 3))
  expect_equal(result$.epred, rep(c(0, 0.5, 1), 2))
})
