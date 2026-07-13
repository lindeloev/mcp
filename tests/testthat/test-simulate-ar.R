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
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  result = fitted(fit, summary = FALSE, probs = FALSE)

  expect_equal(result$.draw, rep(1:2, each = 3))
  expect_equal(result$fitted, rep(c(0, 0.5, 1), 2))
})
