test_that("AR and MA terms are parsed as separate model components", {
  data = data.frame(x = 1:6, y = 1:6)
  predictors = get_predictors(
    list(y ~ 1 + x + ar(2, 1 + x) + ma(1, 0 + x)),
    data,
    mcpfamily(gaussian()),
    par_x = "x"
  )

  expect_setequal(unique(predictors$dpar), c("mu", "sigma", "ar", "ma"))
  expect_equal(unique(predictors$order[predictors$dpar == "ar"]), 1:2)
  expect_equal(unique(predictors$order[predictors$dpar == "ma"]), 1)
  expect_equal(sum(predictors$dpar == "ar"), 4)
  expect_equal(sum(predictors$dpar == "ma"), 1)
  expect_true(all(grepl("^ar[12]_", predictors$code_name[predictors$dpar == "ar"])))
  expect_true(all(grepl("^ma1_", predictors$code_name[predictors$dpar == "ma"])))
})


test_that("AR and MA each allow one term per segment", {
  data = data.frame(x = 1:6, y = 1:6)
  family = mcpfamily(gaussian())

  expect_error(
    get_predictors(list(y ~ ar(1) + ar(2)), data, family, par_x = "x"),
    "Only one of these allowed per segment"
  )
  expect_error(
    get_predictors(list(y ~ ma(1) + ma(2)), data, family, par_x = "x"),
    "Only one of these allowed per segment"
  )
  expect_no_error(
    get_predictors(list(y ~ ar(1) + ma(1)), data, family, par_x = "x")
  )
})


test_that("a lower-order AR/MA declaration turns off higher lags", {
  data = data.frame(x = 1:6, y = 1:6)
  family = mcpfamily(gaussian())

  predictors = get_predictors(
    list(y ~ ar(2), ~ ar(1)), data, family, par_x = "x"
  )
  expect_equal(
    predictors$next_intercept[predictors$code_name == "ar2_1"],
    2L
  )

  fit = mcp(list(y ~ ar(2), ~ ar(1)), data, par_x = "x", sample = FALSE, quiet = TRUE)
  expect_match(
    fit$jags_code,
    "(?s)# Formula for ar2.*x\\[i_\\] < cp_1.*ar2_1",
    perl = TRUE
  )

  zeroed = get_predictors(
    list(y ~ ar(2), ~ ar(2, 0)), data, family, par_x = "x"
  )
  expect_true(all(zeroed$next_intercept[zeroed$dpar == "ar"] == 2L))

  ma_predictors = get_predictors(
    list(y ~ ma(2), ~ ma(1)), data, family, par_x = "x"
  )
  expect_equal(
    ma_predictors$next_intercept[ma_predictors$code_name == "ma2_1"],
    2L
  )
})


test_that("MA order errors describe MA syntax", {
  data = data.frame(x = 1:6, y = 1:6)

  expect_error(
    get_predictors(
      list(y ~ ma(x)),
      data,
      mcpfamily(gaussian()),
      par_x = "x"
    ),
    "Must be ma(order) or ma(order, formula)",
    fixed = TRUE
  )
})


test_that("GARMA boundaries default to 0.1 and can vary by segment", {
  data = data.frame(x = 1:6, y = 1:6)
  predictors = get_predictors(
    list(
      y ~ ar(1),
      ~ ar(1, 1 + x, boundary = 0.2)
    ),
    data,
    mcpfamily(gaussian()),
    par_x = "x"
  )

  boundaries = predictors %>%
    dplyr::filter(.data$dpar == "ar") %>%
    dplyr::distinct(.data$segment, .data$boundary)
  expect_equal(boundaries$boundary, c(0.1, 0.2))

  expect_error(
    get_predictors(list(y ~ ar(1, boundary = 0)), data, mcpfamily(gaussian()), "x"),
    "`boundary` in ar() must be one number between 0 and 1.",
    fixed = TRUE
  )
})


test_that("AR and MA share one boundary within a segment", {
  data = data.frame(x = 1:6, y = 1:6)
  family = mcpfamily(gaussian())

  predictors = get_predictors(
    list(y ~ ar(1, boundary = 0.2) + ma(1)),
    data,
    family,
    par_x = "x"
  )
  expect_equal(unique(predictors$boundary[predictors$dpar == "ar"]), 0.2)
  expect_equal(unique(predictors$boundary[predictors$dpar == "ma"]), 0.2)

  expect_error(
    get_predictors(
      list(y ~ ar(1, boundary = 0.2) + ma(1, boundary = 0.3)),
      data,
      family,
      par_x = "x"
    ),
    "ar() and ma() must use the same `boundary` within a segment.",
    fixed = TRUE
  )
})


test_that("mcp builds MA parameters and JAGS code", {
  data = data.frame(x = 1:6, y = 1:6)
  fit = mcp(list(y ~ 1 + ma(1)), data, par_x = "x", sample = FALSE)

  expect_true("ma1_1" %in% fit$pars$arma)
  expect_match(fit$jags_code, "ma1_\\[i_\\] \\* resid_ma_\\[i_ - 1\\]")
})


test_that("AR and MA defaults regularize direct coefficients", {
  data = data.frame(
    x = 1:12,
    y = rep(0, 12),
    group = rep(c("a", "b"), 6)
  )
  fit = mcp(
    list(y ~ 1 + ar(2, 1 + x) + ma(1, 1 + group)),
    data, sample = FALSE
  )

  expect_equal(fit$prior$ar1_1, "dnorm(0, 0.5) T(-1, 1)")
  expect_equal(fit$prior$ar2_1, "dnorm(0, 0.5) T(-1, 1)")
  expect_equal(fit$prior$ma1_1, "dnorm(0, 0.5) T(-1, 1)")
  expect_equal(fit$prior$ar1_x_1, "dnorm(0, 0.02272727)")
  expect_equal(fit$prior$ma1_groupb_1, "dnorm(0, 0.25)")

  ar_prior = prior_summary(fit)
  ar_prior = ar_prior[ar_prior$parameter == "ar1_1", ]
  expect_match(ar_prior$prior, "normal", fixed = TRUE)
  expect_equal(ar_prior$bounds, "[-1, 1]")
})


test_that("AR/MA root warnings are conditional", {
  expect_equal(
    arma_root_violations(data.frame(ma1_ = c(0.5, 1.1)), "ma"),
    c(FALSE, TRUE)
  )

  data = data.frame(x = 1:6, y = rep(0, 6))
  fit = mcp(list(y ~ 1 + ar(2)), data, par_x = "x", sample = FALSE)
  draws = cbind(
    Intercept_1 = 0,
    sigma_1 = 1,
    ar1_1 = c(0.5, 0.9, 0.5, 0.9),
    ar2_1 = c(0.2, 0.5, 0.2, 0.5)
  )
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  expect_warning(warn_arma_fit(fit, ndraws = 4), "AR: 50%")
  expect_no_warning(warn_arma_fit(fit, ndraws = 4, diagnostics = list(ar = 0.6)))
  expect_no_warning(warn_arma_fit(fit, ndraws = 4, diagnostics = FALSE))

  mcmc_tmp = .subset2(fit, "mcmc_post")
  mcmc_tmp[[1]][, "ar1_1"] = 0.5
  mcmc_tmp[[1]][, "ar2_1"] = 0.2
  fit$mcmc_post = mcmc_tmp
  expect_no_warning(warn_arma_fit(fit, ndraws = 4))

  expect_no_warning(
    suppressMessages(fit$simulate(
      fit, data,
      Intercept_1 = 0, sigma_1 = 1, ar1_1 = 0.5, ar2_1 = 0.2
    ))
  )
  expect_warning(
    suppressMessages(fit$simulate(
      fit, data,
      Intercept_1 = 0, sigma_1 = 1, ar1_1 = 0.9, ar2_1 = 0.5
    )),
    "locally non-stationary AR"
  )

  varying_fit = mcp(
    list(y ~ 1 + ar(1, 1 + x)), data, par_x = "x", sample = FALSE
  )
  expect_warning(
    suppressMessages(varying_fit$simulate(
      varying_fit, data,
      Intercept_1 = 0, sigma_1 = 1, ar1_1 = 0, ar1_x_1 = 0.3
    )),
    "locally non-stationary AR"
  )
})
