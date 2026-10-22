test_that("AR and MA terms are parsed as separate model components", {
  data = data.frame(x = 1:6, y = 1:6)
  predictors = get_predictors(
    list(y ~ 1 + x + ar(2, 1 + x) + ma(1, 0 + x)),
    data,
    as_mcpfamily(gaussian()),
    par_x = "x"
  )

  expect_setequal(unique(predictors$dpar), c("mu", "sigma", "ar", "ma"))
  expect_equal(unique(predictors$order[predictors$dpar == "ar"]), 1:2)
  expect_equal(unique(predictors$order[predictors$dpar == "ma"]), 1)
  expect_equal(sum(predictors$dpar == "ar"), 4)
  expect_equal(sum(predictors$dpar == "ma"), 1)
  expect_true(all(grepl("^ar[12]_", predictors$name[predictors$dpar == "ar"])))
  expect_true(all(grepl("^ma1_", predictors$name[predictors$dpar == "ma"])))
})


test_that("AR and MA each allow one term per segment", {
  data = data.frame(x = 1:6, y = 1:6)
  family = as_mcpfamily(gaussian())

  expect_error(
    get_predictors(list(y ~ ar(1) + ar(2)), data, family, par_x = "x"),
    "Only one `ar()` allowed in each segment", fixed = TRUE
  )
  expect_error(
    get_predictors(list(y ~ ma(1) + ma(2)), data, family, par_x = "x"),
    "Only one `ma()` allowed in each segment", fixed = TRUE
  )
  expect_no_error(
    get_predictors(list(y ~ ar(1) + ma(1)), data, family, par_x = "x")
  )
})


test_that("including a lower-order AR/MA turns off higher lags", {
  data = data.frame(x = 1:6, y = 1:6)
  family = as_mcpfamily(gaussian())

  predictors = get_predictors(
    list(y ~ ar(2), ~ ar(1)), data, family, par_x = "x"
  )
  expect_equal(
    predictors$next_segment[predictors$name == "ar2_1"],
    2L
  )

  joined_lower = get_predictors(
    list(y ~ ar(2), ~ ar(1, 0)), data, family, par_x = "x"
  )
  expect_equal(joined_lower$next_segment[joined_lower$name == "ar1_1"], NA_integer_)
  expect_equal(joined_lower$next_segment[joined_lower$name == "ar2_1"], 2L)

  fit = mcp(list(y ~ ar(2), ~ ar(1)), data, par_x = "x", sample = FALSE, quiet = TRUE)
  expect_match(
    fit$jags_code,
    "(?s)# Formula for ar2.*x\\[i_\\] < cp_1.*ar2_1",
    perl = TRUE
  )

  zeroed = get_predictors(
    list(y ~ ar(2), ~ ar(2, 0)), data, family, par_x = "x"
  )
  expect_true(all(is.na(zeroed$next_segment[zeroed$dpar == "ar"])))

  ma_predictors = get_predictors(
    list(y ~ ma(2), ~ ma(1)), data, family, par_x = "x"
  )
  expect_equal(
    ma_predictors$next_segment[ma_predictors$name == "ma2_1"],
    2L
  )
})


test_that("ar(0) and ma(0) turn off AR/MA components", {
  data = data.frame(x = 1:6, y = 1:6)
  family = as_mcpfamily(gaussian())

  zeroed_ar = get_predictors(
    list(y ~ 1 + x + ar(2), ~ 0 + x + ar(0)), data, family, par_x = "x"
  )
  expect_true(all(zeroed_ar$next_segment[zeroed_ar$dpar == "ar" & zeroed_ar$segment == 1] == 2L))
  expect_equal(sum(zeroed_ar$dpar == "ar" & zeroed_ar$segment == 2), 0)

  zeroed_ma = get_predictors(
    list(y ~ 1 + x + ma(2), ~ 0 + x + ma(0)), data, family, par_x = "x"
  )
  expect_true(all(zeroed_ma$next_segment[zeroed_ma$dpar == "ma" & zeroed_ma$segment == 1] == 2L))
  expect_equal(sum(zeroed_ma$dpar == "ma" & zeroed_ma$segment == 2), 0)

  fit_ar = mcp(list(y ~ 1 + x + ar(1), ~ 0 + x + ar(0)), data, par_x = "x", sample = FALSE, quiet = TRUE)
  expect_match(fit_ar$jags_code, "x\\[i_\\] < cp_1")
  expect_false(grepl("ar1_2", fit_ar$jags_code))

  fit_ma = mcp(list(y ~ 1 + x + ma(1), ~ 0 + x + ma(0)), data, par_x = "x", sample = FALSE, quiet = TRUE)
  expect_match(fit_ma$jags_code, "x\\[i_\\] < cp_1")
  expect_false(grepl("ma1_2", fit_ma$jags_code))
})


test_that("positive-order AR/MA zero formulas are joined to the earlier intercept and x-terms", {
  data = data.frame(x = 1:6, y = 1:6)
  fit = mcp(
    list(y ~ ar(1, 1 + x), ~ ar(1, 0)),
    data, par_x = "x", sample = FALSE, quiet = TRUE
  )
  expect_lifetimes(fit, c(ar1_1 = NA_integer_, ar1_x_1 = NA_integer_))

  ar = fit$simulate(
    fit,
    data,
    cp_1 = 3.5,
    Intercept_1 = 0, Intercept_2 = 0,
    ar1_1 = 0.2, ar1_x_1 = 0.1,
    sigma_1 = 1,
    .type = "fitted", .dpar = "ar1", .scale = "linear"
  )
  expect_equal(as.numeric(ar), 0.2 + 0.1 * pmin(data$x, 3.5))
})


test_that("AR/MA components are off until included again", {
  data = data.frame(x = 1:9, y = rep(c(1, 2, 3), 3))
  predictors = get_predictors(
    list(y ~ ar(1), ~ 0, ~ ar(1)), data, as_mcpfamily(gaussian()), par_x = "x"
  )
  expect_equal(
    stats::setNames(predictors$next_segment[predictors$dpar == "ar"], predictors$name[predictors$dpar == "ar"]),
    c(ar1_1 = 2L, ar1_3 = NA_integer_)
  )

  expect_warning(
    fit <- mcp(list(y ~ ar(1), ~ 0, ~ ar(1)), data, par_x = "x", sample = FALSE, quiet = TRUE),
    "Segment 2 does not include ar()", fixed = TRUE
  )
  expect_false("ar1_2" %in% mcp_pars(fit)$parameter)
  expect_match(fit$jags_code, "x\\[i_\\] < cp_1.*ar1_1", perl = TRUE)
  expect_match(fit$jags_code, "x\\[i_\\] >= cp_2.*ar1_3", perl = TRUE)
})


test_that("MA order errors describe MA syntax", {
  data = data.frame(x = 1:6, y = 1:6)

  expect_error(
    get_predictors(
      list(y ~ ma(x)),
      data,
      as_mcpfamily(gaussian()),
      par_x = "x"
    ),
    "Must be ma(order) or ma(order, formula)",
    fixed = TRUE
  )
})


test_that("One GARMA threshold applies to all segments and defaults to 0.1", {
  data = data.frame(x = 1:6, y = 1:6)
  family = as_mcpfamily(gaussian())
  get_thresholds = function(model) {
    predictors = get_predictors(model, data, family, par_x = "x")
    unique(predictors$threshold[predictors$dpar %in% c("ar", "ma")])
  }

  expect_equal(get_thresholds(list(y ~ ar(1), ~ ar(1, 1 + x))), 0.1)
  expect_equal(get_thresholds(list(y ~ ar(1), ~ ar(1, 1 + x, threshold = 0.2) + ma(1))), 0.2)

  expect_error(
    get_predictors(list(y ~ ar(1, threshold = 0.2) + ma(1, threshold = 0.3)), data, family, "x"),
    "Only one `threshold` value can be specified across all ar() and ma() terms.",
    fixed = TRUE
  )
  expect_error(
    get_predictors(list(y ~ ar(1, threshold = 0.2), ~ ar(1, threshold = 0.3)), data, family, "x"),
    "Only one `threshold` value can be specified across all ar() and ma() terms.",
    fixed = TRUE
  )
  expect_error(
    get_predictors(list(y ~ ar(1, threshold = 0)), data, family, "x"),
    "`threshold` in ar() must be one number between 0 and 1.",
    fixed = TRUE
  )
})


test_that("mcp builds MA parameters and JAGS code", {
  data = data.frame(x = 1:6, y = 1:6)
  fit = mcp(list(y ~ 1 + ma(1)), data, par_x = "x", sample = FALSE)

  expect_true(is_arma(fit))
  expect_match(fit$jags_code, "ma1_\\[i_\\] \\* resid_ma_\\[i_ - 1\\]")
})


test_that("GARMA evaluation keeps group-level draws as separate series", {
  data = data.frame(y = 1:5, x = -1:3, id = letters[1:5])
  fit = mcp(
    list(y ~ 1, ~ cp(1 + (1 | id)) + 1 + ar(1)),
    data, par_x = "x", sample = FALSE, quiet = TRUE
  )
  draws = cbind(
    cp_1 = c(0.5, 0.6),
    cp_1_id_sd = c(0.1, 0.2),
    Intercept_1 = c(0, 0.1),
    Intercept_2 = c(1, 1.1),
    sigma_1 = c(1, 1),
    ar1_2 = c(0.2, 0.3)
  )
  group_draws = matrix(0, nrow = 2, ncol = nrow(data))
  colnames(group_draws) = paste0("cp_1_id[", data$id, "]")
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(cbind(draws, group_draws)))

  loglik = expect_no_error(log_lik(fit))
  expect_equal(dim(loglik), c(2, 5))
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
  expect_prior(fit$prior$ar1_x_1, paste0("dnorm(0, ", format_prior_number(0.25 / stats::sd(data$x)), ")"))
  expect_equal(fit$prior$ma1_groupb_1, "dnorm(0, 0.25)")

  ar_prior = prior_summary(fit)
  ar_prior = ar_prior[ar_prior$parameter == "ar1_1", ]
  expect_equal(ar_prior$prior, "dnorm(0, 0.5) T(-1, 1)")
  expect_equal(ar_prior$bounds, "[-1, 1]")
})


test_that("AR/MA root warnings are conditional", {
  expect_equal(
    arma_root_violations(list(ma1 = matrix(c(0.5, 1.1))), "ma"),
    matrix(c(FALSE, TRUE))
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

  group_fit = mcp(
    list(y ~ 1 + ar(1, 1 + x)), data, par_x = "x", sample = FALSE
  )
  expect_warning(
    suppressMessages(group_fit$simulate(
      group_fit, data,
      Intercept_1 = 0, sigma_1 = 1, ar1_1 = 0, ar1_x_1 = 0.3
    )),
    "locally non-stationary AR"
  )
})
