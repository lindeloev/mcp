test_that("legacy fit detection blocks pre-0.4.0 mcpfit objects", {
  data = data.frame(
    x = 1:6,
    y = c(2, 4, 3, 8, 7, 9),
    id = rep(c("a", "b"), each = 3)
  )
  fit = mcp(list(y ~ 1, ~ cp((1 | id)) + 1), data, par_x = "x", sample = FALSE)
  legacy = fit
  legacy$version = NULL

  expect_error(
    get_fit_model_tables(legacy),
    "install_github",
    fixed = TRUE
  )
})


test_that("legacy mcpfits explain how to reproduce v0.3.4", {
  legacy_fit = demo_fit
  legacy_fit$version = NULL

  expect_true(legacy_is_mcpfit(legacy_fit))
  expect_error(mcp_pars(legacy_fit), "install_github", fixed = TRUE)
  expect_error(summary(legacy_fit), "v0.4.0 bug fixes", fixed = TRUE)
})


test_that("deprecated change points on the left-hand side warn and equal cp()", {
  data = data.frame(x = 1:12, y = c(2, 4, 3, 8, 7, 9, 5, 6, 4, 8, 9, 7), id = rep(c("a", "b"), 6))
  fit_new = mcp(list(y ~ 1 + x, ~ cp(1 + (1 | id)) + 0 + x, ~ 1), data, par_x = "x", sample = FALSE)
  expect_warning(
    expect_warning(
      fit_old <- mcp(list(y ~ 1 + x, 1 + (1 | id) ~ 0 + x, y ~ 1 ~ 1), data, par_x = "x", sample = FALSE),
      "replace `1 + (1 | id) ~ 0 + x` with `~ cp(1 + (1 | id)) + 0 + x` in segment 2",
      fixed = TRUE
    ),
    "replace `y ~ 1 ~ 1` with `~ 1` in segment 3",
    fixed = TRUE
  )
  expect_identical(fit_old$jags_code, fit_new$jags_code)
  expect_identical(mcp_pars(fit_old), mcp_pars(fit_new))
  expect_identical(formula(fit_old), formula(fit_new))
})


test_that("helpful deprecation detections work for old code idioms", {
  data = data.frame(x = 1:6, y = c(2, 4, 3, 8, 7, 9))
  fit = mcp(list(y ~ 1 + x), data, sample = FALSE)

  # data = NULL or missing in mcp()
  expect_error(
    mcp(list(y ~ 1), data = NULL),
    "`data` is required in mcp() since mcp v0.4.0",
    fixed = TRUE
  )

  # plot/fitted with which_y
  expect_warning(
    plot(demo_fit, which_y = "sigma"),
    "`which_y` was deprecated in mcp v0.4.0. Use `plot_dpar()` instead.",
    fixed = TRUE
  )
  expect_warning(
    p <- plot(demo_fit, which_y = "sigma", interval = "fitted"),
    "which_y"
  )
  expect_equal(p$labels$y, "sigma")
  expect_true("quantile" %in% names(p$layers[[2]]$data))
  expect_warning(
    p <- plot(demo_fit, which_y = "ct"),
    "which_y"
  )
  expect_equal(nrow(p$layers[[1]]$data), nrow(demo_fit$data))
  expect_warning(
    expect_equal(fitted(demo_fit, which_y = "ct"), fitted(demo_fit)),
    "which_y"
  )

  # fit$simulate() signature change
  expect_error(
    fit$simulate(1:5),
    "install_github",
    fixed = FALSE
  )

  # v0.3.4 parameter names in common post-fit code
  expect_error(
    hypothesis(fit, "int_1 > 0"),
    "use `Intercept_i`",
    fixed = TRUE
  )
  expect_error(
    plot_pars(fit, pars = "int_1"),
    "install_github",
    fixed = TRUE
  )
  data_id = cbind(data, id = rep(c("a", "b"), each = 3))
  expect_error(
    mcp(list(y ~ 1, ~ cp(1 + (1 | id)) + 1), data_id, prior = list(cp_1_sd = "dnorm(0, 1) T(0, )"), par_x = "x", sample = FALSE),
    "use `cp_i_<group>_sd`",
    fixed = TRUE
  )

  # ex$fit warning
  expect_warning(
    { val = fit$fit },
    "`mcp_example()` now returns an `mcpfit` object directly",
    fixed = TRUE
  )
  expect_s3_class(val, "mcpfit")
})


test_that("cores is deprecated and ignored", {
  model = list(y ~ 1)
  data = data.frame(x = 1:3, y = 1:3)

  expect_warning(
    mcp(model, data, par_x = "x", sample = FALSE, cores = 2),
    "Setting `cores` above one no longer enables parallel processing",
    class = "lifecycle_warning_deprecated"
  )
})


test_that("nsamples is a soft-deprecated alias for ndraws", {
  lifecycle::expect_deprecated({
    result = fitted(demo_fit, nsamples = 2, summary = FALSE)
  })
  expect_equal(length(unique(result$.draw)), 2)
  expect_error(
    suppressWarnings(fitted(demo_fit, ndraws = 2, nsamples = 2)),
    "Use only one of `ndraws` and deprecated `nsamples`.",
    fixed = TRUE
  )
})


test_that("probs is a deprecated alias for central intervals", {
  lifecycle::expect_deprecated({
    result = fitted(demo_fit, probs = c(0.1, 0.9))
  })
  expect_equal(result, fitted(demo_fit, width = 0.8))
  expect_named(suppressWarnings(residuals(demo_fit, probs = FALSE)), c(colnames(demo_fit$data), ".mean", ".sd"))
  expect_error(suppressWarnings(fitted(demo_fit, probs = c(0.1, 0.5, 0.9))), "central interval")
  expect_error(
    suppressWarnings(predict(demo_fit, width = 0.8, probs = c(0.1, 0.9))),
    "Use only one of `width` and deprecated `probs`.",
    fixed = TRUE
  )
})


test_that("deprecated mcpfit accessors warn and still work", {
  lifecycle::expect_deprecated(expect_true(is.mcpfit(demo_fit)))
  lifecycle::expect_deprecated({
    val = demo_fit$mcmc_post
  })
  expect_s3_class(val, "mcmc.list")
})
