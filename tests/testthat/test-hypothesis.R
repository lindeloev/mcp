# Test hypothesis() on demo_fit without requiring prior MCMC sampling
test_that("hypothesis()", {
  # Create a posterior-only fit by stripping prior draws from demo_fit
  fit_post_only = demo_fit
  fit_post_only$mcmc_prior = NULL

  # Use a draw-derived threshold so the hypothesis is neither rare nor certain.
  raw = .subset2(fit_post_only, "mcmc_post")
  cp_draws = unlist(lapply(raw, function(chain) chain[, "cp_1"]))
  threshold = unname(stats::quantile(cp_draws, 0.25))
  threshold_text = format(threshold, digits = 16)
  directional = paste0("cp_1 > ", threshold_text)

  # Directional hypothesis on posterior-only fit returns posterior probability with BF = NA
  expect_message(res_post_only <- hypothesis(fit_post_only, directional), "need prior draws, so `BF` is NA", fixed = TRUE)
  expect_true("prob" %in% names(res_post_only))
  expect_equal(res_post_only$prob, mean(cp_draws > threshold))
  expect_true(is.na(res_post_only$BF))
  expect_equal(res_post_only$mean, mean(cp_draws - threshold))
  expect_equal(res_post_only$sd, stats::sd(cp_draws))

  # Equality hypothesis requires prior draws
  expect_error(
    hypothesis(fit_post_only, paste0("cp_1 = ", threshold_text)),
    "Both prior and posterior draws are needed",
    fixed = TRUE
  )

  fit_asymmetric = fit_post_only
  fit_asymmetric$mcmc_prior = .subset2(fit_post_only, "mcmc_post")

  # Force the prior probability above the threshold to 0.25. This makes the
  # prior odds differ from the posterior odds.
  mcmc_prior = .subset2(fit_asymmetric, "mcmc_prior")
  for (chain in seq_along(mcmc_prior)) {
    n_draws = nrow(mcmc_prior[[chain]])
    mcmc_prior[[chain]][, "cp_1"] = c(
      rep(threshold + 1, floor(n_draws / 4)),
      rep(threshold - 1, n_draws - floor(n_draws / 4))
    )
  }
  fit_asymmetric$mcmc_prior = mcmc_prior

  actual_directional = expect_no_warning(hypothesis(fit_asymmetric, directional))
  p_post = mean(cp_draws > threshold)
  prior_draws = unlist(lapply(.subset2(fit_asymmetric, "mcmc_prior"), function(chain) chain[, "cp_1"]))
  p_prior = mean(prior_draws > threshold)
  expected_BF = (p_post / (1 - p_post)) / (p_prior / (1 - p_prior))
  effect_draws = cp_draws - threshold
  expect_equal(actual_directional$lower, unname(quantile(effect_draws, 0.025)))
  expect_equal(actual_directional$upper, unname(quantile(effect_draws, 0.975)))
  expect_equal(actual_directional$prob, p_post)
  expect_equal(actual_directional$BF, expected_BF)

  # Hypotheses implied or excluded by the prior have no Bayes factor
  expect_warning(implied <- hypothesis(demo_fit, "cp_2 > cp_1"), "Bayes factor is undefined")
  expect_equal(implied$prob, 1)
  expect_true(is.na(implied$BF))
  expect_equal(hypothesis(demo_fit, "cp_1 > 10")$BF, Inf)
  far_equality = suppressWarnings(hypothesis(demo_fit, "time_2 = 0"))  # No posterior draws below 0
  expect_equal(far_equality$BF, 0)

  prior_directional = hypothesis(fit_asymmetric, directional, prior = TRUE)
  expect_equal(prior_directional$prob, p_prior)
  expect_true(is.na(prior_directional$BF))

  fit_prior_only = fit_asymmetric
  fit_prior_only$mcmc_post = NULL
  expect_equal(hypothesis(fit_prior_only, directional, prior = TRUE)$prob, p_prior)

  # Identical prior and posterior draws must give BF = 1, also for intervals and Savage-Dickey equality.
  fit_same = fit_post_only
  fit_same$mcmc_prior = .subset2(fit_post_only, "mcmc_post")
  bounds = stats::quantile(cp_draws, c(0.2, 0.8))
  interval = paste0(
    "cp_1 > ", format(bounds[[1]], digits = 16),
    " & cp_1 < ", format(bounds[[2]], digits = 16)
  )
  actual_interval = hypothesis(fit_same, interval)
  expect_equal(actual_interval$prob, mean(cp_draws > bounds[[1]] & cp_draws < bounds[[2]]))
  expect_equal(actual_interval$BF, 1)

  # Savage-Dickey point equality test (requires prior)
  mid_val = format(mean(cp_draws), digits = 16)
  equality_expr = paste0("cp_1 = ", mid_val)
  expect_warning(
    actual_equality <- hypothesis(fit_same, equality_expr),
    "Savage-Dickey Bayes factor was computed using default prior(s) for `cp_1`",
    fixed = TRUE
  )
  expect_s3_class(actual_equality, "data.frame")
  expect_equal(actual_equality$hypothesis, paste0("cp_1 - ", mid_val, " = 0"))
  expect_true(is.na(actual_equality$prob))
  expect_false(is.na(actual_equality$BF))
  expect_equal(actual_equality$BF, 1, tolerance = 1e-3)

  # With user-specified prior, no default prior warning is emitted
  fit_user_prior = fit_same
  fit_user_prior$.internal$prior_table$source[fit_user_prior$.internal$prior_table$parameter == "cp_1"] = "user"
  expect_no_warning(hypothesis(fit_user_prior, equality_expr))

  prior_equality = hypothesis(fit_same, equality_expr, prior = TRUE)
  expect_true(is.na(prior_equality$BF))

  tail_val = format(stats::quantile(cp_draws, 0.005), digits = 16)
  expect_warning(
    hypothesis(fit_user_prior, paste0("cp_1 = ", tail_val)),
    "tested value is in a sparse tail",
    fixed = TRUE
  )

  # A legitimate boundary test also needs checking: draws cannot establish the bounds.
  fit_bound = demo_fit
  fit_bound$.internal$prior_table$source[fit_bound$.internal$prior_table$parameter == "sigma_1"] = "user"
  expect_warning(hypothesis(fit_bound, "sigma_1 = 0"), "at or outside the prior bounds", fixed = TRUE)

  # Kernel density is positive just outside the bounds, even with dense nearby draws.
  fit_bounded = fit_user_prior
  bounded_draws = .subset2(fit_bounded, "mcmc_prior")
  for (chain in seq_along(bounded_draws)) {
    bounded_draws[[chain]][, "cp_1"] = seq(0.001, 0.999, length.out = nrow(bounded_draws[[chain]]))
  }
  fit_bounded$mcmc_prior = bounded_draws
  fit_bounded$mcmc_post = bounded_draws
  for (value in c(-0.001, 1.001)) {
    expect_warning(
      result <- hypothesis(fit_bounded, paste0("cp_1 = ", value)),
      "misleading Bayes factor for an impossible value", fixed = TRUE
    )
    expect_equal(result$BF, 1)
  }

  # Degenerate hypotheses with constant prior contrasts are rejected
  expect_error(
    hypothesis(fit_same, "Intercept_1 = Intercept_1"),
    "prior contrast is constant",
    fixed = TRUE
  )

  fit_constant = fit_same
  mcmc_prior = .subset2(fit_constant, "mcmc_prior")
  mcmc_post = .subset2(fit_constant, "mcmc_post")
  for (chain in seq_along(mcmc_prior)) {
    mcmc_prior[[chain]][, "sigma_1"] = 1
    mcmc_post[[chain]][, "sigma_1"] = 1
  }
  fit_constant$mcmc_prior = mcmc_prior
  fit_constant$mcmc_post = mcmc_post
  expect_error(
    hypothesis(fit_constant, "sigma_1 = 0.99"),
    "prior contrast is constant",
    fixed = TRUE
  )

  fit_same_draws = fit_same
  mcmc_prior = .subset2(fit_same_draws, "mcmc_prior")
  mcmc_post = .subset2(fit_same_draws, "mcmc_post")
  for (chain in seq_along(mcmc_prior)) {
    mcmc_prior[[chain]][, "Intercept_3"] = mcmc_prior[[chain]][, "Intercept_1"]
    mcmc_post[[chain]][, "Intercept_3"] = mcmc_post[[chain]][, "Intercept_1"]
  }
  fit_same_draws$mcmc_prior = mcmc_prior
  fit_same_draws$mcmc_post = mcmc_post
  expect_error(
    hypothesis(fit_same_draws, "Intercept_1 = Intercept_3"),
    "prior contrast is constant",
    fixed = TRUE
  )

  # Inequality is not an equality test
  expect_error(hypothesis(demo_fit, "cp_1 != 30"), "`!=` and `==` are not supported", fixed = TRUE)
})


test_that("Savage-Dickey hypotheses are affine", {
  parameters = c("x", "y", "group[x]")

  expect_true(validate_savage_dickey_expression("x = 1", parameters))
  expect_true(validate_savage_dickey_expression("2 * x - y = 3", parameters))
  expect_true(validate_savage_dickey_expression("`group[x]` - x = 0", parameters))

  nonlinear = c(
    "x / y = 1",
    "x * y = 1",
    "x^2 = 1",
    "exp(x) = 1",
    "exp(y) * (x - 1) = 0"
  )
  for (expression in nonlinear) {
    expect_error(
      validate_savage_dickey_expression(expression, parameters),
      "named scalar parameter or affine contrast",
      fixed = TRUE
    )
  }
})


test_that("Savage-Dickey density is evaluated directly", {
  x = seq(-3, 3, length.out = 101)
  bandwidth = stats::bw.SJ(x)

  expect_equal(
    get_density(x, 0),
    mean(stats::dnorm(0, mean = x, sd = bandwidth))
  )
  expect_gte(get_density(x, 10), 0)
  expect_false(is_sparse_tail(x, 0))
  expect_true(is_sparse_tail(x, 3))
  expect_error(get_density(rep(1, 101), 0), "Cannot estimate density for constant values", fixed = TRUE)
})
