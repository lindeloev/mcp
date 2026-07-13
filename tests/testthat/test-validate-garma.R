garma_test_fit = function(family, model, data, intercept, phi = 0.5) {
  fit = mcp(model, data, family = family, par_x = "x", sample = FALSE)
  args = stats::setNames(as.list(rep(0, length(fit$pars$population))), fit$pars$population)
  args$Intercept_1 = intercept
  args$ar1_1 = phi
  if ("sigma_1" %in% names(args))
    args$sigma_1 = 1

  list(fit = fit, args = args)
}


evaluate_garma_test_fit = function(built, data, ...) {
  predictors = add_rhs_predictors(data, built$fit)
  args = lapply(built$args, rep, nrow(data))
  do.call(
    simulate_vectorized,
    c(list(fit = built$fit), as.list(predictors), args, list(...))
  )
}


test_that("R evaluates GARMA AR on the link scale", {
  phi = 0.5
  count_data = data.frame(x = 1:4, y = c(0, 1, 4, 2))
  binomial_data = transform(count_data, N = 4)
  gaussian_data = data.frame(x = 1:4, y = c(1, 2, 3, 4))
  specs = list(
    gaussian = list(
      family = gaussian(),
      model = list(y ~ 1 + ar(1)),
      data = gaussian_data,
      intercept = 2,
      garma_y = gaussian_data$y
    ),
    binomial = list(
      family = binomial(),
      model = list(y | trials(N) ~ 1 + ar(1)),
      data = binomial_data,
      intercept = 0,
      garma_y = pmin(pmax(binomial_data$y, 0.1), binomial_data$N - 0.1) / binomial_data$N
    ),
    poisson = list(
      family = poisson(),
      model = list(y ~ 1 + ar(1)),
      data = count_data,
      intercept = log(2),
      garma_y = pmax(count_data$y, 0.1)
    ),
    negbinomial = list(
      family = negbinomial(),
      model = list(y ~ 1 + ar(1)),
      data = count_data,
      intercept = log(2),
      garma_y = pmax(count_data$y, 0.1)
    )
  )

  for (spec in specs) {
    built = garma_test_fit(spec$family, spec$model, spec$data, spec$intercept, phi)
    expected_link = spec$intercept + phi * c(
      0,
      head(spec$family$linkfun(spec$garma_y) - spec$intercept, -1)
    )
    expected = spec$family$linkinv(expected_link)
    actual = evaluate_garma_test_fit(
      built, spec$data,
      .type = "fitted", .rate = spec$family$family == "binomial"
    )
    expect_equal(as.numeric(actual), expected)
    actual_linear = evaluate_garma_test_fit(
      built, spec$data,
      .type = "fitted", .scale = "linear"
    )
    expect_equal(as.numeric(actual_linear), expected_link)
  }
})


test_that("R-side GARMA log likelihood uses its conditional mean", {
  data = data.frame(x = 1:4, y = c(0, 1, 4, 2))
  built = garma_test_fit(poisson(), list(y ~ 1 + ar(1)), data, log(2))
  expected_link = log(2) + 0.5 * c(0, head(log(pmax(data$y, 0.1)) - log(2), -1))
  expected = stats::dpois(data$y, exp(expected_link), log = TRUE)
  actual = evaluate_garma_test_fit(built, data, .type = "loglik")
  expect_equal(as.numeric(actual), expected)
})


test_that("fit$simulate generates GARMA series recursively", {
  data = data.frame(x = 1:6, y = 0)
  built = garma_test_fit(poisson(), list(y ~ 1 + ar(1)), data, log(2))

  set.seed(42)
  expected = numeric(nrow(data))
  residual = numeric(nrow(data))
  for (i in seq_len(nrow(data))) {
    lagged = if (i == 1) 0 else 0.5 * residual[i - 1]
    expected[i] = stats::rpois(1, exp(log(2) + lagged))
    residual[i] = log(max(expected[i], 0.1)) - log(2)
  }

  set.seed(42)
  call = c(
    list(fit = built$fit, newdata = data),
    built$args,
    list(.type = "predict")
  )
  actual = do.call(built$fit$simulate, call)

  expect_equal(as.numeric(actual), expected)
})


test_that("generated JAGS uses the same bounded link-scale residuals", {
  count_data = data.frame(x = 1:4, y = c(0, 1, 4, 2))
  binomial_data = transform(count_data, N = 4)
  poisson_fit = mcp(
    list(y ~ 1 + ar(1, boundary = 0.2)),
    count_data,
    family = poisson(),
    par_x = "x",
    sample = FALSE
  )
  binomial_fit = mcp(
    list(y | trials(N) ~ 1 + ar(1)),
    binomial_data,
    family = binomial(),
    par_x = "x",
    sample = FALSE
  )
  segmented_fit = mcp(
    list(
      y ~ 1 + ar(1),
      ~ 1 + ar(1, boundary = 0.2)
    ),
    count_data,
    family = poisson(),
    par_x = "x",
    sample = FALSE
  )

  expect_match(poisson_fit$jags_code, "garma_boundary_\\[i_\\] =\\s+0\\.2")
  expect_match(poisson_fit$jags_code, "garma_y_\\[i_\\] = max\\(y\\[i_\\], garma_boundary_\\[i_\\]\\)")
  expect_match(poisson_fit$jags_code, "resid_abs_\\[i_\\] = log\\(garma_y_\\[i_\\]\\) - link_mu_\\[i_\\]")
  expect_match(binomial_fit$jags_code, "garma_y_\\[i_\\] = min\\(max\\(y\\[i_\\], garma_boundary_\\[i_\\]\\), N\\[i_\\] - garma_boundary_\\[i_\\]\\) / N\\[i_\\]")
  expect_match(binomial_fit$jags_code, "resid_abs_\\[i_\\] = logit\\(garma_y_\\[i_\\]\\) - link_mu_\\[i_\\]")
  expect_match(segmented_fit$jags_code, "\\(x\\[i_\\] < cp_1\\) \\* 0\\.1")
  expect_match(segmented_fit$jags_code, "\\(x\\[i_\\] >= cp_1\\) \\* 0\\.2")
})
