garma_test_fit = function(family, model, data, intercept, phi = 0.5, theta = 0.25) {
  fit = mcp(model, data, family = family, par_x = "x", sample = FALSE)
  population = mcp_pars(fit, scope = "population")$name
  args = stats::setNames(as.list(rep(0, length(population))), population)
  args$Intercept_1 = intercept
  if ("ar1_1" %in% names(args))
    args$ar1_1 = phi
  if ("ma1_1" %in% names(args))
    args$ma1_1 = theta
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


expected_garma_link = function(link_y, base_link_mu, phi = 0, theta = 0) {
  link_mu = numeric(length(link_y))
  resid_abs = link_y - base_link_mu
  resid_ma = numeric(length(link_y))
  for (i in seq_along(link_y)) {
    lagged_ar = if (i == 1) 0 else phi * resid_abs[i - 1]
    lagged_ma = if (i == 1) 0 else theta * resid_ma[i - 1]
    link_mu[i] = base_link_mu + lagged_ar + lagged_ma
    resid_ma[i] = link_y[i] - link_mu[i]
  }
  link_mu
}


test_that("R evaluates combined AR and MA on the link scale", {
  phi = 0.5
  theta = 0.25
  count_data = data.frame(x = 1:4, y = c(0, 1, 4, 2))
  binomial_data = transform(count_data, N = 4)
  gaussian_data = data.frame(x = 1:4, y = c(1, 2, 3, 4))
  specs = list(
    gaussian = list(
      family = gaussian(),
      model = list(y ~ 1 + ar(1) + ma(1)),
      data = gaussian_data,
      intercept = 2,
      garma_y = gaussian_data$y
    ),
    binomial = list(
      family = binomial(),
      model = list(y | trials(N) ~ 1 + ar(1) + ma(1)),
      data = binomial_data,
      intercept = 0,
      garma_y = pmin(pmax(binomial_data$y, 0.1), binomial_data$N - 0.1) / binomial_data$N
    ),
    poisson = list(
      family = poisson(),
      model = list(y ~ 1 + ar(1) + ma(1)),
      data = count_data,
      intercept = log(2),
      garma_y = pmax(count_data$y, 0.1)
    ),
    negbinomial = list(
      family = negbinomial(),
      model = list(y ~ 1 + ar(1) + ma(1)),
      data = count_data,
      intercept = log(2),
      garma_y = pmax(count_data$y, 0.1)
    )
  )

  for (spec in specs) {
    built = garma_test_fit(spec$family, spec$model, spec$data, spec$intercept, phi, theta)
    expected_link = expected_garma_link(
      spec$family$linkfun(spec$garma_y), spec$intercept, phi, theta
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
  built = garma_test_fit(poisson(), list(y ~ 1 + ar(1) + ma(1)), data, log(2))
  expected_link = expected_garma_link(log(pmax(data$y, 0.1)), log(2), 0.5, 0.25)
  expected = stats::dpois(data$y, exp(expected_link), log = TRUE)
  actual = evaluate_garma_test_fit(built, data, .type = "loglik")
  expect_equal(as.numeric(actual), expected)
})


test_that("conditional GARMA evaluation resets vectorized ordered series", {
  y = c(1, 2, 4, 10, 20, 40)
  series_id = rep(c("a", "b"), each = 3)
  base_link_mu = rep(0, length(y))
  expected = numeric(length(y))
  for (series in unique(series_id)) {
    rows = which(series_id == series)
    expected[rows] = expected_garma_link(y[rows], 0, phi = 0.5, theta = 0.25)
  }

  args = list(
    base_link_mu = base_link_mu,
    ar_list = list(ar1_ = rep(0.5, length(y))),
    ma_list = list(ma1_ = rep(0.25, length(y))),
    boundary = rep(0.1, length(y)),
    family = mcpfamily(gaussian()),
    dpars = list(mu = base_link_mu, sigma = rep(1, length(y))),
    y = y,
    series_id = series_id
  )
  actual = do.call(simulate_garma, args)

  expect_equal(actual$link_mu, expected)
  args$series_id = rep(c("a", "b"), 3)
  expect_error(
    do.call(simulate_garma, args),
    "Rows belonging to each series_id must be contiguous.",
    fixed = TRUE
  )
})


test_that("fit$simulate generates GARMA series recursively", {
  data = data.frame(x = 1:6, y = 0)
  built = garma_test_fit(poisson(), list(y ~ 1 + ar(1) + ma(1)), data, log(2))

  set.seed(42)
  expected = numeric(nrow(data))
  resid_abs = numeric(nrow(data))
  resid_ma = numeric(nrow(data))
  for (i in seq_len(nrow(data))) {
    lagged = if (i == 1) 0 else 0.5 * resid_abs[i - 1] + 0.25 * resid_ma[i - 1]
    expected[i] = stats::rpois(1, exp(log(2) + lagged))
    link_y = log(max(expected[i], 0.1))
    resid_abs[i] = link_y - log(2)
    resid_ma[i] = link_y - log(2) - lagged
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


test_that("posterior_predict generates fresh GARMA series recursively", {
  data = data.frame(x = 1:3, y = c(10, 20, 30))
  fit = mcp(
    list(y ~ 1 + ar(1)), data,
    par_x = "x", sample = FALSE, quiet = TRUE
  )
  fixed_draws = matrix(
    rep(c(0, 0, 0.5), 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("Intercept_1", "sigma_1", "ar1_1"))
  )
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(fixed_draws))

  set.seed(42)
  conditional = predict(fit, summary = FALSE, probs = FALSE)
  set.seed(42)
  replicated = posterior_predict.mcpfit(fit)
  set.seed(42)
  predictor_only = posterior_predict.mcpfit(
    fit,
    newdata = data.frame(x = data$x)
  )

  expect_equal(conditional$.prediction, rep(c(0, 5, 10), 2), tolerance = 0.01)
  expect_equal(unname(replicated), matrix(0, nrow = 2, ncol = 3), tolerance = 0.01)
  expect_equal(predictor_only, replicated)

  novel_missing = data
  novel_missing$y[2] = NA
  for (method in list(fitted, predict, residuals)) {
    expect_error(
      method(fit, newdata = novel_missing),
      "supported only for the original fitted data",
      fixed = TRUE
    )
  }

  fit$data$y[2] = NA
  set.seed(42)
  expect_equal(posterior_predict.mcpfit(fit), replicated)
})


test_that("generated JAGS uses the same bounded GARMA residuals", {
  count_data = data.frame(x = 1:4, y = c(0, 1, 4, 2))
  binomial_data = transform(count_data, N = 4)
  poisson_fit = mcp(
    list(y ~ 1 + ar(1, boundary = 0.2) + ma(1)),
    count_data,
    family = poisson(),
    par_x = "x",
    sample = FALSE
  )
  binomial_fit = mcp(
    list(y | trials(N) ~ 1 + ar(1) + ma(1)),
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

  expect_match(poisson_fit$jags_code, "garma_y_\\[i_\\] = max\\(y\\[i_\\], garma_boundary_\\[i_\\]\\)")
  expect_match(poisson_fit$jags_code, "garma_link_y_\\[i_\\] = log\\(garma_y_\\[i_\\]\\)")
  expect_match(poisson_fit$jags_code, "resid_abs_\\[i_\\] = garma_link_y_\\[i_\\] - link_mu_\\[i_\\]")
  expect_match(poisson_fit$jags_code, "resid_ma_\\[i_\\] = garma_link_y_\\[i_\\] - link_mu_\\[i_\\] - resid_garma_\\[i_\\]")
  expect_match(poisson_fit$jags_code, "ma1_\\[i_\\] \\* resid_ma_\\[i_ - 1\\]")
  expect_match(binomial_fit$jags_code, "garma_y_\\[i_\\] = min\\(max\\(y\\[i_\\], garma_boundary_\\[i_\\]\\), N\\[i_\\] - garma_boundary_\\[i_\\]\\) / N\\[i_\\]")
  expect_match(binomial_fit$jags_code, "garma_link_y_\\[i_\\] = logit\\(garma_y_\\[i_\\]\\)")
  expect_match(segmented_fit$jags_code, "\\(x\\[i_\\] < cp_1\\) \\* 0\\.1")
  expect_match(segmented_fit$jags_code, "\\(x\\[i_\\] >= cp_1\\) \\* 0\\.2")
})


test_that("missing GARMA responses stay paired with posterior draws", {
  data = data.frame(x = 1:5, y = c(1, NA, 2, 3, 4))
  fit = suppressWarnings(mcp(
    list(y ~ 1 + ar(1)), data, par_x = "x",
    chains = 1, iter = 30, warmup = 20, quiet = TRUE, seed = 42
  ))

  expect_equal(fit$.internal$imputed_response_rows, 2L)
  expect_equal(ncol(fit$.internal$imputed_response[[1]]), 1L)
  expect_false("y[2]" %in% colnames(.subset2(fit, "mcmc_post")[[1]]))

  fitted_draws = fitted(fit, summary = FALSE, probs = FALSE)
  fitted_again = fitted(fit, summary = FALSE, probs = FALSE)
  prediction_draws = predict(fit, summary = FALSE, probs = FALSE)
  imputed = get_imputed_response_draws(fit, prediction_draws)

  expect_equal(fitted_draws$.epred, fitted_again$.epred)
  expect_true(all(is.na(fitted_draws$y[fitted_draws$data_row == 2])))
  expect_equal(
    prediction_draws$.prediction[prediction_draws$data_row == 2],
    imputed[prediction_draws$data_row == 2]
  )

  row3 = fitted_draws$data_row == 3
  imputed_by_draw = prediction_draws$.prediction[
    prediction_draws$data_row == 2
  ][match(
    fitted_draws$.draw[row3],
    prediction_draws$.draw[prediction_draws$data_row == 2]
  )]
  expected_row3 = fitted_draws$Intercept_1[row3] +
    fitted_draws$ar1_1[row3] * (imputed_by_draw - fitted_draws$Intercept_1[row3])
  expect_equal(fitted_draws$.epred[row3], expected_row3)

  fitted_summary = fitted(fit)
  prediction_summary = predict(fit)
  expect_false(is.na(fitted_summary$fitted[2]))
  expect_false(is.na(prediction_summary$predict[2]))
  expect_true(is.na(residuals(fit)$residuals[2]))

  expect_error(log_lik(fit), "missing response occurs in the history")
  expect_error(loo(fit), "missing response occurs in the history")
  expect_error(waic(fit), "missing response occurs in the history")
  expect_equal(ncol(log_lik(fit, arma = FALSE)), 4L)
})
