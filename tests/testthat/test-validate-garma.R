garma_test_fit = function(family, model, data, intercept, phi = 0.5, theta = 0.25) {
  fit = mcp(model, data, family = family, par_x = "x", sample = FALSE)
  population = mcp_pars(fit, scope = "population")$parameter
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
  model_data = get_model_data(built$fit, data)
  evaluate_model(built$fit, model_data, built$args, n_draws = 1, y = data$y, ...)$value
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
      type = "fitted", rate = spec$family$family == "binomial"
    )
    expect_equal(as.numeric(actual), expected)
    actual_linear = evaluate_garma_test_fit(
      built, spec$data,
      type = "fitted", scale = "linear"
    )
    expect_equal(as.numeric(actual_linear), expected_link)
  }
})


test_that("R-side GARMA log likelihood uses its conditional mean", {
  data = data.frame(x = 1:4, y = c(0, 1, 4, 2))
  built = garma_test_fit(poisson(), list(y ~ 1 + ar(1) + ma(1)), data, log(2))
  expected_link = expected_garma_link(log(pmax(data$y, 0.1)), log(2), 0.5, 0.25)
  expected = stats::dpois(data$y, exp(expected_link), log = TRUE)
  actual = evaluate_garma_test_fit(built, data, type = "loglik")
  expect_equal(as.numeric(actual), expected)
})


test_that("conditional GARMA evaluation resets vectorized ordered series", {
  y = c(1, 2, 4, 10, 20, 40)
  series = rep(c("a", "b"), each = 3)
  expected = numeric(length(y))
  for (this_series in unique(series)) {
    rows = which(series == this_series)
    expected[rows] = expected_garma_link(y[rows], 0, phi = 0.5, theta = 0.25)
  }

  # Two identical draws (columns)
  cells = function(value) matrix(value, length(y), 2)
  args = list(
    link_mu = cells(0),
    ar = list(ar1 = cells(0.5)),
    ma = list(ma1 = cells(0.25)),
    threshold = 0.1,
    family = as_mcpfamily(gaussian()),
    dpars = list(mu = cells(0), sigma = cells(1)),
    y = y,
    series = series
  )
  actual = do.call(simulate_garma, args)

  expect_equal(actual$link_mu, cells(expected))
  args$series = rep(c("a", "b"), 3)
  expect_error(
    do.call(simulate_garma, args),
    "Rows belonging to each series must be contiguous.",
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


test_that("prediction interfaces share conditional and replicated GARMA behavior", {
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
  conditional = predict(fit, summary = FALSE, width = FALSE)
  set.seed(42)
  conditional_matrix = posterior_predict.mcpfit(fit)
  expect_equal(unname(conditional_matrix), matrix(conditional$.prediction, nrow = 2))
  expect_equal(posterior_epred.mcpfit(fit), fitted(fit, summary = FALSE, draws_format = "matrix"))
  set.seed(42)
  replicated = posterior_predict.mcpfit(fit, arma = "replicate")
  set.seed(42)
  expect_equal(predict(fit, arma = "replicate", summary = FALSE, draws_format = "matrix"), replicated)
  set.seed(42)
  predictor_only = posterior_predict.mcpfit(
    fit,
    newdata = data.frame(x = data$x), arma = "replicate"
  )

  expect_equal(conditional$.prediction, rep(c(0, 5, 10), each = 2), tolerance = 0.01)
  expect_equal(unname(replicated), matrix(0, nrow = 2, ncol = 3), tolerance = 0.01)
  expect_equal(predictor_only, replicated)

  # Replication generates new responses, so it applies only to predictions
  for (method in list(fitted, residuals, log_lik))
    expect_error(method(fit, arma = "replicate"), "applies only to predictions", fixed = TRUE)
  expect_error(predict(fit, arma = "yes"), "must be TRUE, FALSE, or \"replicate\"", fixed = TRUE)

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
  expect_equal(posterior_predict.mcpfit(fit, arma = "replicate"), replicated)
})


test_that("prior predict uses the same explicit replication setting", {
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
  fit$mcmc_prior = coda::mcmc.list(coda::mcmc(fixed_draws))

  conditional = predict(fit, prior = TRUE, summary = FALSE)
  expect_equal(conditional$.prediction, rep(c(0, 5, 10), each = 2), tolerance = 0.01)
  set.seed(42)
  predicted = predict(fit, prior = TRUE, arma = "replicate", summary = FALSE, width = FALSE)
  set.seed(42)
  predictor_only = predict(
    fit,
    newdata = data.frame(x = data$x),
    prior = TRUE,
    arma = "replicate",
    summary = FALSE,
    width = FALSE
  )

  expect_equal(predicted$.prediction, rep(0, 6), tolerance = 0.01)
  expect_equal(predictor_only$.prediction, predicted$.prediction)
})


test_that("generated JAGS uses the same bounded GARMA residuals", {
  count_data = data.frame(x = 1:4, y = c(0, 1, 4, 2))
  binomial_data = transform(count_data, N = 4)
  poisson_fit = mcp(
    list(y ~ 1 + ar(1, threshold = 0.2) + ma(1)),
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
      ~ 1 + ar(1, threshold = 0.2)
    ),
    count_data,
    family = poisson(),
    par_x = "x",
    sample = FALSE
  )

  expect_match(poisson_fit$jags_code, "garma_y_\\[i_\\] = max\\(y\\[i_\\], garma_threshold_\\[i_\\]\\)")
  expect_match(poisson_fit$jags_code, "garma_link_y_\\[i_\\] = log\\(garma_y_\\[i_\\]\\)")
  expect_match(poisson_fit$jags_code, "resid_abs_\\[i_\\] = garma_link_y_\\[i_\\] - link_mu_\\[i_\\]")
  expect_match(poisson_fit$jags_code, "resid_ma_\\[i_\\] = garma_link_y_\\[i_\\] - link_mu_\\[i_\\] - resid_garma_\\[i_\\]")
  expect_match(poisson_fit$jags_code, "ma1_\\[i_\\] \\* resid_ma_\\[i_ - 1\\]")
  expect_match(binomial_fit$jags_code, "garma_y_\\[i_\\] = min\\(max\\(y\\[i_\\], garma_threshold_\\[i_\\]\\), N\\[i_\\] - garma_threshold_\\[i_\\]\\) / N\\[i_\\]")
  expect_match(binomial_fit$jags_code, "garma_link_y_\\[i_\\] = logit\\(garma_y_\\[i_\\]\\)")
  expect_match(segmented_fit$jags_code, "garma_threshold_\\[i_\\] = 0\\.2\n")
})


test_that("missing GARMA responses stay paired with posterior draws", {
  data = data.frame(x = 1:5, y = c(1, NA, 2, 3, 4))
  expect_message(
    fit <- suppressWarnings(mcp(
      list(y ~ 1 + ar(1)), data, par_x = "x",
      chains = 1, iter = 30, warmup = 20, quiet = TRUE, seed = 42
    )),
    "NA values detected in 'y'"
  )

  expect_equal(fit$.internal$imputed_response_rows, 2L)
  expect_equal(ncol(fit$.internal$imputed_response[[1]]), 1L)
  expect_false("y[2]" %in% colnames(.subset2(fit, "mcmc_post")[[1]]))

  fitted_draws = fitted(fit, summary = FALSE, width = FALSE)
  fitted_again = fitted(fit, summary = FALSE, width = FALSE)
  prediction_draws = predict(fit, summary = FALSE, width = FALSE)
  imputed = as.vector(t(get_imputed_response(fit, fit$data$y, unique(prediction_draws$.draw))))  # Tidy order

  expect_equal(fitted_draws$.epred, fitted_again$.epred)
  expect_true(all(is.na(fitted_draws$y[fitted_draws$.row == 2])))
  expect_false(isTRUE(all.equal(
    prediction_draws$.prediction[prediction_draws$.row == 2],
    imputed[prediction_draws$.row == 2]
  )))

  row3 = fitted_draws$.row == 3
  imputed_by_draw = imputed[prediction_draws$.row == 2][match(
    fitted_draws$.draw[row3],
    prediction_draws$.draw[prediction_draws$.row == 2]
  )]
  expected_row3 = fitted_draws$Intercept_1[row3] +
    fitted_draws$ar1_1[row3] * (imputed_by_draw - fitted_draws$Intercept_1[row3])
  expect_equal(fitted_draws$.epred[row3], expected_row3)

  fitted_summary = fitted(fit)
  prediction_summary = predict(fit)
  expect_false(is.na(fitted_summary$.mean[2]))
  expect_false(is.na(prediction_summary$.mean[2]))
  expect_true(is.na(residuals(fit)$.mean[2]))

  expect_error(log_lik(fit), "missing response occurs in the history")
  expect_error(loo(fit), "missing response occurs in the history")
  expect_error(waic(fit), "missing response occurs in the history")
  expect_equal(ncol(log_lik(fit, arma = FALSE)), 4L)
})


test_that("one GARMA recurrence handles generated and supplied histories", {
  # Unequal lengths exercise series starts and unavailable higher-order lags.
  # Row-varying coefficients and segment transitions represent segment/group changes.
  # Columns are draws with different coefficients.
  n = 12
  cells = function(value) matrix(value, n, 2)
  series = rep(c("01", "1", "site,a", "a[b]"), c(1, 4, 2, 5))
  base = cells(seq(-0.2, 0.2, length.out = n))
  ar = list(ar1 = cbind(seq(0.1, 0.4, length.out = n), 0.3), ar2 = cells(0.1))
  ma = list(ma1 = cells(-0.2), ma2 = cbind(seq(0, 0.1, length.out = n), -0.1))
  threshold = 0.2
  for (family in list(gaussian(), poisson(), binomial(), negbinomial())) {
    family = as_mcpfamily(family)
    data = if (family$family == "binomial") list(trials = rep(c(3, 7), 6)) else list()
    args = list(link_mu = base, ar = ar, ma = ma,
                threshold = threshold, family = family, series = series,
                dpars = list(mu = family$linkinv(base), sigma = cells(1), shape = cells(2)),
                data = data)
    set.seed(21)
    generated = suppressMessages(do.call(simulate_garma, args))
    supplied = do.call(simulate_garma, c(args, list(y = generated$y)))
    expect_equal(supplied, generated)

    # Independent row-by-row reference for each draw: coefficients belong to the
    # current row; history crosses segment changes but never independent series.
    link_y = family$linkfun(family$garma$observed_r(generated$y, data, threshold))
    expected = matrix(0, n, 2)
    for (draw in 1:2) {
      for (rows in split(seq_len(n), series)) {
        for (position in seq_along(rows)) {
          row = rows[position]
          previous = rev(head(rows, position - 1))
          previous = head(previous, 2)
          phi = c(ar$ar1[row, draw], ar$ar2[row, draw])[seq_along(previous)]
          theta = c(ma$ma1[row, draw], ma$ma2[row, draw])[seq_along(previous)]
          expected[row, draw] = base[row, draw] + sum(phi * (link_y[previous, draw] - base[previous, draw])) +
            sum(theta * (link_y[previous, draw] - expected[previous, draw]))
        }
      }
    }
    expect_equal(generated$link_mu, expected)
  }
})


test_that("GARMA thresholds enforce family-appropriate limits during model validation", {
  data_bern = data.frame(x = 1:6, y = c(0, 1, 0, 1, 0, 1))

  # Bernoulli: threshold must be strictly below 0.5
  expect_error(
    mcp(list(y ~ 1 + ar(1, threshold = 0.75)), data_bern, family = bernoulli(), par_x = "x", sample = FALSE),
    "`threshold` for family = bernoulli() must be strictly between 0 and 0.5.",
    fixed = TRUE
  )
  expect_error(
    mcp(list(y ~ 1 + ar(1, threshold = 0.5)), data_bern, family = bernoulli(), par_x = "x", sample = FALSE),
    "`threshold` for family = bernoulli() must be strictly between 0 and 0.5.",
    fixed = TRUE
  )
  expect_silent(
    mcp(list(y ~ 1 + ar(1, threshold = 0.1)), data_bern, family = bernoulli(), par_x = "x", sample = FALSE)
  )
  expect_silent(
    mcp(list(y ~ 1 + ar(1, threshold = 0.49)), data_bern, family = bernoulli(), par_x = "x", sample = FALSE)
  )

  # Binomial: threshold must be strictly below half the smallest trial count
  data_binom_single = data.frame(x = 1:6, y = c(0, 1, 0, 1, 0, 1), N = 1)
  expect_error(
    mcp(list(y | trials(N) ~ 1 + ar(1, threshold = 0.75)), data_binom_single, family = binomial(), par_x = "x", sample = FALSE),
    "`threshold` for family = binomial() must be strictly between 0 and half the smallest trial count (0.5).",
    fixed = TRUE
  )
  expect_error(
    mcp(list(y | trials(N) ~ 1 + ar(1, threshold = 0.5)), data_binom_single, family = binomial(), par_x = "x", sample = FALSE),
    "`threshold` for family = binomial() must be strictly between 0 and half the smallest trial count (0.5).",
    fixed = TRUE
  )

  # Binomial with varying trials: smallest trial count governs the limit
  data_binom_var = data.frame(x = 1:6, y = c(0, 1, 0, 1, 0, 1), N = c(1, 10, 10, 10, 10, 10))
  expect_error(
    mcp(list(y | trials(N) ~ 1 + ar(1, threshold = 0.75)), data_binom_var, family = binomial(), par_x = "x", sample = FALSE),
    "`threshold` for family = binomial() must be strictly between 0 and half the smallest trial count (0.5).",
    fixed = TRUE
  )

  # Binomial with larger trials allows threshold = 0.75 (0.75 < 10 / 2)
  data_binom_large = data.frame(x = 1:6, y = c(0, 1, 0, 1, 0, 1), N = 10)
  expect_silent(
    mcp(list(y | trials(N) ~ 1 + ar(1, threshold = 0.75)), data_binom_large, family = binomial(), par_x = "x", sample = FALSE)
  )

  # simulate_garma directly enforces thresholds
  expect_error(
    simulate_garma(
      link_mu = matrix(0, 4, 1),
      ar = list(ar1 = matrix(0.5, 4, 1)),
      ma = list(),
      threshold = 0.75,
      family = bernoulli(),
      dpars = list(mu = matrix(0.5, 4, 1))
    ),
    "`threshold` for family = bernoulli() must be strictly between 0 and 0.5.",
    fixed = TRUE
  )
})
