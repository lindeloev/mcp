test_that("missing response draws follow covariates and group-level effects", {
  data = data.frame(
    x = c(0:4, 0:3, 20),
    id = factor(c(rep("a", 5), rep("b", 5))),
    y = c(2 + 3 * (0:4), 52 + 3 * (0:3), NA)
  )
  expect_message(
    fit <- suppressWarnings(mcp(
      list(y ~ 1 + x + (1 | id)), data, par_x = "x",
      chains = 1, iter = 40, warmup = 20, quiet = TRUE,
      seed = 8, diagnostics = FALSE
    )),
    "NA values detected in 'y'"
  )

  fitted_draws = fitted(fit, summary = FALSE, probs = FALSE)
  prediction_draws = predict(fit, summary = FALSE, probs = FALSE)
  missing_row = prediction_draws$data_row == 10
  imputed = get_imputed_response_draws(fit, prediction_draws)

  expect_equal(fit$.internal$imputed_response_rows, 10L)
  expect_false(isTRUE(all.equal(prediction_draws$.prediction[missing_row], imputed[missing_row])))
  retained = imputed_draws(fit)
  expect_equal(retained$.imputed, imputed[missing_row])
  expect_equal(unique(retained$data_row), 10L)
  expect_equal(nrow(imputed_draws(fit, ndraws = 2)), 2L)
  expect_true(all(is.na(prediction_draws$y[missing_row])))
  expect_true(all(prediction_draws$id[missing_row] == "b"))
  expect_gt(
    mean(fitted_draws$.epred[fitted_draws$data_row == 10]),
    mean(data$y, na.rm = TRUE) + 50
  )
  expect_gt(
    mean(prediction_draws$.prediction[missing_row]),
    mean(data$y, na.rm = TRUE) + 50
  )
  prediction_summary = predict(fit)
  expect_true(all(is.finite(unlist(prediction_summary[10, c("Q2.5", "Q97.5")]))))
  expect_no_warning(ggplot2::ggplot_build(plot(fit, lines = 0, cp_dens = FALSE)))
})


test_that("missing-response predictions target the conditional response distribution", {
  data = data.frame(x = 1:3, y = c(0, NA, 10))
  expect_message(
    fit <- suppressWarnings(mcp(
      list(y ~ 1 + ar(1)), data, par_x = "x",
      prior = list(Intercept_1 = 0, ar1_1 = 0.8, sigma_1 = 1),
      chains = 2, iter = 4000, warmup = 100, quiet = TRUE, seed = 42
    )),
    "NA values detected in 'y'"
  )
  epred = posterior_epred.mcpfit(fit)
  predicted = posterior_predict.mcpfit(fit, seed = 17)
  retained = imputed_draws(fit)
  expect_equal(unname(epred[, 2]), rep(0, nrow(epred)))
  expect_equal(mean(predicted[, 2]), 0, tolerance = 0.06)
  expect_equal(stats::sd(predicted[, 2]), 1, tolerance = 0.06)
  expect_equal(mean(retained$.imputed), 8 / 1.64, tolerance = 0.06)
  expect_equal(stats::sd(retained$.imputed), sqrt(1 / 1.64), tolerance = 0.06)
  expect_equal(unname(epred[, 3]), 0.8 * retained$.imputed)
  expect_equal(colMeans(predicted), colMeans(epred), tolerance = 0.06)

  # Mixture intervals use the same conditional distribution at missing rows.
  pred_sum = predict(fit)
  expect_equal(unname(unlist(pred_sum[2, c("Q2.5", "Q97.5")])), qnorm(c(0.025, 0.975)), tolerance = 1e-5)
  expect_error(imputed_draws(demo_fit), "No retained posterior imputations")
})


test_that("change point densities are shown with missing responses", {
  fit = demo_fit
  fit$data$y[c(1, 5, 10)] = NA

  plotted = plot(fit, lines = 0)
  density_layer = which(vapply(
    plotted$layers,
    function(layer) inherits(layer$geom, "GeomPolygon"),
    logical(1)
  ))
  density = ggplot2::ggplot_build(plotted)$data[[density_layer]]

  expect_length(density_layer, 1)
  expect_gt(nrow(density), 0)
  expect_true(all(is.finite(density$y)))
})


test_that("completed histories stay paired across draws, groups, series, and segments", {
  data = data.frame(x = rep(1:6, 2), id = rep(c("01", "site,a"), each = 6),
                    y = c(1, NA, 2, 3, 4, 3, 4, 5, NA, 6, 5, 7))
  expect_message(
    fit <- suppressWarnings(mcp(
      list(y ~ 1 + (1 | id) + ar(2, series = id) + ma(1), ~ 1 + ar(1, series = id) + ma(2)),
      data, par_x = "x", prior = list(cp_1 = 3.5),
      chains = 2, iter = 30, warmup = 20, seed = 51, quiet = TRUE, diagnostics = FALSE
    )),
    "NA values detected"
  )
  expected = posterior_epred.mcpfit(fit)
  imputed = imputed_draws(fit)
  raw = tibble::as_tibble(as_draws_df(fit))
  parameter_names = setdiff(names(raw), c(".chain", ".iteration", ".draw"))
  for (draw in c(1, nrow(raw))) {
    completed = data
    use = imputed$.draw == raw$.draw[draw]
    completed$y[imputed$data_row[use]] = imputed$.imputed[use]
    single = fit
    single$mcmc_post = coda::mcmc.list(coda::mcmc(as.matrix(raw[draw, parameter_names])))
    expect_equal(unname(posterior_epred.mcpfit(single, newdata = completed)[1, ]), unname(expected[draw, ]))
  }
  set.seed(42)
  selected = fitted(fit, summary = FALSE, ndraws = 2)
  expect_equal(selected$.epred, unname(expected[cbind(match(selected$.draw, raw$.draw), selected$data_row)]))
  expect_error(fitted(fit, varying = FALSE), "requires `varying = TRUE`", fixed = TRUE)
  fit$mcmc_prior = .subset2(fit, "mcmc_post")
  expect_error(predict(fit, prior = TRUE), "Missing GARMA histories require posterior draws", fixed = TRUE)
  expect_no_error(predict(fit, prior = TRUE, conditional = FALSE, summary = FALSE, ndraws = 2))
})
