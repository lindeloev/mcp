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
      seed = 42, diagnostics = FALSE
    )),
    "NA values detected in 'y'"
  )

  fitted_draws = fitted(fit, summary = FALSE, width = FALSE)
  prediction_draws = predict(fit, summary = FALSE, width = FALSE)
  missing_row = prediction_draws$.row == 10
  completed_values = as.vector(t(get_imputed_response(fit, fit$data$y, unique(prediction_draws$.draw))))  # Tidy order

  expect_equal(fit$.internal$imputed_response_rows, 10L)
  expect_false(isTRUE(all.equal(prediction_draws$.prediction[missing_row], completed_values[missing_row])))
  imputations = imputed(fit, summary = FALSE)
  expect_named(imputations, c(names(fit$data), ".row", ".chain", ".iteration", ".draw", ".imputed"))
  expect_equal(imputations$.imputed, completed_values[missing_row])
  expect_equal(unique(imputations$.row), 10L)
  expect_equal(nrow(imputed(fit, summary = FALSE, ndraws = 2)), 2L)
  expect_error(imputed(list()), "Must inherit from class 'mcpfit'")
  expect_error(imputed(fit, ndraws = 0), "Element 1 is not >= 1")
  expect_true(all(is.na(prediction_draws$y[missing_row])))
  expect_true(all(prediction_draws$id[missing_row] == "b"))
  expect_gt(
    mean(fitted_draws$.epred[fitted_draws$.row == 10]),
    mean(data$y, na.rm = TRUE) + 50
  )
  expect_gt(
    mean(prediction_draws$.prediction[missing_row]),
    mean(data$y, na.rm = TRUE) + 50
  )
  prediction_summary = predict(fit)
  expect_identical(prediction_summary$.row, seq_len(nrow(data)))
  expect_identical(fitted(fit)$.row, prediction_summary$.row)
  expect_true(all(is.finite(unlist(prediction_summary[10, c(".lower", ".upper")]))))
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
  imputations = imputed(fit, summary = FALSE)
  expect_equal(unname(epred[, 2]), rep(0, nrow(epred)))
  expect_equal(mean(predicted[, 2]), 0, tolerance = 0.06)
  expect_equal(stats::sd(predicted[, 2]), 1, tolerance = 0.06)
  expect_equal(mean(imputations$.imputed), 8 / 1.64, tolerance = 0.06)
  expect_equal(stats::sd(imputations$.imputed), sqrt(1 / 1.64), tolerance = 0.06)
  expect_equal(unname(epred[, 3]), 0.8 * imputations$.imputed)
  expect_equal(colMeans(predicted), colMeans(epred), tolerance = 0.06)

  # Mixture intervals use the same conditional distribution at missing rows.
  pred_sum = predict(fit)
  expect_equal(unname(unlist(pred_sum[2, c(".lower", ".upper")])), qnorm(c(0.025, 0.975)), tolerance = 1e-5)
  expect_identical(log_lik(fit, arma = FALSE, summary = TRUE)$.row, c(1L, 3L))
  expect_identical(colnames(log_lik(fit, arma = FALSE)), c("1", "3"))
  expect_equal(nrow(imputed(demo_fit)), 0L)

  fit$data$.imputed = "user column"
  expect_error(imputed(fit, summary = FALSE), "data contain a column named `.imputed`", fixed = TRUE)
  expect_identical(imputed(fit)$.imputed, "user column")
  expect_identical(colnames(imputed(fit, summary = FALSE, draws_format = "matrix")), "2")
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
      chains = 2, iter = 30, warmup = 20, seed = 42, quiet = TRUE, diagnostics = FALSE
    )),
    "NA values detected"
  )
  expected = posterior_epred.mcpfit(fit)
  imputations = imputed(fit, summary = FALSE)
  raw = tibble::as_tibble(as_draws_df(fit))
  parameter_names = setdiff(names(raw), c(".chain", ".iteration", ".draw"))
  for (draw in c(1, nrow(raw))) {
    completed = data
    use = imputations$.draw == raw$.draw[draw]
    completed$y[imputations$.row[use]] = imputations$.imputed[use]
    single = fit
    single$mcmc_post = coda::mcmc.list(coda::mcmc(as.matrix(raw[draw, parameter_names])))
    expect_equal(unname(posterior_epred.mcpfit(single, newdata = completed)[1, ]), unname(expected[draw, ]))
  }
  set.seed(42)
  selected = fitted(fit, summary = FALSE, ndraws = 2)
  expect_equal(selected$.epred, unname(expected[cbind(match(selected$.draw, raw$.draw), selected$.row)]))
  expect_identical(imputations$.row, rep(which(is.na(data$y)), each = nrow(raw)))
  completed_draws = imputed(fit, complete = TRUE, summary = FALSE)
  completed_matrix = imputed(fit, complete = TRUE, summary = FALSE, draws_format = "matrix")
  expect_equal(as.vector(completed_matrix), completed_draws$.imputed)
  expect_equal(completed_draws$.imputed[is.na(completed_draws$y)], imputations$.imputed)
  expect_equal(completed_draws$.imputed[!is.na(completed_draws$y)], completed_draws$y[!is.na(completed_draws$y)])
  one_draw = imputed(fit, complete = TRUE, ndraws = 1, width = FALSE)
  expect_equal(one_draw$.mean[!is.na(data$y)], data$y[!is.na(data$y)])
  expect_true(all(one_draw$.sd[!is.na(data$y)] == 0))
  expect_true(all(is.na(one_draw$.sd[is.na(data$y)])))
  expect_error(fitted(fit, group_effects = FALSE), "requires `group_effects = TRUE`", fixed = TRUE)
  fit$mcmc_prior = .subset2(fit, "mcmc_post")
  expect_error(predict(fit, prior = TRUE), "Missing GARMA histories require posterior draws", fixed = TRUE)
  expect_no_error(predict(fit, prior = TRUE, arma = "replicate", summary = FALSE, ndraws = 2))
  expect_no_error(ggplot2::ggplot_build(plot(fit, facet_by = "id")))
})
