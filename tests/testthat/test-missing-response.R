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
  expect_equal(prediction_draws$.prediction[missing_row], imputed[missing_row])
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
  expect_no_warning(ggplot2::ggplot_build(plot(fit, lines = 0, cp_dens = FALSE)))
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
