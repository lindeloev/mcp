test_that("plot_pars() paginates without changing its single-page return", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  plots = lapply(c(5, 7), function(nvariables) {
    plot_pars(
      demo_fit, type = "dens_overlay",
      nvariables = nvariables, ask = FALSE
    )
  })
  pages = plots[[1]]
  page_sizes = vapply(
    pages,
    function(page) length(unique(page[[1]]$data$Parameter)),
    integer(1)
  )

  expect_equal(unname(page_sizes), c(5L, 2L))
  expect_true(all(vapply(pages, inherits, logical(1), what = "ggplot")))
  expect_s3_class(plots[[2]], "ggplot")

  # nvariables = NULL and nvariables = Inf plot all parameters on one page
  p_null = plot_pars(demo_fit, type = "dens_overlay", nvariables = NULL)
  p_inf = plot_pars(demo_fit, type = "dens_overlay", nvariables = Inf)
  expect_s3_class(p_null, "ggplot")
  expect_s3_class(p_inf, "ggplot")
  expect_equal(length(unique(p_null[[1]]$data$Parameter)), 7L)
  expect_equal(length(unique(p_inf[[1]]$data$Parameter)), 7L)
})

test_that("plot_pars() prefers the group selector and supports its deprecated alias", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  fit = unclass(demo_fit)
  fit$.internal$model_tables$parameters = dplyr::bind_rows(
    fit$.internal$model_tables$parameters,
    tibble::tibble(
      name = "group_par", part = "predictor", scope = "group",
      role = "group_deviation", segment = 1L, dpar = "mu", order = NA_integer_,
      group_col = "group", population_name = "Intercept_1"
    )
  )
  for (chain in seq_along(fit$mcmc_post))
    colnames(fit$mcmc_post[[chain]])[1] = "group_par[1]"
  class(fit) = class(demo_fit)

  expect_no_warning(
    group_plot <- plot_pars(fit, pars = "group", type = "dens", ask = FALSE)
  )
  expect_warning(
    varying_plot <- plot_pars(fit, pars = "varying", type = "dens", ask = FALSE),
    "deprecated"
  )
  expect_s3_class(group_plot, "ggplot")
  expect_s3_class(varying_plot, "ggplot")
})

test_that("geom_cp_density draws a vertical spike for fixed change points", {
  # Unclass because accessing fit$mcmc_* elicits a deprecation warning
  fit = unclass(demo_fit)
  for (chain in seq_along(fit$mcmc_post)) {
    fit$mcmc_post[[chain]][, "cp_1"] = 50
  }
  class(fit) = class(demo_fit)
  p = plot(fit)
  gb = ggplot2::ggplot_build(p)
  cp_data = gb$data[[length(gb$data)]]

  # cp_1 is fixed at 50, so its x coordinates should be strictly 50
  cp_1_data = cp_data[cp_data$group %in% c(1, 2), ]
  expect_true(all(cp_1_data$x == 50))

  # Peak height should match varying change points (e.g. cp_2)
  cp_2_data = cp_data[!cp_data$group %in% c(1, 2), ]
  expect_equal(max(cp_1_data$y), max(cp_2_data$y))
})


test_that("plot() supports exact q_predict, deprecates samples, and checks dots", {
  # Exact q_predict quantiles
  p = plot(demo_fit, q_fit = TRUE, q_predict = TRUE, lines = 0)
  expect_s3_class(p, "ggplot")

  # Without lines and using all draws, plot is completely deterministic across consecutive calls
  p1 = plot(demo_fit, q_fit = TRUE, q_predict = TRUE, lines = 0, ndraws = NULL)
  p2 = plot(demo_fit, q_fit = TRUE, q_predict = TRUE, lines = 0, ndraws = NULL)
  expect_equal(p1$layers[[2]]$data, p2$layers[[2]]$data)
  expect_equal(p1$layers[[3]]$data, p2$layers[[3]]$data)

  # Deprecation warning for samples and nsamples
  expect_warning(
    plot(demo_fit, samples = 10),
    "deprecated"
  )
  expect_warning(
    plot(demo_fit, nsamples = 10),
    "deprecated"
  )

  # scale = "linear" in plot.mcpfit
  p_lin = plot(demo_fit, scale = "linear")
  expect_s3_class(p_lin, "ggplot")

  # Unknown dots trigger an error
  expect_error(
    plot(demo_fit, weird_arg = 123),
    "must be empty"
  )
  expect_error(
    plot_dpar(demo_fit, weird_arg = 123),
    "must be empty"
  )
})

