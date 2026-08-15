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
