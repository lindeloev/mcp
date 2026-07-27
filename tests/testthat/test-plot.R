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
