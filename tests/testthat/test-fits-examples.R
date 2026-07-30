# COMMENT THIS LINE TO ENABLE EXTENSIVE TESTS. THEN REVERT. DO NOT COMMIT.
# It's a slow test so we do it rarely.
testthat::skip("This time-consuming test is only run locally before release.")
testthat::skip_if_not_installed("vdiffr")

# Runs through all examples in mcp::mcp_example() and verifies the fits and plots.
# Set MCP_MAKE_PLOT_TEST_FILE to also write all plots to a multi-page PDF for review.
test_fits_examples = function() {
  examples = c("ar", "binomial", "demo", "group", "intercepts", "multiple", "quadratic", "variance", "varying")
  snapshot_examples = c("demo", "group", "variance")
  plot_file = Sys.getenv("MCP_MAKE_PLOT_TEST_FILE")

  grDevices::pdf(
    file = if (nzchar(plot_file)) path.expand(plot_file) else NULL,
    width = 10,
    height = 8,
    onefile = TRUE
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  for (example in examples) {
    message("Now running example mcp_example('", example, "')")
    fit = mcp_example(example, plot = TRUE)
    example_plot = ggplot2::last_plot()

    test_matches_simulated(fit)
    test_example_plot(example_plot, example)
    if (example %in% snapshot_examples)
      snapshot_example_plot(example_plot, example)
  }
}

testthat::test_that("mcp examples recover their parameters and plots render", {
  test_fits_examples()
})
