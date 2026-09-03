expect_lifetimes = function(x, expected) {
  if (inherits(x, "mcpfit"))
    x = get_fit_model_tables(x)$predictors
  actual = stats::setNames(x$next_segment, x$code_name)[names(expected)]
  testthat::expect_equal(unname(actual), unname(expected))
}
