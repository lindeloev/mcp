expect_lifetimes = function(x, expected) {
  if (inherits(x, "mcpfit"))
    x = get_fit_model_tables(x)$coefs
  actual = stats::setNames(x$next_segment, x$name)[names(expected)]
  testthat::expect_equal(unname(actual), unname(expected))
}


# The population coefficient table of a model
get_coefs = function(model, data, family, par_x) {
  get_coef_tables(model, data, family, par_x)$coefs
}
