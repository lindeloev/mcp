expect_lifetimes = function(x, expected) {
  if (inherits(x, "mcpfit"))
    x = get_fit_model_tables(x)$predictors
  actual = stats::setNames(x$next_segment, x$name)[names(expected)]
  testthat::expect_equal(unname(actual), unname(expected))
}


# The population predictor table of a model
get_predictors = function(model, data, family, par_x) {
  get_predictor_tables(model, data, family, par_x)$predictors
}
