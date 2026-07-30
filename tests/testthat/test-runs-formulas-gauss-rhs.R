#################
# TEST RESPONSE #
#################
bad_y = list(
  list( ~ 1),  # No y
  list((1|id) ~ 1),  # y cannot be varying
  list(1 ~ 1),  # 1 is not y
  list(y ~ 1,  # Two y
       a ~ 1 ~ 1),
  list(y ~ 1,  # Intercept y
       1 ~ 1 ~ 1),
  list(bad_y_char ~ 1),  # Character y
  list(bad_y_factor ~ 1)  # Factor y
)

test_bad(bad_y)


good_y = list(
  list(y ~ 1),  # Regular
  list(y ~ 1,  # Explicit and implicit y and cp
       y ~ 1 ~ 1,
       1 + (1|id) ~ 1 + x,
       ~ 1),
  list(ok_y ~ 1)  # decimal y
)

test_good(good_y)



###################
# TEST INTERCEPTS #
###################
bad_intercepts = list(
  list(y ~ 2)  # 2 not supported
)

test_bad(bad_intercepts)


good_intercepts = list(
  #list(y ~ 0),  # would be nice if it worked, but mcmc.list does not behave well with just one variable
  list(ok_y ~ 1),  # y can be called whatever
  list(y ~ 0,  # Multiple segments
       ~ 1,
       ~ 0,
       ~ 1),
  list(y ~ 1 + (1 | id),  # Predictor group intercept carries until explicit turn-off
       ~ 1,
       ~ 1 + (0 | id))
)

test_good(good_intercepts)


###############
# TEST SLOPES #
###############
bad_slopes = list(
  list(y ~ 1,
       ~ log(x)),  # should fail explicitly because negative x
  list(y ~ 1,
       ~ sqrt(x))  # should fail explicitly because negative x
)

test_bad(bad_slopes)

test_that("formula functions reject multiple terms containing par_x", {
  expect_no_warning(
    expect_error(
      mcp(list(y ~ I(x:ok_x)), data_gauss, par_x = "x", sample = FALSE),
      "does not currently support 2\\+ terms within a formula function"
    )
  )
})

test_that("transformations use the original par_x while segment bases stay local", {
  data = data.frame(x = 0:10, y = 0)
  fit = mcp(
    list(y ~ 1, ~ 1 + sin(x) + x + I(x^2)),
    data,
    sample = FALSE
  )
  table = get_fit_model_tables(fit)$predictors

  expect_equal(table$x_factor[table$matrix_name == "sin(x)"], "1")
  expect_equal(
    unname(table$matrix_data[[which(table$matrix_name == "sin(x)")]]),
    sin(data$x)
  )
  expect_equal(table$x_factor[table$matrix_name == "x"], "x")
  expect_equal(table$x_factor[table$matrix_name == "I(x^2)"], "x^2")

  fitted = fit$simulate(
    fit,
    data,
    cp_1 = 5,
    Intercept_1 = 0,
    Intercept_2 = 0,
    sinx_2 = 2,
    x_2 = 0,
    xE2_2 = 0,
    sigma_1 = 1,
    .type = "fitted"
  )
  expected = ifelse(data$x < 5, 0, 2 * sin(data$x))
  expect_equal(as.numeric(fitted), expected)
})



good_slopes = list(
  list(y ~ 0 + x),  # Regular
  list(y ~ 0 + x,  # Multiple on/off
       ~ 0,
       ~ 1 + x),
  list(y ~ 0 + x + I(x^2) + I(x^3),  # Test "non-linear" x
       ~ 0 + exp(x) + abs(x),
       ~ 0 + sin(x) + cos(x) + tan(x)),
  list(y ~ ok_x)  # alternative x
)

test_good(good_slopes, par_x = NULL)
