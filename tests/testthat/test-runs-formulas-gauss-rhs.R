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
  list(y ~ 0),
  list(ok_y ~ 1),  # y can be called whatever
  list(y ~ 0,  # Multiple segments
       ~ 1,
       ~ 0,
       ~ 1),
  list(y ~ 1 + (1 | id),  # Predictor group intercept carries until explicit turn-off
       ~ 1,
       ~ 1 + (0 | id)),
  list(y ~ 1 + (ok_id_factor || id))  # Independent intercept and factor deviations
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

test_that("formula offsets work in mcp formulas", {
  d = data.frame(
    x = 1:20,
    z = rep(c(2, 5), each = 10),
    pop = rep(c(10, 100), each = 10),
    y = 5 + 0.5 * 1:20 + rep(c(2, 5), each = 10)
  )

  # Gaussian with offset
  fit_gauss = mcp(
    list(y ~ 1 + x + offset(z), ~ 1 + x + offset(z)),
    d, sample = FALSE
  )
  expect_true(grepl("offset_mu_1_\\[i_\\]", fit_gauss$jags_code))
  expect_true(grepl("offset_mu_2_\\[i_\\]", fit_gauss$jags_code))
  expect_false("z" %in% mcp_pars(fit_gauss)$name)
  expect_false("offset_mu_1_" %in% mcp_pars(fit_gauss)$name)

  # Poisson with offset(log(pop))
  d_pois = data.frame(
    x = 1:20,
    pop = rep(c(10, 100), each = 10),
    y = rpois(20, rep(c(10, 100), each = 10) * exp(0.05 * 1:20))
  )
  fit_pois = mcp(
    list(y ~ 1 + x + offset(log(pop)), ~ 1 + x + offset(log(pop))),
    d_pois, family = poisson(), sample = FALSE
  )
  expect_true(grepl("offset_mu_1_\\[i_\\]", fit_pois$jags_code))
  expect_true(grepl("offset_mu_2_\\[i_\\]", fit_pois$jags_code))

  # fit$simulate with offset
  sim_out = fit_gauss$simulate(
    fit_gauss,
    data.frame(x = 1:4, z = c(10, 10, 20, 20)),
    cp_1 = 2.5, Intercept_1 = 5, x_1 = 1, Intercept_2 = 10, x_2 = 2, sigma_1 = 0.001
  )
  # Segment 1 (x=1, z=10): mu = 5 + 1*1 + 10 = 16
  # Segment 1 (x=2, z=10): mu = 5 + 1*2 + 10 = 17
  # Segment 2 (x=3, z=20): mu = 10 + 2*(3-2.5) + 20 = 31
  # Segment 2 (x=4, z=20): mu = 10 + 2*(4-2.5) + 20 = 33
  expect_equal(unname(attr(sim_out, "dpars")$mu), c(16, 17, 31, 33))

  # Offset in sigma
  fit_sigma = mcp(list(y ~ 1 + x + sigma(1 + offset(z))), d, sample = FALSE)
  expect_true(grepl("offset_sigma_1_\\[i_\\]", fit_sigma$jags_code))

  # Warning when par_x is in offset()
  expect_warning(
    mcp(list(y ~ 1 + offset(log(x)), ~ 1 + x), d_pois, family = poisson(), sample = FALSE),
    regexp = "The change-point predictor 'x' is also used in offset\\(\\)"
  )
  expect_warning(
    mcp(list(y ~ 1 + offset(x), ~ 1), d, par_x = "x", sample = FALSE),
    regexp = "The change-point predictor 'x' is also used in offset\\(\\)"
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

test_that("par_x is recognized as an exact formula symbol", {
  data = data.frame(a.b = 0:4, axb = c(2, 4, 6, 8, 10), y = 0)
  fit = mcp(list(y ~ 0 + axb), data, par_x = "a.b", sample = FALSE)
  predictor = dplyr::filter(get_fit_model_tables(fit)$predictors, .data$dpar == "mu")

  expect_equal(predictor$x_factor, "1")
  expect_equal(unname(predictor$matrix_data[[1]]), data$axb)
})

test_that("x-free interactions are evaluated directly at zero", {
  data = data.frame(x = 0:4, z = c(2, 0, 0, 4, 5), y = 0)
  fit = mcp(list(y ~ 0 + x:z), data, par_x = "x", sample = FALSE)
  predictor = dplyr::filter(get_fit_model_tables(fit)$predictors, .data$dpar == "mu")

  expect_equal(predictor$x_factor, "x")
  expect_equal(unname(predictor$matrix_data[[1]]), data$z)

  segmented = mcp(list(y ~ 0, ~ 0 + x:z), data, par_x = "x", sample = FALSE)
  fitted = segmented$simulate(
    segmented, data, cp_1 = 0, xz_2 = 1, sigma_1 = 1, .type = "fitted"
  )
  expect_equal(as.numeric(fitted), data$x * data$z)
})

test_that("factor punctuation produces safe parameter names", {
  data = data.frame(
    x = 1:4,
    condition = factor(
      c("reference", "a/b", "reference", "a/b"),
      levels = c("reference", "a/b")
    ),
    y = 0
  )
  fit = mcp(list(y ~ condition), data, par_x = "x", sample = FALSE)
  predictor = get_fit_model_tables(fit)$predictors
  slash_row = grepl("a/b", predictor$matrix_name, fixed = TRUE)
  safe_name = predictor$code_name[slash_row]

  expect_true(length(safe_name) == 1 && grepl("^[A-Za-z][A-Za-z0-9_.]*$", safe_name))
  args = list(fit = fit, newdata = data, Intercept_1 = 0, sigma_1 = 1, .type = "fitted")
  args[[safe_name]] = 1
  expect_no_error(do.call(fit$simulate, args))
})

test_that("prediction reuses fitted factor encodings", {
  data = data.frame(
    x = 1:6,
    y = 1:6,
    condition = ordered(rep(c("low", "middle", "high"), 2))
  )
  ordered_fit = mcp(list(y ~ condition), data, par_x = "x", sample = FALSE)
  ordered_matrix = get_predictor_matrix(
    get_fit_model_tables(ordered_fit)$predictors,
    get_fit_model_tables(ordered_fit)$group_effects
  )
  ordered_new = add_rhs_predictors(
    transform(data, condition = as.character(condition)), ordered_fit
  )
  expect_equal(
    unname(as.matrix(ordered_new[, paste0(".pred_", colnames(ordered_matrix))])),
    unname(ordered_matrix)
  )

  custom_data = transform(data, condition = factor(condition))
  contrasts(custom_data$condition) = stats::contr.sum(3)
  custom_fit = suppressWarnings(mcp(list(y ~ condition), custom_data, par_x = "x", sample = FALSE))
  custom_matrix = get_predictor_matrix(
    get_fit_model_tables(custom_fit)$predictors,
    get_fit_model_tables(custom_fit)$group_effects
  )
  custom_new = add_rhs_predictors(
    transform(custom_data, condition = as.character(condition)), custom_fit
  )
  expect_equal(
    unname(as.matrix(custom_new[, paste0(".pred_", colnames(custom_matrix))])),
    unname(custom_matrix)
  )
})

test_that("prediction is independent of later contrast options", {
  old_options = options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old_options), add = TRUE)
  data = data.frame(
    x = 1:6,
    y = 1:6,
    condition = factor(rep(c("a", "b", "c"), 2))
  )
  fit = mcp(list(y ~ condition), data, par_x = "x", sample = FALSE)
  fitted_matrix = get_predictor_matrix(
    get_fit_model_tables(fit)$predictors,
    get_fit_model_tables(fit)$group_effects
  )

  options(contrasts = c("contr.sum", "contr.poly"))
  new_matrix = add_rhs_predictors(data, fit)
  expect_equal(
    unname(as.matrix(new_matrix[, paste0(".pred_", colnames(fitted_matrix))])),
    unname(fitted_matrix)
  )
})

test_that("data-derived bases reuse their fitted specification", {
  data = data.frame(x = 1:8, y = 1:8)
  newdata = data.frame(x = 9:10)

  scale_fit = mcp(list(y ~ scale(x)), data, sample = FALSE)
  scale_matrix = add_rhs_predictors(newdata, scale_fit)
  scale_pars = get_fit_model_tables(scale_fit)$predictors$code_name[
    get_fit_model_tables(scale_fit)$predictors$dpar == "mu"
  ]
  expect_equal(
    scale_matrix[[paste0(".pred_", setdiff(scale_pars, "Intercept_1"))]],
    as.numeric(scale(newdata$x, center = mean(data$x), scale = stats::sd(data$x)))
  )

  poly_basis = stats::poly(data$x, 2)
  poly_fit = mcp(list(y ~ poly(x, 2)), data, sample = FALSE)
  poly_matrix = add_rhs_predictors(newdata, poly_fit)
  poly_pars = get_fit_model_tables(poly_fit)$predictors$code_name[
    get_fit_model_tables(poly_fit)$predictors$dpar == "mu"
  ]
  expect_equal(
    unname(as.matrix(poly_matrix[, paste0(".pred_", setdiff(poly_pars, "Intercept_1"))])),
    unname(stats::predict(poly_basis, newdata$x))
  )

  spline_basis = splines::ns(data$x, df = 3)
  spline_fit = mcp(list(y ~ splines::ns(x, df = 3)), data, sample = FALSE)
  spline_matrix = add_rhs_predictors(newdata, spline_fit)
  spline_pars = get_fit_model_tables(spline_fit)$predictors$code_name[
    get_fit_model_tables(spline_fit)$predictors$dpar == "mu"
  ]
  expect_equal(
    unname(as.matrix(spline_matrix[, paste0(".pred_", setdiff(spline_pars, "Intercept_1"))])),
    matrix(
      as.numeric(stats::predict(spline_basis, newdata$x)),
      nrow = nrow(newdata)
    )
  )

  bs_newdata = data.frame(x = 2:3)
  bs_basis = splines::bs(data$x, df = 3)
  bs_fit = mcp(list(y ~ splines::bs(x, df = 3)), data, sample = FALSE)
  bs_matrix = add_rhs_predictors(bs_newdata, bs_fit)
  bs_pars = get_fit_model_tables(bs_fit)$predictors$code_name[
    get_fit_model_tables(bs_fit)$predictors$dpar == "mu"
  ]
  expect_equal(
    unname(as.matrix(bs_matrix[, paste0(".pred_", setdiff(bs_pars, "Intercept_1"))])),
    matrix(
      as.numeric(stats::predict(bs_basis, bs_newdata$x)),
      nrow = nrow(bs_newdata)
    )
  )
})


test_that("formulas preserve local transformation environments", {
  data = data.frame(
    x = 1:8,
    y = c(1, 2, 1, 2, 5, 6, 5, 6),
    id = rep(c("a", "b"), 4)
  )
  local_fits = function(data) {
    sq_local = function(x) x^2
    list(
      mu = mcp(list(y ~ sq_local(x)), data, par_x = "x", sample = FALSE, quiet = TRUE),
      sigma = mcp(list(y ~ 1 + sigma(sq_local(x))), data, par_x = "x", sample = FALSE, quiet = TRUE),
      group = mcp(list(y ~ 1 + (sq_local(x) || id)), data, par_x = "x", sample = FALSE, quiet = TRUE),
      ar = mcp(list(y ~ 1 + ar(1, sq_local(x))), data, par_x = "x", sample = FALSE, quiet = TRUE)
    )
  }

  fits = local_fits(data)
  expect_true(all(vapply(fits, function(fit) {
    exists("sq_local", environment(fit$model[[1]]), inherits = FALSE)
  }, logical(1))))
  new_matrix = add_rhs_predictors(data.frame(x = 9:10), fits$mu)
  expect_equal(new_matrix$.pred_sq_localx_1, c(81, 100))
})



good_slopes = list(
  list(y ~ 0 + x),  # Regular
  list(y ~ 0 + x,  # Multiple on/off
       ~ 0,
       ~ 1 + x),
  list(y ~ 0 + x + I(x^2) + I(x^3),  # Test "non-linear" x
       ~ 0 + exp(x) + abs(x),
       ~ 0 + sin(x) + cos(x) + tan(x)),
  list(y ~ poly(x, 2)),
  list(y ~ splines::ns(x, df = 3)),
  list(y ~ ok_x)  # alternative x
)

test_good(good_slopes, par_x = NULL)
