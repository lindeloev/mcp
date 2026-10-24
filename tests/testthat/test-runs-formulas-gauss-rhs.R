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
  list(y ~ 1 + (1 | id),  # Predictor group intercept is included in this segment
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

test_that("population terms apply only where included while x-terms join", {
  data = data.frame(
    x = 1:8,
    z = c(0, 1, 0, 2, 1, 3, 2, 4),
    w = c(1, 1, 2, 1, 3, 2, 4, 3),
    y = 0
  )
  fit = mcp(
    list(y ~ 1 + x + z + w, ~ 0 + x + z),
    data,
    par_x = "x",
    sample = FALSE
  )
  expect_lifetimes(fit, c(z_1 = 2L, w_1 = 2L, x_1 = NA_integer_))
  expect_match(
    fit$jags_code,
    "(?s)x\\[i_\\] < cp_1.*c\\(z_1, w_1\\)",
    perl = TRUE
  )

  fitted = fit$simulate(
    fit,
    data,
    cp_1 = 4.5,
    Intercept_1 = 1, x_1 = 2, z_1 = 10, w_1 = 3, sigma_1 = 1,
    x_2 = -1, z_2 = 20,
    .type = "fitted"
  )
  expected = 1 + 2 * pmin(data$x, 4.5) +
    ifelse(data$x < 4.5, 10 * data$z + 3 * data$w, -(data$x - 4.5) + 20 * data$z)
  expect_equal(as.numeric(fitted), expected)
})

test_that("x-terms end at a disjoined segment even without an earlier intercept", {
  data = data.frame(x = 1:8, y = 0)
  fit = mcp(list(y ~ 0 + x, ~ 0 + x, ~ 1 + sigma(1)), data, par_x = "x", sample = FALSE)
  expect_lifetimes(fit, c(x_1 = 3L, x_2 = 3L))

  fitted = fit$simulate(fit, data, cp_1 = 3.5, cp_2 = 5.5, x_1 = 1, x_2 = 2, Intercept_3 = 10,
    sigma_Intercept_1 = 1, sigma_Intercept_3 = 1, .type = "fitted")
  expected = ifelse(data$x < 5.5, pmin(data$x, 3.5) + 2 * pmax(data$x - 3.5, 0), 10)
  expect_equal(as.numeric(fitted), expected)
})

test_that("redefinitions replace complete formula terms", {
  data = data.frame(
    x = 1:12,
    z = seq(-1, 1, length.out = 12),
    group = factor(rep(c("a", "b", "c"), 4)),
    y = 0
  )
  predictors = get_predictors(
    list(
      y ~ 1 + group + poly(z, 2),
      ~ 0 + x + group + poly(z, 2)
    ),
    data,
    as_mcpfamily(gaussian()),
    par_x = "x",
    check_rank = FALSE
  )
  first_terms = predictors$segment == 1 & predictors$par_type != "Intercept" & predictors$x_degree == 0

  expect_true(all(predictors$next_segment[first_terms] == 2L))
  expect_equal(
    length(unique(predictors$term_key[grepl("poly", predictors$matrix_name)])),
    1L
  )

  reset = get_predictors(
    list(y ~ 1 + z, ~ 1, ~ 0 + z),
    data,
    as_mcpfamily(gaussian()),
    par_x = "x"
  )
  expect_lifetimes(reset, c(z_1 = 2L))
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
  expect_false("z" %in% mcp_pars(fit_gauss)$parameter)
  expect_false("offset_mu_1_" %in% mcp_pars(fit_gauss)$parameter)

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

  # Offsets apply only in segments that include them.
  offset_lifetime = mcp(
    list(
      y ~ 1 + x + offset(log(pop)),
      ~ 1 + x + offset(0),
      ~ 0 + x + offset(log(pop)),
      ~ 0 + x + offset(0)
    ),
    d_pois,
    family = poisson(),
    par_x = "x",
    sample = FALSE
  )
  offset_mu = offset_lifetime$simulate(
    offset_lifetime,
    d_pois,
    cp_1 = 7.5, cp_2 = 11.5, cp_3 = 15.5,
    Intercept_1 = 0, x_1 = 0,
    Intercept_2 = 0, x_2 = 0,
    x_3 = 0,
    x_4 = 0,
    .type = "fitted"
  )
  expect_equal(
    as.numeric(offset_mu),
    ifelse(d_pois$x < 7.5 | (d_pois$x >= 11.5 & d_pois$x < 15.5), d_pois$pop, 1)
  )
  expect_match(
    offset_lifetime$jags_code,
    "(?s)x\\[i_\\] < cp_1.*offset_mu_1_",
    perl = TRUE
  )

  # fit$simulate with offset
  sim_out = fit_gauss$simulate(
    fit_gauss,
    data.frame(x = 1:4, z = c(10, 10, 20, 20)),
    cp_1 = 2.5, Intercept_1 = 5, x_1 = 1, Intercept_2 = 10, x_2 = 2, sigma_1 = 0.001,
    .type = "fitted", .dpar = "mu"
  )
  # Segment 1 (x=1, z=10): mu = 5 + 1*1 + 10 = 16
  # Segment 1 (x=2, z=10): mu = 5 + 1*2 + 10 = 17
  # Segment 2 (x=3, z=20): mu = 10 + 2*(3-2.5) + 20 = 31
  # Segment 2 (x=4, z=20): mu = 10 + 2*(4-2.5) + 20 = 33
  expect_equal(as.numeric(sim_out), c(16, 17, 31, 33))

  # Offset in sigma
  fit_sigma = mcp(
    list(y ~ 1 + x + sigma(1 + offset(z)), ~ 0 + x + sigma(0 + offset(0))),
    d,
    sample = FALSE
  )
  expect_true(grepl("offset_sigma_1_\\[i_\\]", fit_sigma$jags_code))
  sigma_linear = fit_sigma$simulate(
    fit_sigma,
    d,
    cp_1 = 10.5,
    Intercept_1 = 0, x_1 = 0, sigma_Intercept_1 = 0,
    x_2 = 0,
    .type = "fitted", .dpar = "sigma", .scale = "linear"
  )
  expect_equal(as.numeric(sigma_linear), ifelse(d$x < 10.5, d$z, 0))

  offset_only = mcp(
    list(y ~ 0 + offset(z)), d, par_x = "x", sample = FALSE
  )
  offset_only_mu = offset_only$simulate(
    offset_only, d, sigma_1 = 1, .type = "fitted"
  )
  expect_equal(as.numeric(offset_only_mu), d$z)

  # Warning when par_x is in offset()
  expect_warning(
    mcp(list(y ~ 1 + offset(log(x)), ~ 1 + x + offset(log(x))), d_pois, family = poisson(), sample = FALSE),
    regexp = "The change-point predictor 'x' is also used in offset\\(\\)"
  )
  expect_warning(
    mcp(list(y ~ 1 + offset(x), ~ 1 + offset(x)), d, par_x = "x", sample = FALSE),
    regexp = "The change-point predictor 'x' is also used in offset\\(\\)"
  )
})

test_that("offset() is required in all segments where its dpar applies", {
  d = data.frame(x = 1:20, y = rnorm(20), z = rnorm(20))
  expect_error(
    mcp(list(y ~ 1 + offset(z), ~ 1 + x, ~ 1), d, par_x = "x", sample = FALSE),
    "`offset()` must be included in all segments of `mu` but is missing in segment 2 and 3",
    fixed = TRUE
  )
  expect_error(
    mcp(list(y ~ 1 + sigma(1 + offset(z)), ~ 1 + x), d, par_x = "x", sample = FALSE),
    "`offset()` must be included in all segments of `sigma` but is missing in segment 2",
    fixed = TRUE
  )
  # AR applies only in segments that include it
  expect_warning(
    mcp(list(y ~ 1 + ar(1, 1 + offset(z)), ~ 1 + x), d, par_x = "x", sample = FALSE),
    "Segment 2 does not include ar()", fixed = TRUE
  )
  expect_error(
    mcp(list(y ~ 1 + ar(1, 1 + offset(z)), ~ 1 + x + ar(1)), d, par_x = "x", sample = FALSE),
    "`offset()` must be included in all segments of `ar` but is missing in segment 2",
    fixed = TRUE
  )
})

test_that("formula variables are classified by role for par_x detection", {
  vars = get_formula_vars(list(
    y | weights(w) ~ 1 + x + offset(log(t)) + sigma(1 + z) + (1 + u || id),
    ~ cp((1 | id)) + 0 + ar(2, 1 + v, series = s)
  ))
  expect_equal(vars, list(
    response = "y", auxiliary = "w", predictor = c("x", "z", "u", "v"),
    group = "id", offset = "t", series = "s"
  ))

  # Numeric series columns are not change-point candidates
  data = data.frame(x = 1:12, id = rep(1:3, each = 4), y = sin(1:12))
  fit = mcp(list(y ~ 1 + x + ar(1, series = id)), data, sample = FALSE)
  expect_equal(mcp_columns(fit)$par_x, "x")
  expect_error(
    mcp(list(y ~ 1 + ar(1, series = id)), data, sample = FALSE),
    "No continuous column for change points"
  )
})

test_that("transformations use the original par_x while segment bases stay local", {
  data = data.frame(x = 0:10, y = 0)
  fit = mcp(
    list(y ~ 1, ~ 1 + sin(x) + x + I(x^2)),
    data,
    sample = FALSE
  )
  table = get_predictors(fit$model, data, fit$family, par_x = "x")

  expect_equal(table$x_degree[table$matrix_name == "sin(x)"], 0L)
  expect_equal(
    unname(table$matrix_data[[which(table$matrix_name == "sin(x)")]]),
    sin(data$x)
  )
  expect_equal(table$x_degree[table$matrix_name == "x"], 1L)
  expect_equal(table$x_degree[table$matrix_name == "I(x^2)"], 2L)

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
  predictor = dplyr::filter(get_predictors(fit$model, data, fit$family, par_x = "a.b"), .data$dpar == "mu")

  expect_equal(predictor$x_degree, 0L)
  expect_equal(unname(predictor$matrix_data[[1]]), data$axb)
})

test_that("x-free interactions are evaluated directly at zero", {
  data = data.frame(x = 0:4, z = c(2, 0, 0, 4, 5), y = 0)
  fit = mcp(list(y ~ 0 + x:z), data, par_x = "x", sample = FALSE)
  predictor = dplyr::filter(get_predictors(fit$model, data, fit$family, par_x = "x"), .data$dpar == "mu")

  expect_equal(predictor$x_degree, 1L)
  expect_equal(unname(predictor$matrix_data[[1]]), data$z)

  segmented = mcp(list(y ~ 0, ~ 0 + x:z), data, par_x = "x", sample = FALSE)
  fitted = segmented$simulate(
    segmented, data, cp_1 = 0, x.z_2 = 1, sigma_1 = 1, .type = "fitted"
  )
  expect_equal(as.numeric(fitted), data$x * data$z)
})

test_that("parameter names that end in other parameter names are evaluated correctly", {
  # gB.x_1 ends in x_1
  data = data.frame(x = c(1, 2, 3, 4, 5, 6), g = c("A", "B", "C", "A", "B", "C"), y = 0)
  fit = mcp(list(y ~ 1 + g * x), data, par_x = "x", sample = FALSE)
  fitted = fit$simulate(
    fit, data, Intercept_1 = 1, gB_1 = 2, gC_1 = 3, x_1 = 10, gB.x_1 = 20, gC.x_1 = 30,
    sigma_1 = 1, .type = "fitted"
  )
  slope = c(A = 10, B = 30, C = 40)[data$g]
  expect_equal(as.numeric(fitted), c(A = 1, B = 3, C = 4)[data$g] + slope * data$x, ignore_attr = TRUE)
})

test_that("factor punctuation produces safe parameter names", {
  data = data.frame(
    x = 1:4,
    state = factor(
      c("reference", "a/b", "reference", "a/b"),
      levels = c("reference", "a/b")
    ),
    y = 0
  )
  fit = mcp(list(y ~ state), data, par_x = "x", sample = FALSE)
  predictor = get_fit_model_tables(fit)$predictors
  slash_row = grepl("a/b", predictor$matrix_name, fixed = TRUE)
  safe_name = predictor$name[slash_row]

  expect_true(length(safe_name) == 1 && grepl("^[A-Za-z][A-Za-z0-9_.]*$", safe_name))
  args = list(fit = fit, newdata = data, Intercept_1 = 0, sigma_1 = 1, .type = "fitted")
  args[[safe_name]] = 1
  expect_no_error(do.call(fit$simulate, args))
})

# Columns of the design matrix evaluated on `newdata` for the coefficients `names`
design_columns = function(fit, newdata, names) {
  predictors = get_fit_model_tables(fit)$predictors
  get_model_data(fit, newdata)$matrix[, predictors$matrix_col[match(names, predictors$name)], drop = FALSE]
}

test_that("prediction reuses fitted factor encodings", {
  data = data.frame(
    x = 1:6,
    y = 1:6,
    state = ordered(rep(c("low", "middle", "high"), 2))
  )
  ordered_fit = mcp(list(y ~ state), data, par_x = "x", sample = FALSE)
  ordered_tables = get_predictor_tables(ordered_fit$model, data, ordered_fit$family, par_x = "x")
  ordered_matrix = get_predictor_matrix(ordered_tables$predictors)  # No group-level effects
  ordered_new = get_model_data(ordered_fit, transform(data, state = as.character(state)))$matrix
  expect_equal(unname(ordered_new), unname(ordered_matrix))

  custom_data = transform(data, state = factor(state))
  contrasts(custom_data$state) = stats::contr.sum(3)
  custom_fit = suppressWarnings(mcp(list(y ~ state), custom_data, par_x = "x", sample = FALSE))
  custom_tables = suppressWarnings(get_predictor_tables(custom_fit$model, custom_data, custom_fit$family, par_x = "x"))
  custom_matrix = get_predictor_matrix(custom_tables$predictors)  # No group-level effects
  custom_new = get_model_data(custom_fit, transform(custom_data, state = as.character(state)))$matrix
  expect_equal(unname(custom_new), unname(custom_matrix))
})

test_that("prediction is independent of later contrast options", {
  old_options = options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old_options), add = TRUE)
  data = data.frame(
    x = 1:6,
    y = 1:6,
    state = factor(rep(c("a", "b", "c"), 2))
  )
  fit = mcp(list(y ~ state), data, par_x = "x", sample = FALSE)
  fitted_tables = get_predictor_tables(fit$model, data, fit$family, par_x = "x")
  fitted_matrix = get_predictor_matrix(fitted_tables$predictors)  # No group-level effects

  options(contrasts = c("contr.sum", "contr.poly"))
  new_matrix = get_model_data(fit, data)$matrix
  expect_equal(unname(new_matrix), unname(fitted_matrix))
})

test_that("data-derived bases reuse their fitted specification", {
  data = data.frame(x = 1:8, y = 1:8)
  newdata = data.frame(x = 9:10)
  coefficients = function(fit) setdiff(get_fit_model_tables(fit)$predictors$name, c("Intercept_1", "sigma_1"))

  scale_fit = mcp(list(y ~ scale(x)), data, sample = FALSE)
  expect_equal(
    as.numeric(design_columns(scale_fit, newdata, coefficients(scale_fit))),
    as.numeric(scale(newdata$x, center = mean(data$x), scale = stats::sd(data$x)))
  )

  poly_basis = stats::poly(data$x, 2)
  poly_fit = mcp(list(y ~ poly(x, 2)), data, sample = FALSE)
  expect_equal(
    unname(design_columns(poly_fit, newdata, coefficients(poly_fit))),
    unname(stats::predict(poly_basis, newdata$x))
  )

  spline_basis = splines::ns(data$x, df = 3)
  spline_fit = mcp(list(y ~ splines::ns(x, df = 3)), data, sample = FALSE)
  expect_equal(
    unname(design_columns(spline_fit, newdata, coefficients(spline_fit))),
    matrix(
      as.numeric(stats::predict(spline_basis, newdata$x)),
      nrow = nrow(newdata)
    )
  )

  bs_newdata = data.frame(x = 2:3)
  bs_basis = splines::bs(data$x, df = 3)
  bs_fit = mcp(list(y ~ splines::bs(x, df = 3)), data, sample = FALSE)
  expect_equal(
    unname(design_columns(bs_fit, bs_newdata, coefficients(bs_fit))),
    matrix(
      as.numeric(stats::predict(bs_basis, bs_newdata$x)),
      nrow = nrow(bs_newdata)
    )
  )
})


test_that("fits do not include the calling environment, so formula functions must be reachable from .GlobalEnv", {
  data = data.frame(
    x = 1:8,
    y = c(1, 2, 1, 2, 5, 6, 5, 6),
    id = rep(c("a", "b"), 4)
  )
  local_fit = function(formula) {
    sq_local = function(x) x^2
    big_local = numeric(1e6)
    mcp(list(stats::as.formula(formula, env = environment())), data, par_x = "x", sample = FALSE, quiet = TRUE)
  }

  # Locally defined functions err early in mcp(), not later in fitted() etc.
  for (formula in c("y ~ sq_local(x)", "y ~ 1 + sigma(sq_local(x))", "y ~ 1 + (sq_local(x) || id)", "y ~ 1 + ar(1, sq_local(x))"))
    expect_error(local_fit(formula), "not available from the global environment: `sq_local`")

  # Base and namespaced functions work, and the calling environment is not included
  fit = local_fit("y ~ log(x) + splines::ns(x, 2)")
  expect_true(identical(environment(fit$model[[1]]), globalenv()))
  expect_true(all(vapply(get_fit_model_tables(fit)$design_specs, function(spec) identical(environment(spec$terms), globalenv()), logical(1))))
  expect_lt(length(serialize(fit, NULL)), 1e6)
  expect_equal(as.numeric(design_columns(fit, data.frame(x = 9:10), "logx_1")), log(9:10))
})


test_that("each fitted design specification is evaluated once in get_model_data", {
  d = data.frame(x = 1:10, y = 1:10, z = 1:10)
  fit = mcp(list(y ~ 1 + x + offset(z), ~ 1 + x + offset(2 * z)), data = d, family = poisson(), sample = FALSE)
  expect_equal(length(get_fit_model_tables(fit)$design_specs), 2)

  eval_count = 0
  trace_env = environment(get_model_data)
  orig_get_fitted_design = trace_env$get_fitted_design
  unlockBinding("get_fitted_design", trace_env)
  on.exit({
    trace_env$get_fitted_design = orig_get_fitted_design
    lockBinding("get_fitted_design", trace_env)
  }, add = TRUE)
  trace_env$get_fitted_design = function(...) {
    eval_count <<- eval_count + 1
    orig_get_fitted_design(...)
  }

  res = get_model_data(fit, d)
  expect_equal(eval_count, 2)
  expect_equal(res$offsets, list(offset_mu_1_ = d$z, offset_mu_2_ = 2 * d$z))
  expect_equal(ncol(res$matrix), 4)  # Intercept_1, x_1, Intercept_2, x_2
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
