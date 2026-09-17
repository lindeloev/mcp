test_that("shared slopes retain implicit intercepts and shared intercepts retain factor coding", {
  data = data.frame(x = 1:12, y = 0, g = factor(rep(c("a", "b", "c"), 4)))
  fit = mcp(list(y ~ 1 + x, ~ same(x)), data, par_x = "x", sample = FALSE)
  expect_equal(as.numeric(fit$simulate(fit, data, cp_1 = 6, Intercept_1 = 1,
    Intercept_2 = 10, x_1 = 2, sigma_1 = 1, .type = "fitted")),
    ifelse(data$x < 6, 1 + 2 * data$x, 10 + 2 * (data$x - 6)))

  for (rhs in list(~ same(1) + g, ~ 0 + same(1) + g)) {
    fit = mcp(list(y ~ 1, rhs), data, par_x = "x", sample = FALSE)
    expect_false("ga_2" %in% mcp_pars(fit)$name)
    expect_equal(as.numeric(fit$simulate(fit, data, cp_1 = 6, Intercept_1 = 2,
      gb_2 = 3, gc_2 = 5, sigma_1 = 1, .type = "fitted")),
      2 + (data$x >= 6) * c(0, 3, 5)[as.integer(data$g)])
  }
  expect_error(mcp(list(y ~ 1, ~ 1 + same(1)), data, par_x = "x", sample = FALSE), "competing")
})


test_that("AR sharing resolves every lag and preserves current segment transitions", {
  data = data.frame(x = 1:12, y = sin(1:12))
  fit = mcp(list(y ~ 1 + ar(2), ~ 1 + ar(2, same(1), threshold = 0.2)),
    data, par_x = "x", sample = FALSE)
  args = c(as.list(add_rhs_predictors(data, fit)),
    lapply(list(cp_1 = 6, Intercept_1 = 1, Intercept_2 = 2, sigma_1 = 1, ar1_1 = 0.3, ar2_1 = -0.1), rep, 12))
  tables = get_fit_model_tables(fit)
  values = evaluate_model_dpars(fit, args,
    paste0(".pred_", get_predictor_design_names(tables$predictors, tables$group_effects)))
  expect_equal(values$ar1_, rep(0.3, 12))
  expect_equal(values$ar2_, rep(-0.1, 12))
  expect_equal(values$garma_threshold_, ifelse(data$x < 6, 0.1, 0.2))
  expect_error(mcp(list(y ~ 1 + ar(1), ~ 1 + ar(2, same(1))),
    data, par_x = "x", sample = FALSE), "lag 2")

  # A declared no-intercept lag keeps its earlier level and local-x endpoint
  fit = mcp(list(y ~ 1 + ar(1, 1 + x), ~ 0 + ar(1, 0 + same(x))),
    data, par_x = "x", sample = FALSE)
  values = fit$simulate(fit, data, cp_1 = 6, Intercept_1 = 0, sigma_1 = 1,
    ar1_1 = 0.1, ar1_x_1 = 0.01, .type = "fitted", .dpar = "ar1", .arma = FALSE)
  expect_equal(as.numeric(values), 0.1 + 0.01 * data$x)
})


test_that("same() reuses population definitions with current local x coordinates", {
  data = data.frame(
    x = 1:8,
    z = c(0, 1, 0, 1, 0, 1, 0, 1),
    y = 0
  )
  fit = mcp(
    list(y ~ 1 + x + z, ~ 0 + same(x) + same(z)),
    data, par_x = "x", sample = FALSE
  )
  tables = get_fit_model_tables(fit)

  expect_setequal(mcp_pars(fit)$name, c("cp_1", "Intercept_1", "x_1", "z_1", "sigma_1"))
  expect_equal(
    tables$predictors[tables$predictors$segment == 2, c("code_name", "definition_segment")],
    tibble::tibble(code_name = c("x_1", "z_1"), definition_segment = c(1L, 1L))
  )
  expect_equal(
    as.numeric(fit$simulate(
      fit, data, cp_1 = 4, Intercept_1 = 1, x_1 = 2, z_1 = 3, sigma_1 = 1,
      .type = "fitted"
    )),
    c(3, 8, 7, 12, 11, 16, 15, 20)
  )
})


test_that("a fitted shared-coefficient model evaluates its likelihood", {
  set.seed(1)
  data = data.frame(x = 1:12, z = rep(c(0, 1), 6))
  data$y = 1 + 0.2 * data$x + 0.5 * data$z + stats::rnorm(nrow(data), 0, 0.2)
  fit = suppressWarnings(mcp(
    list(y ~ 1 + x + z, ~ 0 + same(x) + same(z)),
    data, par_x = "x", iter = 40, warmup = 20, chains = 1, quiet = TRUE
  ))

  likelihood = log_lik(fit, summary = FALSE)
  expect_equal(dim(likelihood), c(40L, nrow(data)))
  expect_true(all(is.finite(likelihood)))
  draws = as.matrix(coda::as.mcmc(fit))
  expected = vapply(seq_len(nrow(data)), function(i) {
    stats::dnorm(data$y[i], draws[, "Intercept_1"] + draws[, "x_1"] * data$x[i] +
      draws[, "z_1"] * data$z[i], draws[, "sigma_1"], log = TRUE)
  }, numeric(nrow(draws)))
  expect_equal(unname(likelihood), expected)
})


test_that("sharing chains reuse fitted bases and explicit sources bridge gaps", {
  data = data.frame(x = 1:20, z = sin(1:20), y = 0)
  fit = mcp(list(y ~ 1 + poly(z, 2), ~ 0 + same(poly(z, 2)),
    ~ 0 + same(poly(z, 2))), data, par_x = "x", sample = FALSE)
  tables = get_fit_model_tables(fit)
  source = tables$predictors[tables$predictors$term_key == "poly(z, 2)" & tables$predictors$segment == 1, ]
  newdata = data.frame(x = c(3, 9, 17), z = c(-0.4, 0.2, 0.7))
  args = c(list(fit = fit, newdata = newdata, cp_1 = 6, cp_2 = 14,
    Intercept_1 = 2, sigma_1 = 1, .type = "fitted"),
    stats::setNames(list(3, -1), source$code_name))
  basis = predict(poly(data$z, 2), newdata$z)
  expect_equal(as.numeric(do.call(fit$simulate, args)), as.numeric(2 + basis %*% c(3, -1)))

  expect_error(mcp(list(y ~ 1 + z, ~ 0, ~ 0 + same(z)), data,
    par_x = "x", sample = FALSE), "cannot find")
  fit = mcp(list(y ~ 1 + z, ~ 0, ~ 0 + same(z, as = 1)), data,
    par_x = "x", sample = FALSE)
  expect_equal(as.numeric(fit$simulate(fit, data, cp_1 = 6, cp_2 = 14,
    Intercept_1 = 2, z_1 = 3, sigma_1 = 1, .type = "fitted")),
    2 + 3 * data$z * (data$x < 6 | data$x >= 14))
})


test_that("same() keeps fitted contrasts and current segment offsets", {
  data = data.frame(x = 1:12, y = 0, g = factor(rep(c("a", "b", "c"), 4)), z = 1:12 / 10)
  old_options = options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old_options), add = TRUE)
  fit = mcp(list(y ~ 1 + g + offset(z), ~ 0 + same(g), ~ same(1, as = 1) + offset(2 * z)),
    data, par_x = "x", sample = FALSE)
  options(old_options)
  values = fit$simulate(fit, data, cp_1 = 4, cp_2 = 8, Intercept_1 = 2,
    g1_1 = 3, g2_1 = 5, sigma_1 = 1, .type = "fitted")
  expected = 2 + (data$x < 8) * c(3, 5, -8)[as.integer(data$g)] +
    ifelse(data$x < 4, data$z, ifelse(data$x >= 8, 2 * data$z, 0))
  expect_equal(as.numeric(values), expected)
})


test_that("same() preserves source factor coding and supports explicit sources", {
  data = data.frame(
    x = 1:9,
    group = factor(rep(c("a", "b", "c"), 3)),
    y = 0
  )
  fit = mcp(
    list(y ~ 1 + group, ~ 0 + same(group), ~ 0 + same(group, as = 1)),
    data, par_x = "x", sample = FALSE
  )
  tables = get_fit_model_tables(fit)

  expect_setequal(mcp_pars(fit)$name, c("cp_1", "cp_2", "Intercept_1", "groupb_1", "groupc_1", "sigma_1"))
  expect_true(all(tables$predictors$code_name[tables$predictors$segment > 1 & tables$predictors$term_key == "group"] %in% c("groupb_1", "groupc_1")))
  expect_equal(unique(tables$predictors$definition_segment[tables$predictors$segment == 3]), 1L)
  expect_no_error(add_rhs_predictors(data[1:3, ], fit))

  interaction_fit = mcp(
    list(y ~ 1 + x * group, ~ 0 + same(x * group)),
    data, par_x = "x", sample = FALSE
  )
  expect_false(any(grepl("_2$", mcp_pars(interaction_fit)$name)))
})


test_that("same() supports distributional and AR/MA coefficient formulas", {
  data = data.frame(x = 1:20, y = rnorm(20))
  dpar_fit = mcp(
    list(y ~ 1 + sigma(1 + x), ~ 1 + sigma(0 + same(x))),
    data, par_x = "x", sample = FALSE
  )
  arma_fit = mcp(
    list(y ~ 1 + ar(1, 1 + x), ~ 0 + ar(1, 0 + same(x))),
    data, par_x = "x", sample = FALSE
  )

  expect_false("sigma_x_2" %in% mcp_pars(dpar_fit)$name)
  expect_equal(
    unique(get_fit_model_tables(dpar_fit)$predictors$definition_segment[
      get_fit_model_tables(dpar_fit)$predictors$dpar == "sigma" &
        get_fit_model_tables(dpar_fit)$predictors$term_key == "x"
    ]),
    1L
  )
  expect_false("ar1_x_2" %in% mcp_pars(arma_fit)$name)
})


test_that("same() validates sources, competing declarations, and group blocks", {
  data = data.frame(x = 1:8, y = 0, id = rep(1:2, each = 4))

  expect_error(mcp(list(y ~ same(1)), data, par_x = "x", sample = FALSE), "segment 1")
  expect_error(mcp(list(y ~ 1, ~ same(x)), data, par_x = "x", sample = FALSE), "cannot find")
  expect_error(mcp(list(y ~ 1 + x, ~ x + same(x)), data, par_x = "x", sample = FALSE), "both bare and shared")
  expect_error(
    mcp(list(y ~ 1 + x, ~ 0 + same(x)), data, par_x = "x", prior = list(x_2 = "dnorm(0, 1)"), sample = FALSE),
    "Set prior on `x_1` instead of `x_2`"
  )
  expect_no_error(mcp(list(y ~ 1 + (1 | id), ~ same((1 | id))), data, par_x = "x", sample = FALSE))
})


test_that("explicit mu() syntax produces identical model tables to standard bare syntax", {
  data = data.frame(x = 1:10, y = 1:10, id = rep(1:2, 5))
  fit_explicit = mcp(
    list(y ~ mu(1 + x + (1 | id)) + sigma(1 + x), ~ mu(0 + same(x))),
    data = data, par_x = "x", sample = FALSE
  )
  fit_bare = mcp(
    list(y ~ 1 + x + (1 | id) + sigma(1 + x), ~ 0 + same(x)),
    data = data, par_x = "x", sample = FALSE
  )
  expect_identical(mcp_pars(fit_explicit), mcp_pars(fit_bare))
})


test_that("shared group blocks define one hierarchy and simulate one deviation per group", {
  data = data.frame(x = 1:12, y = 0, id = rep(c("a", "b"), 6))
  fit = mcp(list(y ~ 1 + (1 | id), ~ 0 + same((1 || id))), data, par_x = "x", sample = FALSE)
  expect_equal(mcp_pars(fit, scope = "group")$name, "Intercept_1_id")
  expect_equal(mcp_pars(fit, role = "group_sd")$name, "Intercept_1_id_sd")
  for (ids in list(data$id, rep(c("new-a", "new-b"), 6))) {
    newdata = data
    newdata$id = ids
    y = fit$simulate(fit, newdata, cp_1 = 6, Intercept_1 = 2, Intercept_1_id_sd = 1,
      sigma_1 = 1, .type = "fitted")
    expect_equal(length(unique(y)), 2L)
    expect_equal(as.numeric(y[1:2]), as.numeric(y[7:8]))
  }
  expect_error(mcp(fit$model, data, par_x = "x", sample = FALSE,
    prior = list(Intercept_2_id_sd = "dnorm(0,1)")), "Set prior on `Intercept_1_id_sd` instead of `Intercept_2_id_sd`")
  bare = mcp(list(y ~ 1 + (1 | id), ~ 0 + (1 | id)), data, par_x = "x", sample = FALSE)
  expect_equal(nrow(mcp_pars(bare, role = "group_sd")), 2L)
})


test_that("whole group slopes reuse local coordinates, endpoints, and explicit sources", {
  data = data.frame(x = 1:12, y = 0, id = rep(c("a", "b"), 6))
  fit = mcp(list(y ~ 1 + (0 + x || id), ~ 0 + same((0 + x || id)), ~ 0),
    data, par_x = "x", sample = FALSE)
  values = fit$simulate(fit, data, cp_1 = 4, cp_2 = 8, Intercept_1 = 1,
    x_1_id_sd = 1, sigma_1 = 1, .type = "fitted")
  b = attr(values, "simulated")$x_1_id
  expect_equal(as.numeric(values), 1 + b * pmin(data$x, 8))
  expect_error(mcp(list(y ~ 1 + (1 + x || id), ~ 0 + same((0 + x || id))),
    data, par_x = "x", sample = FALSE), "complete active group block")
  expect_error(mcp(list(y ~ 1 + (1 | id), ~ 0, ~ same((1 | id))),
    data, par_x = "x", sample = FALSE), "complete active group block")
  expect_no_error(mcp(list(y ~ 1 + (1 | id), ~ 0, ~ same((1 | id), as = 1)),
    data, par_x = "x", sample = FALSE))
  expect_error(mcp(list(y ~ 1 + (1 | id), ~ (0 + x || id) + same((1 | id))),
    data, par_x = "x", sample = FALSE), "Only one predictor group-level term")
})


test_that("posterior prediction and likelihood reuse fitted and new group deviations", {
  data = data.frame(x = 1:18, id = rep(c("a", "b", "c"), 6), y = sin(1:18))
  fit = suppressWarnings(mcp(list(y ~ 1 + (1 | id), ~ 0 + same((1 | id))),
    data, par_x = "x", iter = 40, warmup = 30, chains = 1, diagnostics = FALSE, quiet = TRUE))
  draws = as.matrix(coda::as.mcmc(fit))
  expected = vapply(seq_len(nrow(data)), function(i) stats::dnorm(data$y[i],
    draws[, "Intercept_1"] + draws[, paste0("Intercept_1_id[", data$id[i], "]")],
    draws[, "sigma_1"], log = TRUE), numeric(nrow(draws)))
  expect_equal(unname(log_lik(fit, summary = FALSE)), expected)
  for (id in c("a", "new")) {
    nd = data.frame(x = c(2, 17), id = id)
    pred = fitted(fit, newdata = nd, summary = FALSE)
    expect_equal(pred$.epred[pred$data_row == 1], pred$.epred[pred$data_row == 2])
  }
})

