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
    "use `x_1`"
  )
  expect_error(mcp(list(y ~ 1 + (1 | id), ~ same((1 | id))), data, par_x = "x", sample = FALSE), "group-effect blocks")
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
