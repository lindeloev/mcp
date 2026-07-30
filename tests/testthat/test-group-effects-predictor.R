group_data = data.frame(
  x = 1:8,
  y = c(1, 2, 1, 2, 5, 6, 5, 6),
  id = rep(c("a", "b"), 4),
  site = rep(c("north", "south"), each = 4)
)


test_that("predictor group intercepts carry, replace, and turn off", {
  fit = mcp(
    list(
      y ~ 1 + (1 | id),
      ~ 1,
      ~ 1 + (0 | id)
    ),
    group_data,
    par_x = "x",
    sample = FALSE
  )
  effect = get_fit_model_tables(fit)$group_effects %>%
    dplyr::filter(.data$part == "predictor")

  expect_equal(effect$name, "Intercept_1_id")
  expect_equal(effect$population_name, "Intercept_1")
  expect_equal(effect$next_segment, 3L)
  expect_equal(effect$dpar, "mu")
  expect_equal(effect$matrix_col, nrow(get_fit_model_tables(fit)$predictors) + 1L)

  simulated = fit$simulate(
    fit,
    group_data,
    cp_1 = 2.5,
    cp_2 = 5.5,
    Intercept_1 = 0,
    Intercept_2 = 100,
    Intercept_3 = 200,
    sigma_1 = 1,
    Intercept_1_id = 10,
    .type = "fitted"
  )
  expect_equal(
    as.numeric(simulated),
    c(10, 10, 110, 110, 110, 200, 200, 200)
  )

  replacement = mcp(
    list(y ~ 1 + (1 | id), ~ 1 + (1 || id)),
    group_data,
    par_x = "x",
    sample = FALSE
  )
  effects = get_fit_model_tables(replacement)$group_effects %>%
    dplyr::filter(.data$part == "predictor")
  expect_equal(effects$name, c("Intercept_1_id", "Intercept_2_id"))
  expect_equal(effects$next_segment, c(2L, NA_integer_))
  expect_equal(effects$correlated, c(TRUE, FALSE))
})


test_that("predictor group intercepts work for family dpars and group-only formulas", {
  fit = mcp(
    list(y ~ 0 + (1 | id) + sigma(0 + (1 || id))),
    group_data,
    par_x = "x",
    sample = FALSE
  )
  effects = get_fit_model_tables(fit)$group_effects

  expect_equal(effects$name, c("Intercept_1_id", "sigma_1_id"))
  expect_true(all(is.na(effects$population_name)))
  expect_equal(effects$dpar, c("mu", "sigma"))
  expect_equal(fit$pars$varying, effects$name)
  expect_equal(
    fit$pars$population,
    c("Intercept_1_id_sd", "sigma_1_id_sd")
  )
  expect_equal(unpack_varying(fit, "predictor")$pars, effects$name)
  expect_null(unpack_varying(fit, "cp")$pars)

  expect_named(
    fit$prior,
    c(
      "Intercept_1_id_sd", "sigma_1_id_sd",
      "Intercept_1_id", "sigma_1_id"
    ),
    ignore.order = TRUE
  )
  expect_match(fit$jags_code, "Intercept_1_id\\[id_\\] ~")
  expect_match(fit$jags_code, "sigma_1_id\\[id_\\] ~")
  expect_false(grepl("Intercept_1_id_uncentered", fit$jags_code, fixed = TRUE))

  simulated = fit$simulate(
    fit,
    group_data,
    Intercept_1_id = rep(c(1, 2), 4),
    sigma_1_id = rep(log(c(1, 2)), 4),
    .type = "fitted"
  )
  expect_equal(as.numeric(simulated), rep(c(1, 2), 4))
})


test_that("predictor group effects use existing prediction selectors", {
  set.seed(42)
  data = data.frame(
    x = rep(1:6, 2),
    id = rep(c("a", "b"), each = 6)
  )
  data$y = 3 + ifelse(data$id == "a", -1, 1) + stats::rnorm(12, 0, 0.2)
  fit = suppressWarnings(mcp(
    list(y ~ 1 + (1 | id)),
    data,
    par_x = "x",
    chains = 1,
    adapt = 20,
    iter = 20
  ))

  varying = fitted(
    fit, summary = FALSE, varying = "predictor", ndraws = 2
  )
  expect_true(all(c("id", "Intercept_1_id") %in% names(varying)))

  population = fitted(
    fit,
    newdata = data.frame(x = 1:3),
    summary = FALSE,
    varying = FALSE,
    ndraws = 2
  )
  expect_false("id" %in% names(population))
  expect_true(all(population$Intercept_1_id == 0))

  effects = ranef(fit)
  expect_equal(nrow(effects), length(unique(data$id)))
  expect_true(all(grepl("^Intercept_1_id\\[", effects$name)))
})


test_that("unsupported predictor group structures fail explicitly", {
  expect_error(
    mcp(
      list(y ~ 1 + (1 + x | id)),
      group_data,
      par_x = "x",
      sample = FALSE
    ),
    "intercepts only"
  )
  expect_error(
    mcp(
      list(y ~ 1 + (site || id)),
      group_data,
      par_x = "x",
      sample = FALSE
    ),
    "intercepts only"
  )
  expect_error(
    mcp(
      list(y ~ 1 + ar(1, 1 + (1 | id))),
      group_data,
      par_x = "x",
      sample = FALSE
    ),
    "inside ar\\(\\)"
  )
  expect_error(
    mcp(
      list(y ~ 1 + (1 | id:site)),
      group_data,
      par_x = "x",
      sample = FALSE
    ),
    "plain data-column names"
  )
})
