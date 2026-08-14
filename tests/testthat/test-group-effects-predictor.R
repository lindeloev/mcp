group_data = data.frame(
  x = 1:8,
  y = c(1, 2, 1, 2, 5, 6, 5, 6),
  id = rep(c("a", "b"), 4),
  site = rep(c("north", "south"), each = 4)
)

coefficient_data = data.frame(
  x = 1:12,
  y = rep(c(2, 4, 6), 4),
  id = rep(c("a", "b", "c"), 4),
  condition = factor(rep(c("A", "B", "C"), each = 4)),
  z = seq(-1, 1, length.out = 12)
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
  expect_true(effect$sd_name %in% names(formals(fit$simulate)))
  expect_false(effect$name %in% names(formals(fit$simulate)))

  set.seed(42)
  simulated = fit$simulate(
    fit,
    group_data,
    cp_1 = 2.5,
    cp_2 = 5.5,
    Intercept_1 = 0,
    Intercept_2 = 100,
    Intercept_3 = 200,
    sigma_1 = 1,
    Intercept_1_id_sd = 10,
    .type = "fitted"
  )
  deviation = attr(simulated, "simulated")$Intercept_1_id
  expect_equal(
    as.numeric(simulated),
    c(deviation[1:2], 100 + deviation[3:5], rep(200, 3))
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
  expect_equal(fit$pars$group, effects$name)
  expect_equal(
    fit$pars$population,
    c("Intercept_1_id_sd", "sigma_1_id_sd")
  )
  expect_equal(unpack_group_effects(fit, "predictor")$pars, effects$name)
  expect_null(unpack_group_effects(fit, "cp")$pars)

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

  set.seed(42)
  simulated = fit$simulate(
    fit,
    group_data,
    Intercept_1_id_sd = 1,
    sigma_1_id_sd = 0.5,
    .type = "fitted"
  )
  simulation = attr(simulated, "simulated")
  expect_equal(as.numeric(simulated), simulation$Intercept_1_id)
  expect_equal(simulation$Intercept_1_id_sd, 1)
  expect_equal(simulation$sigma_1_id_sd, 0.5)
  expect_equal(
    simulation$Intercept_1_id[group_data$id == "a"],
    rep(simulation$Intercept_1_id[1], 4)
  )
})


test_that("double-bar terms expand into independent group coefficients", {
  fit = mcp(
    list(
      y ~ 1 + (condition || id),
      ~ 1,
      ~ 1 + (0 | id)
    ),
    coefficient_data,
    par_x = "x",
    sample = FALSE
  )
  effects = get_fit_model_tables(fit)$group_effects %>%
    dplyr::filter(.data$part == "predictor")

  expect_equal(
    effects$name,
    c("Intercept_1_id", "conditionB_1_id", "conditionC_1_id")
  )
  expect_equal(
    effects$population_name,
    c("Intercept_1", NA_character_, NA_character_)
  )
  expect_equal(effects$par_type, c("Intercept", "dummy", "dummy"))
  expect_equal(effects$next_segment, rep(3L, 3))
  expect_false(any(effects$correlated))
  expect_true(all(effects$sd_name %in% names(fit$prior)))
  expect_match(fit$jags_code, "conditionB_1_id\\[id_\\] ~")
  expect_match(fit$jags_code, "conditionC_1_id\\[id_\\] ~")

  set.seed(42)
  simulated = fit$simulate(
    fit,
    coefficient_data,
    cp_1 = 4.5,
    cp_2 = 8.5,
    Intercept_1 = 1,
    Intercept_2 = 2,
    Intercept_3 = 3,
    sigma_1 = 1,
    Intercept_1_id_sd = 1,
    conditionB_1_id_sd = 10,
    conditionC_1_id_sd = 100,
    .type = "fitted"
  )
  simulation = attr(simulated, "simulated")
  intercept_by_id = simulation$Intercept_1_id
  condition_b_by_id = simulation$conditionB_1_id
  condition_c_by_id = simulation$conditionC_1_id
  expected = ifelse(
    coefficient_data$x < 4.5,
    1 + intercept_by_id +
      (coefficient_data$condition == "B") * condition_b_by_id +
      (coefficient_data$condition == "C") * condition_c_by_id,
    ifelse(
      coefficient_data$x < 8.5,
      2 + intercept_by_id +
        (coefficient_data$condition == "B") * condition_b_by_id +
        (coefficient_data$condition == "C") * condition_c_by_id,
      3
    )
  )
  expect_equal(as.numeric(simulated), as.numeric(expected))
})


test_that("double-bar terms support no-intercept factors and numeric slopes", {
  factor_fit = mcp(
    list(y ~ 0 + (0 + condition || id)),
    coefficient_data,
    par_x = "x",
    sample = FALSE
  )
  factor_effects = get_fit_model_tables(factor_fit)$group_effects
  expect_equal(
    factor_effects$name,
    c("conditionA_1_id", "conditionB_1_id", "conditionC_1_id")
  )
  expect_true(all(is.na(factor_effects$population_name)))
  expect_true(all(factor_effects$par_type == "dummy"))

  slope_fit = mcp(
    list(y ~ 1 + (0 + z || id)),
    coefficient_data,
    par_x = "x",
    sample = FALSE
  )
  slope_effect = get_fit_model_tables(slope_fit)$group_effects
  expect_equal(slope_effect$name, "z_1_id")
  expect_equal(slope_effect$par_type, "slope")
  expect_true(is.na(slope_effect$population_name))
  expect_match(slope_fit$jags_code, "z_1_id\\[id_\\]")
})


test_that("double-bar coefficient blocks are replaced together", {
  fit = mcp(
    list(
      y ~ 1 + (condition || id),
      ~ 1 + (1 || id)
    ),
    coefficient_data,
    par_x = "x",
    sample = FALSE
  )
  effects = get_fit_model_tables(fit)$group_effects

  expect_equal(
    effects$name,
    c(
      "Intercept_1_id", "conditionB_1_id", "conditionC_1_id",
      "Intercept_2_id"
    )
  )
  expect_equal(effects$next_segment, c(2L, 2L, 2L, NA_integer_))
})


test_that("factor group coefficients sample and predict by level", {
  set.seed(43)
  data = expand.grid(
    id = c("a", "b", "c"),
    condition = factor(c("A", "B", "C")),
    replicate = 1:2
  )
  data$x = seq_len(nrow(data))
  data$y = 2 +
    c(a = -1, b = 0, c = 1)[data$id] +
    (data$condition == "B") * c(a = 0.5, b = 1, c = 1.5)[data$id] +
    stats::rnorm(nrow(data), 0, 0.3)

  expect_warning({
    fit = mcp(
      list(y ~ 1 + (condition || id)),
      data,
      par_x = "x",
      chains = 1,
      adapt = 20,
      iter = 20,
      diagnostics = FALSE,
      quiet = TRUE
    )
  }, "Adaptation incomplete", fixed = TRUE)

  effects = ranef(fit)
  expect_equal(nrow(effects), 3 * length(unique(data$id)))
  expect_true(all(
    c("Intercept_1_id[a]", "conditionB_1_id[a]", "conditionC_1_id[a]") %in%
      effects$name
  ))

  fitted_values = fitted(
    fit, summary = FALSE, varying = "predictor", ndraws = 2
  )
  expect_true(all(
    c(
      "id", "condition", "Intercept_1_id",
      "conditionB_1_id", "conditionC_1_id"
    ) %in% names(fitted_values)
  ))
})


test_that("predictor group effects use existing prediction selectors", {
  set.seed(42)
  data = data.frame(
    x = rep(1:6, 2),
    id = rep(c("a", "b"), each = 6)
  )
  data$y = 3 + ifelse(data$id == "a", -1, 1) + stats::rnorm(12, 0, 0.2)
  expect_warning({
    fit = mcp(
      list(y ~ 1 + (1 | id)),
      data,
      par_x = "x",
      chains = 1,
      adapt = 20,
      iter = 20,
      diagnostics = FALSE,
      quiet = TRUE
    )
  }, "Adaptation incomplete", fixed = TRUE)

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
    "require `\\|\\|`"
  )
  expect_error(
    mcp(
      list(y ~ 1 + (site | id)),
      group_data,
      par_x = "x",
      sample = FALSE
    ),
    "require `\\|\\|`"
  )
  expect_error(
    mcp(
      list(y ~ 1 + (1 || id) + (0 + site || id)),
      group_data,
      par_x = "x",
      sample = FALSE
    ),
    "Only one predictor group-level term"
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
