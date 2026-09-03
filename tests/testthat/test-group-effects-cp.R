cp_group_data = data.frame(
  x = rep(1:10, 3),
  y = seq_len(30),
  id = factor(rep(letters[1:3], each = 10)),
  site = factor(rep(LETTERS[1:3], 10))
)


test_that("group-level change points use one grouping factor", {
  expect_error(
    mcp(
      list(y ~ 1, 1 + (1 | id) ~ 1, 1 + (1 | site) ~ 1),
      cp_group_data,
      par_x = "x",
      sample = FALSE
    ),
    "must use the same grouping factor"
  )
})


test_that("group-level change points use absolute locations and sequential bounds", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    sample = FALSE
  )

  expect_match(fit$prior$cp_2_id, "cp_1_id\\[id_\\]")
  expect_match(fit$jags_code, "cp_1_id_location\\[id_\\] ~ dnorm\\(cp_1")
  expect_match(fit$jags_code, "cp_1_id\\[id_\\] = cp_1_id_location\\[id_\\] - cp_1")
  expect_match(fit$jags_code, "cp_2_id_location\\[id_\\].*cp_1_id_location\\[id_\\]")
  expect_false(grepl("_uncentered", fit$jags_code, fixed = TRUE))
  expect_false(grepl("cp_order_", fit$jags_code, fixed = TRUE))
})


test_that("custom normal change-point hierarchies retain model bounds", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    prior = list(cp_1_id = "dnorm(0, cp_1_sd)"),
    sample = FALSE
  )
  expect_match(fit$prior$cp_1_id, "T\\(1 - cp_1, 10 - cp_1\\)")
  expect_match(fit$jags_code, "cp_1_id_location\\[id_\\] ~ dnorm\\(cp_1")
})


test_that("simulation generates bounded hierarchical change-point locations", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    sample = FALSE
  )

  expect_true("cp_1_sd" %in% names(formals(fit$simulate)))
  expect_false("cp_1_id" %in% names(formals(fit$simulate)))
  set.seed(42)
  simulated = fit$simulate(
    fit,
    cp_group_data,
    cp_1 = 5,
    cp_1_sd = 1,
    Intercept_1 = 0,
    Intercept_2 = 0,
    sigma_1 = 1,
    .type = "fitted"
  )
  deviations = attr(simulated, "simulated")$cp_1_id
  by_group = vapply(split(deviations, cp_group_data$id), unique, numeric(1))
  expect_false(isTRUE(all.equal(mean(by_group), 0)))
  expect_true(all(5 + by_group >= 1 & 5 + by_group <= 10))
  cp_prior = prior_summary(fit) %>% dplyr::filter(.data$parameter == "cp_1_id")
  expect_equal(cp_prior$bounds, "[min(x) - cp_1, max(x) - cp_1]")
})


test_that("prior draws keep a single varying change point in range", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    sample = "prior",
    chains = 1, warmup = 100, iter = 200,
    seed = 42, quiet = TRUE
  )
  draws = as.matrix(.subset2(fit, "mcmc_prior"))
  deviations = draws[, grep("^cp_1_id\\[", colnames(draws)), drop = FALSE]
  locations = sweep(deviations, 1, draws[, "cp_1"], "+")
  expect_true(all(locations >= 1 & locations <= 10))
})


test_that("simulation samples adjacent varying change points in realized order", {
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1, 1 + (1 | id) ~ 1),
    cp_group_data,
    par_x = "x",
    sample = FALSE
  )

  set.seed(1)
  result = fit$simulate(
    fit,
    cp_group_data,
    cp_1 = 3,
    cp_2 = 6,
    cp_1_sd = 10,
    cp_2_sd = 10,
    Intercept_1 = 0,
    Intercept_2 = 0,
    Intercept_3 = 0,
    sigma_1 = 1,
    .type = "fitted"
  )
  simulation = attr(result, "simulated")
  cp_1 = 3 + vapply(split(simulation$cp_1_id, cp_group_data$id), unique, numeric(1))
  cp_2 = 6 + vapply(split(simulation$cp_2_id, cp_group_data$id), unique, numeric(1))
  expect_true(all(cp_1 < cp_2))
  expect_true(all(cp_1 >= 1 & cp_2 <= 10))
})


test_that("sampled change points are checked for range and ordering", {
  cps = tibble::tibble(
    name = c("cp_1", "cp_2"),
    varying = c(TRUE, TRUE),
    group_name = c("cp_1_id", "cp_2_id")
  )
  draws = cbind(
    cp_1 = c(3, 4), cp_2 = c(7, 8),
    `cp_1_id[A]` = c(1, 1), `cp_2_id[A]` = c(-1, -1),
    `cp_1_id[B]` = c(-1, -1), `cp_2_id[B]` = c(1, 1)
  )
  samples = coda::mcmc.list(coda::mcmc(draws))

  expect_no_error(assert_ordered_cp_draws(samples, cps, 0:10))

  unordered_population = samples
  unordered_population[[1]][1, "cp_2"] = 2
  expect_error(
    assert_ordered_cp_draws(unordered_population, cps, 0:10),
    "population-level.*strictly ordered"
  )

  unordered_group = samples
  unordered_group[[1]][1, "cp_1_id[A]"] = 5
  expect_error(
    assert_ordered_cp_draws(unordered_group, cps, 0:10),
    "group-level.*strictly ordered"
  )

  out_of_range = samples
  out_of_range[[1]][1, "cp_1"] = -1
  expect_error(
    assert_ordered_cp_draws(out_of_range, cps, 0:10),
    "population-level.*observed range"
  )

  out_of_range_group = samples
  out_of_range_group[[1]][1, "cp_1_id[A]"] = -4
  expect_error(
    assert_ordered_cp_draws(out_of_range_group, cps, 0:10),
    "group-level.*observed range"
  )
})
