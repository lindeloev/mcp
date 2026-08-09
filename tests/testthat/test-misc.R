# Link functions
test_that("Link functions", {
  expect_equal(mcp:::ilogit(1), 0.73105858)
  expect_equal(mcp:::logit(0.7), 0.84729786)
  expect_equal(mcp:::phi(1), 0.84134475)
  expect_equal(mcp:::probit(0.7), 0.52440051)
  expect_false(any(c("ilogit", "logit", "phi", "probit") %in% getNamespaceExports("mcp")))
})

test_that("families", {
  expect_true(is.mcpfamily(bernoulli()))
  expect_false(is.mcpfamily(gaussian()))
})


test_that("model metadata uses aligned table names and varying selectors", {
  data = data.frame(
    x = 1:6,
    y = c(2, 4, 3, 8, 7, 9),
    id = rep(c("a", "b"), each = 3)
  )
  fit = mcp(list(y ~ 1, (1 | id) ~ 1), data, par_x = "x", sample = FALSE)
  tables = get_fit_model_tables(fit)

  expect_named(
    tables,
    c("segments", "cps", "predictors", "group_effects", "pars")
  )
  expect_named(
    tables$group_effects,
    c(
      "population_name", "name", "part", "group_col", "segment", "dpar",
      "sd_name", "par_type", "matrix_name", "display_name", "order",
      "x_factor", "matrix_col", "matrix_data", "next_segment", "correlated"
    )
  )
  expect_true(all(tables$pars$part %in% c("cp", "predictor")))
  expect_true(all(tables$pars$scope %in% c("population", "group")))

  expect_equal(unpack_varying(fit, pars = "cp")$pars, fit$pars$group)
  expect_null(unpack_varying(fit, pars = "predictor")$pars)
  expect_equal(
    unpack_varying(fit, pars = fit$pars$group)$pars,
    fit$pars$group
  )
  expect_error(unpack_varying(fit, pars = "unknown"), "Unknown `varying`")
})


test_that("legacy fit detection blocks pre-0.4.0 mcpfit objects", {
  data = data.frame(
    x = 1:6,
    y = c(2, 4, 3, 8, 7, 9),
    id = rep(c("a", "b"), each = 3)
  )
  fit = mcp(list(y ~ 1, (1 | id) ~ 1), data, par_x = "x", sample = FALSE)
  legacy = fit
  legacy$family = gaussian()

  expect_error(
    get_fit_model_tables(legacy),
    "created with an older version of mcp (< 0.4.0)",
    fixed = TRUE
  )
})


test_that("helpful deprecation detections work for old code idioms", {
  data = data.frame(x = 1:6, y = c(2, 4, 3, 8, 7, 9))
  fit = mcp(list(y ~ 1 + x), data, sample = FALSE)

  # data = NULL or missing in mcp()
  expect_error(
    mcp(list(y ~ 1), data = NULL),
    "`data` is required in mcp() since mcp v0.4.0",
    fixed = TRUE
  )

  # plot/fitted with which_y
  expect_warning(
    plot(demo_fit, which_y = "sigma"),
    "`which_y` was deprecated in mcp v0.4.0. Use `plot_dpar()` instead.",
    fixed = TRUE
  )

  # fit$simulate() signature change
  expect_error(
    fit$simulate(1:5),
    "breaking changes",
    fixed = FALSE
  )

  # ex$fit warning
  expect_warning(
    { val = fit$fit },
    "`mcp_example()` now returns an `mcpfit` object directly",
    fixed = TRUE
  )
  expect_s3_class(val, "mcpfit")
})


# Code
test_that("formula tools", {
  expect_equal(sd_to_prec("dnorm(5, 1)"), "dnorm(5, 1/(1)^2) ")
  expect_equal(sd_to_prec("dt(3,2,1)"), "dt(3, 1/(2)^2, 1) ")
})

test_that("priors are resolved without changing their parameterization", {
  data = data.frame(x = 1:6, y = c(2, 4, 3, 8, 7, 9))
  default_fit = mcp(list(y ~ 1 + x, ~ 1 + x), data, sample = FALSE)
  expect_equal(
    unclass(default_fit$prior),
    list(
      cp_1 = "dunif(1, 6)",
      Intercept_1 = "dt(5.5, 3.7, 3)",
      x_1 = "dt(0, 1.48, 3)",
      Intercept_2 = "dt(5.5, 3.7, 3)",
      x_2 = "dt(0, 1.48, 3)",
      sigma_1 = "dt(0, 3.7, 3) T(0, )"
    )
  )

  fit = mcp(
    list(y ~ 1 + x, ~ 1 + x),
    data,
    prior = list(
      Intercept_1 = "dt(median(y), mad(y), 3)T(, max(y))",
      x_1 = 5,
      Intercept_2 = "Intercept_1",
      x_2 = "x_1 / (max(x) - min(x))"
    ),
    sample = FALSE
  )

  expect_equal(fit$prior$Intercept_1, "dt(5.5, 3.7065, 3) T(, 9)")
  expect_equal(fit$prior$x_2, "x_1/5")
  expect_false(grepl("MINX|MAXX|N_CP|LINKY", fit$jags_code))
  expect_match(fit$jags_code, "# User-specified prior", fixed = TRUE)

  compact = prior_summary(fit)
  verbose = prior_summary(fit, verbose = TRUE)
  expect_named(compact, c("parameter", "segment", "dpar", "prior", "bounds"))
  expect_named(
    verbose,
    c("parameter", "segment", "dpar", "prior", "bounds", "rule", "description", "source", "kind")
  )
  expect_equal(
    verbose$kind[match(c("Intercept_1", "x_1", "Intercept_2", "x_2"), verbose$parameter)],
    c("distribution", "constant", "alias", "expression")
  )
  expect_equal(verbose$source[verbose$parameter == "Intercept_1"], "user")

  legacy_fit = NULL
  expect_warning(
    {
      legacy_fit = mcp(
        list(y ~ 1, ~1), data,
        par_x = "x", prior = list(cp_1 = "dunif(MINX, MAXX)"), sample = FALSE
      )
    },
    "Deprecated prior data constant"
  )
  expect_equal(legacy_fit$prior$cp_1, "dunif(1, 6)")

  legacy_fit$.internal$prior_table = NULL
  legacy_fit$prior$cp_1 = "dunif(MINX, MAXX)"
  expect_equal(prior_summary(legacy_fit)$prior[1], "uniform(min = 1, max = 6)")
})


test_that("Gaussian defaults use coherent response and link scales", {
  log_data = data.frame(
    x = 1:6,
    y = c(-10, 80, 90, 100, 110, 120)
  )
  log_fit = mcp(
    list(y ~ 1 + x), log_data,
    family = gaussian(link = "log"), sample = FALSE
  )

  # Non-positive responses are valid: the log link applies to mu, not y.
  expect_equal(log_fit$prior$Intercept_1, "dt(0, 2.5, 3)")
  expect_equal(log_fit$prior$x_1, "dt(0, 0.5, 3)")
  expect_equal(log_fit$prior$sigma_1, "dt(0, 22.2, 3) T(0, )")
  expect_false(grepl("log\\(y\\)", log_fit$jags_code))

  rules = prior_summary(log_fit, verbose = TRUE)
  sigma_rule = rules$rule[rules$parameter == "sigma_1"]
  expect_match(sigma_rule, "mad(y)", fixed = TRUE)
  expect_false(grepl("mad(log(y))", sigma_rule, fixed = TRUE))

  wide_data = data.frame(x = 1:6, y = exp(c(-10, -5, 0, 5, 10, 15)))
  wide_fit = mcp(
    list(y ~ 1 + x), wide_data,
    family = gaussian(link = "log"), sample = FALSE
  )
  expect_equal(wide_fit$prior$Intercept_1, "dt(2.5, 11.1, 3)")
  expect_equal(wide_fit$prior$x_1, "dt(0, 2.22, 3)")
  expect_equal(wide_fit$prior$sigma_1, "dt(0, 110.8, 3) T(0, )")

  zero_data = data.frame(x = 1:6, y = 0:5)
  zero_fit = mcp(
    list(y ~ 1 + x), zero_data,
    family = gaussian(link = "log"), sample = FALSE
  )
  expect_equal(zero_fit$prior$Intercept_1, "dt(0.9, 2.5, 3)")

  small_data = data.frame(x = 1:4, y = c(0.01, 0.02, 0.03, 0.04))
  small_fit = mcp(list(y ~ 1 + x), small_data, sample = FALSE)
  expect_equal(small_fit$prior$Intercept_1, "dt(0, 2.5, 3)")
  expect_equal(small_fit$prior$sigma_1, "dt(0, 2.5, 3) T(0, )")
})

test_that("parameter-name collisions give a useful error", {
  data = data.frame(
    y = 1:6,
    x = 1:6,
    a = c(0, 0, 0, 1, 1, 1),
    b = c(0, 1, 2, 0, 1, 2),
    ab = c(0, 1, 0, 1, 0, 1)
  )

  expect_error(
    mcp(list(y ~ a:b + ab), data, par_x = "x", sample = FALSE),
    "`ab_1`: `ab` (mu, segment 1) and `a:b` (mu, segment 1)",
    fixed = TRUE
  )
})


########################
# MCPFIT CLASS-METHODS #
########################
# Test on new fit
demo_settings = mcp_example("demo", sample = FALSE, plot = FALSE)
demo_fit_iter = 50  # only niterations()/nchains() metadata is checked below, not recovery
demo_fit2 = suppressWarnings(mcp(demo_settings$model, demo_settings$data, adapt = 50, iter = demo_fit_iter, warn = FALSE, quiet = TRUE))

test_that("binomial example can be constructed without sampling", {
  fit = mcp_example("binomial", sample = FALSE, plot = FALSE)
  expect_s3_class(fit, "mcpfit")
  expect_true(all(fit$data$y <= fit$data$N))
})

test_that("group_mu example contains independent factor effects", {
  fit = mcp_example("group_mu", sample = FALSE, plot = FALSE)
  expect_s3_class(fit, "mcpfit")
  expect_equal(fit$pars$cp, "cp_1")
  expect_equal(
    fit$pars$group,
    c("Intercept_1_id", "conditionB_1_id")
  )
  expect_match(fit$call, "condition || id", fixed = TRUE)
  expect_equal(length(unique(fit$data$id)), 9)
  simulated = attr(fit$data$y, "simulated")
  expect_equal(simulated$Intercept_1_id_sd, 2)
  expect_equal(simulated$conditionB_1_id_sd, 2)
  expect_length(unique(simulated$Intercept_1_id), 9)
})

test_that("cores is deprecated and ignored", {
  model = list(y ~ 1)
  data = data.frame(x = 1:3, y = 1:3)

  expect_warning(
    mcp(model, data, par_x = "x", sample = FALSE, cores = 2),
    "Setting `cores` above one no longer enables parallel processing",
    class = "lifecycle_warning_deprecated"
  )
})

test_that("Simple mcpfit methods", {
  expect_equal(ndraws(demo_fit2), demo_fit_iter * 3)
  expect_equal(niterations(demo_fit2), demo_fit_iter)
  expect_equal(nchains(demo_fit2), 3)

  expect_true(is.mcpfit(demo_fit))
  expect_false(is.mcpfit(mtcars))
})

test_that("posterior draws accessor preserves the stored chains", {
  expect_identical(as_draws, posterior::as_draws)
  expect_identical(as_draws_df, posterior::as_draws_df)
  expect_identical(as_draws_array, posterior::as_draws_array)
  expect_identical(as_draws_matrix, posterior::as_draws_matrix)
  expect_identical(as_draws_rvars, posterior::as_draws_rvars)

  draws_array = as_draws_array(demo_fit)
  draws_df = as_draws_df(demo_fit)
  draws_matrix = as_draws_matrix(demo_fit)
  mcmc = coda::as.mcmc(demo_fit)

  expect_s3_class(draws_array, "draws_array")
  expect_s3_class(draws_df, "draws_df")
  expect_s3_class(draws_matrix, "draws_matrix")
  expect_s3_class(mcmc, "mcmc.list")

  raw = .subset2(demo_fit, "mcmc_post")
  expect_equal(
    dim(draws_array),
    c(nrow(raw[[1]]), length(raw), ncol(raw[[1]]))
  )
  expect_equal(posterior::variables(draws_array), colnames(raw[[1]]))

  # Accessing mcmc_post directly should soft-deprecate
  lifecycle::expect_deprecated({
    val = demo_fit$mcmc_post
  })
  expect_s3_class(val, "mcmc.list")

  if (requireNamespace("tidybayes", quietly = TRUE)) {
    td = tidybayes::tidy_draws(demo_fit)
    expect_s3_class(td, "tbl_df")
    expect_true("Intercept_1" %in% names(td))

    sd_df = tidybayes::spread_draws(demo_fit, Intercept_1, cp_1)
    expect_s3_class(sd_df, "tbl_df")
    expect_true(all(c("Intercept_1", "cp_1") %in% names(sd_df)))
  }
})

test_that("summaries use central intervals and posterior diagnostics", {
  width = 0.8
  result = fixef(demo_fit, width = width)
  parameter = result$name[[1]]
  raw = .subset2(demo_fit, "mcmc_post")
  values = unlist(lapply(raw, function(chain) chain[, parameter]))
  parameter_matrix = posterior::extract_variable_matrix(
    posterior_draws(demo_fit),
    variable = parameter
  )

  expect_equal(result$lower[[1]], unname(quantile(values, 0.1)))
  expect_equal(result$upper[[1]], unname(quantile(values, 0.9)))
  expect_equal(result$Rhat[[1]], posterior::rhat(parameter_matrix))
  expect_equal(result$ess_bulk[[1]], round(posterior::ess_bulk(parameter_matrix)))
  expect_equal(result$ess_tail[[1]], round(posterior::ess_tail(parameter_matrix)))
  expect_true(all(c("Rhat", "ess_bulk", "ess_tail") %in% names(result)))
  expect_false("n.eff" %in% names(result))
})

test_that("PPC and LOO draws stay aligned", {
  expect_error(
    pp_check(demo_fit, ndraws = 2, not_a_bayesplot_argument = TRUE),
    "not_a_bayesplot_argument",
    fixed = TRUE
  )
  expect_s3_class(
    pp_check(demo_fit, type = "ribbon", ndraws = NULL, alpha = 0.2),
    "ggplot"
  )
  expect_error(
    pp_check(demo_fit, type = "loo_intervals", prior = TRUE, ndraws = 2),
    "LOO predictive checks require posterior draws",
    fixed = TRUE
  )

  fit = demo_fit
  fit$data[[fit$pars$y]][2] = NA_real_
  fit$data$facet = factor(rep(1:2, length.out = nrow(fit$data)))
  fit$loglik = NULL
  fit$loo = NULL

  expect_s3_class(pp_check(fit, ndraws = 5), "ggplot")
  fit = add_loglik(fit, ndraws = 10)
  expect_equal(dim(fit$loglik), c(10, nrow(fit$data) - 1))
  expect_false(anyNA(fit$loglik))

  loo_result = suppressWarnings(loo(fit, ndraws = 10, save_psis = TRUE))
  expect_equal(dim(loo_result$psis_object), dim(fit$loglik))
  expect_equal(attr(loo_result, "mcp_settings")$ndraws, 10L)
  loo_changed = suppressWarnings(loo(fit, ndraws = 10, varying = FALSE, arma = FALSE))
  expect_false(attr(loo_changed, "mcp_settings")$varying)
  expect_false(attr(loo_changed, "mcp_settings")$arma)

  expect_s3_class(
    suppressWarnings(suppressMessages(pp_check(
      fit,
      type = "loo_intervals",
      facet_by = "facet",
      ndraws = 5
    ))),
    "patchwork"
  )
})


test_that("nsamples is a soft-deprecated alias for ndraws", {
  lifecycle::expect_deprecated({
    result = fitted(demo_fit, nsamples = 2, summary = FALSE)
  })
  expect_equal(length(unique(result$.draw)), 2)
  expect_error(
    suppressWarnings(fitted(demo_fit, ndraws = 2, nsamples = 2)),
    "Use only one of `ndraws` and deprecated `nsamples`.",
    fixed = TRUE
  )
})

# hypothesis()
test_that("hypothesis()", {
  # Use a draw-derived threshold so the hypothesis is neither rare nor certain.
  raw = .subset2(demo_fit2, "mcmc_post")
  cp_draws = unlist(lapply(raw, function(chain) chain[, "cp_1"]))
  threshold = unname(stats::quantile(cp_draws, 0.25))
  threshold_text = format(threshold, digits = 16)
  directional = paste0("cp_1 > ", threshold_text)

  expect_error(
    hypothesis(demo_fit2, directional),
    "Directional Bayes factors require both prior and posterior samples",
    fixed = TRUE
  )

  fit_asymmetric = demo_fit2
  fit_asymmetric$mcmc_prior = .subset2(demo_fit2, "mcmc_post")

  # Force the prior probability above the threshold to 0.25. This makes the
  # prior odds differ from the posterior odds.
  mcmc_prior = .subset2(fit_asymmetric, "mcmc_prior")
  for (chain in seq_along(mcmc_prior)) {
    n_draws = nrow(mcmc_prior[[chain]])
    mcmc_prior[[chain]][, "cp_1"] = c(
      rep(threshold + 1, floor(n_draws / 4)),
      rep(threshold - 1, n_draws - floor(n_draws / 4))
    )
  }
  fit_asymmetric$mcmc_prior = mcmc_prior

  actual_directional = expect_no_warning(hypothesis(fit_asymmetric, directional))
  p_post = mean(cp_draws > threshold)
  prior_draws = unlist(lapply(.subset2(fit_asymmetric, "mcmc_prior"), function(chain) chain[, "cp_1"]))
  p_prior = mean(prior_draws > threshold)
  expected_BF = (p_post / (1 - p_post)) / (p_prior / (1 - p_prior))
  effect_draws = cp_draws - threshold
  expect_equal(actual_directional$lower, unname(quantile(effect_draws, 0.025)))
  expect_equal(actual_directional$upper, unname(quantile(effect_draws, 0.975)))
  expect_equal(actual_directional$p, p_post)
  expect_equal(actual_directional$BF, expected_BF)

  # Identical prior and posterior draws must give BF = 1, also for intervals and Savage-Dickey equality.
  fit_same = demo_fit2
  fit_same$mcmc_prior = .subset2(demo_fit2, "mcmc_post")
  bounds = stats::quantile(cp_draws, c(0.2, 0.8))
  interval = paste0(
    "cp_1 > ", format(bounds[[1]], digits = 16),
    " & cp_1 < ", format(bounds[[2]], digits = 16)
  )
  actual_interval = hypothesis(fit_same, interval)
  expect_equal(actual_interval$p, mean(cp_draws > bounds[[1]] & cp_draws < bounds[[2]]))
  expect_equal(actual_interval$BF, 1)

  # Savage-Dickey point equality test (requires prior)
  mid_val = format(mean(cp_draws), digits = 16)
  equality_expr = paste0("cp_1 = ", mid_val)
  actual_equality = hypothesis(fit_same, equality_expr)
  expect_s3_class(actual_equality, "data.frame")
  expect_equal(actual_equality$hypothesis, paste0("cp_1 - ", mid_val, " = 0"))
  expect_false(is.na(actual_equality$p))
  expect_false(is.na(actual_equality$BF))
  expect_equal(actual_equality$BF, 1, tolerance = 1e-3)
})


test_that("warn parameter and summary convergence footer work as expected", {
  data = data.frame(x = 1:20, y = rnorm(20))
  model = list(y ~ 1, ~ 1)

  # 1. warn = FALSE suppresses sampling convergence warning
  fit_nowarn = suppressWarnings(mcp(model, data, par_x = "x", iter = 50, adapt = 50, warn = FALSE, quiet = TRUE))
  expect_equal(fit_nowarn$.internal$warn, FALSE)

  # 2. summary() includes warning footer if convergence is poor
  # Force high Rhat by modifying posterior draws
  fit_bad = fit_nowarn
  raw_post = .subset2(fit_bad, "mcmc_post")
  raw_post[[1]][, "cp_1"] = 5
  raw_post[[2]][, "cp_1"] = 15
  raw_post[[3]][, "cp_1"] = 25
  fit_bad$mcmc_post = raw_post

  sum_out = capture.output(summary(fit_bad))
  expect_true(any(grepl("Warning: .* parameter.* show.* poor convergence", sum_out)))

  # 3. mcp_example defaults to warn = FALSE
  ex_fit = mcp_example("intercepts", sample = FALSE)
  expect_s3_class(ex_fit, "mcpfit")
})
