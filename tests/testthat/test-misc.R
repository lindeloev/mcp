# Link functions
test_that("Link functions", {
  expect_equal(ilogit(1), 0.73105858)
  expect_equal(logit(0.7), 0.84729786)
  expect_equal(phi(1), 0.84134475)
  expect_equal(probit(0.7), 0.52440051)
})

test_that("families", {
  expect_true(is.mcpfamily(bernoulli()))
  expect_false(is.mcpfamily(gaussian()))
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
    legacy_fit <- mcp(
      list(y ~ 1, ~ 1), data,
      par_x = "x", prior = list(cp_1 = "dunif(MINX, MAXX)"), sample = FALSE
    ),
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
  expect_equal(log_fit$prior$Intercept_1, "dt(4.6, 2.5, 3)")
  expect_equal(log_fit$prior$x_1, "dt(0, 0.5, 3)")
  expect_equal(log_fit$prior$sigma_1, "dt(0, 22.2, 3) T(0, )")
  expect_false(grepl("log\\(y\\)", log_fit$jags_code))

  rules = prior_summary(log_fit, verbose = TRUE)
  sigma_rule = rules$rule[rules$parameter == "sigma_1"]
  expect_match(sigma_rule, "mad(y)", fixed = TRUE)
  expect_false(grepl("mad(log(y))", sigma_rule, fixed = TRUE))

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
demo_settings = mcp_example("demo", sample = FALSE)
demo_fit2 = quiet_mcp(demo_settings$model, demo_settings$data, adapt = 2500, iter = 4000)

test_that("binomial example can be constructed without sampling", {
  fit = mcp_example("binomial", sample = FALSE)
  expect_s3_class(fit, "mcpfit")
  expect_true(all(fit$data$y <= fit$data$N))
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
  expect_equal(niterations(demo_fit2), 12000)
  expect_equal(nchains(demo_fit), 3)

  expect_true(is.mcpfit(demo_fit))
  expect_false(is.mcpfit(mtcars))

})

test_that("posterior draws accessor preserves the stored chains", {
  draws = posterior_draws(demo_fit)

  expect_s3_class(draws, "draws_array")
  expect_equal(
    dim(draws),
    c(
      nrow(demo_fit$mcmc_post[[1]]),
      length(demo_fit$mcmc_post),
      ncol(demo_fit$mcmc_post[[1]])
    )
  )
  expect_equal(posterior::variables(draws), colnames(demo_fit$mcmc_post[[1]]))
  expect_s3_class(demo_fit$mcmc_post, "mcmc.list")
})

test_that("summaries use central intervals and posterior diagnostics", {
  width = 0.8
  result = fixef(demo_fit, width = width)
  parameter = result$name[[1]]
  values = unlist(lapply(
    demo_fit$mcmc_post,
    function(chain) chain[, parameter]
  ))
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

  loo_result = suppressWarnings(loo(
    fit, ndraws = 10, save_psis = TRUE
  ))
  expect_equal(dim(loo_result$psis_object), dim(fit$loglik))
  expect_equal(attr(loo_result, "mcp_settings")$ndraws, 10L)
  loo_changed = suppressWarnings(loo(
    fit, ndraws = 10, varying = FALSE, arma = FALSE
  ))
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
  lifecycle::expect_deprecated(
    result <- fitted(demo_fit, nsamples = 2, summary = FALSE)
  )
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
  cp_draws = unlist(lapply(demo_fit2$mcmc_post, function(chain) chain[, "cp_1"]))
  threshold = unname(stats::quantile(cp_draws, 0.25))
  threshold_text = format(threshold, digits = 16)
  directional = paste0("cp_1 > ", threshold_text)

  expect_error(
    hypothesis(demo_fit2, directional),
    "Directional Bayes factors require both prior and posterior samples",
    fixed = TRUE
  )

  fit_asymmetric = demo_fit2
  fit_asymmetric$mcmc_prior = demo_fit2$mcmc_post

  # Force the prior probability above the threshold to 0.25. This makes the
  # prior odds differ from the posterior odds.
  for (chain in seq_along(fit_asymmetric$mcmc_prior)) {
    n_draws = nrow(fit_asymmetric$mcmc_prior[[chain]])
    fit_asymmetric$mcmc_prior[[chain]][, "cp_1"] = c(
      rep(threshold + 1, floor(n_draws / 4)),
      rep(threshold - 1, n_draws - floor(n_draws / 4))
    )
  }

  actual_directional = expect_no_warning(hypothesis(fit_asymmetric, directional))
  p_post = mean(cp_draws > threshold)
  prior_draws = unlist(lapply(fit_asymmetric$mcmc_prior, function(chain) chain[, "cp_1"]))
  p_prior = mean(prior_draws > threshold)
  expected_BF = (p_post / (1 - p_post)) / (p_prior / (1 - p_prior))
  effect_draws = cp_draws - threshold
  expect_equal(actual_directional$lower, unname(quantile(effect_draws, 0.025)))
  expect_equal(actual_directional$upper, unname(quantile(effect_draws, 0.975)))
  expect_equal(actual_directional$p, p_post)
  expect_equal(actual_directional$BF, expected_BF)

  # Identical prior and posterior draws must give BF = 1, also for intervals.
  fit_same = demo_fit2
  fit_same$mcmc_prior = demo_fit2$mcmc_post
  bounds = stats::quantile(cp_draws, c(0.2, 0.8))
  interval = paste0(
    "cp_1 > ", format(bounds[[1]], digits = 16),
    " & cp_1 < ", format(bounds[[2]], digits = 16)
  )
  actual_interval = hypothesis(fit_same, interval)
  expect_equal(actual_interval$p, mean(cp_draws > bounds[[1]] & cp_draws < bounds[[2]]))
  expect_equal(actual_interval$BF, 1)
})
