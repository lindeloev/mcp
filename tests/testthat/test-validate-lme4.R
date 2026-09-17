if (Sys.getenv("MCP_TEST_LEVEL") != "release") {
  testthat::skip("Time-consuming validation tests against reference implementations are only run when MCP_TEST_LEVEL='release'.")
}

test_that("Gaussian random-intercept model agrees with lme4::lmer()", {
  testthat::skip_if_not_installed("lme4")

  set.seed(42)
  N_groups = 10
  N_per_group = 25
  id = factor(rep(paste0("g", 1:N_groups), each = N_per_group), levels = paste0("g", 1:N_groups))
  x = rep(seq(0, 20, length.out = N_per_group), N_groups)
  id_effects = rnorm(N_groups, mean = 0, sd = 1.5)
  names(id_effects) = paste0("g", 1:N_groups)
  y = 2 + 0.5 * x + id_effects[as.character(id)] + rnorm(length(x), sd = 1.0)
  df = data.frame(x = x, y = y, id = id)

  # Fit with lmer
  fit_lmer = lme4::lmer(y ~ x + (1 | id), data = df, REML = FALSE)

  # Fit with mcp
  model = list(y ~ 1 + x + (1 | id))
  fit_mcp = mcp(model, df, family = gaussian(), warmup = 500, iter = 2000, chains = 2, seed = 42, diagnostics = FALSE, quiet = TRUE)

  # Extract mcp estimates
  params_mcp = fixef(fit_mcp)
  capture.output({ summary_mcp = summary(fit_mcp) })
  sigma_mcp = summary_mcp$mean[summary_mcp$variable == "sigma_1"]
  sd_mcp = summary_mcp$mean[summary_mcp$variable == "Intercept_1_id_sd"]

  # Extract lmer estimates
  params_lmer = lme4::fixef(fit_lmer)
  sigma_lmer = stats::sigma(fit_lmer)
  sd_lmer = as.data.frame(lme4::VarCorr(fit_lmer))$sdcor[1]

  # Compare fixed effects (within 0.1)
  expect_equal(params_mcp$mean[params_mcp$variable == "Intercept_1"], unname(params_lmer[1]), tolerance = 0.1)
  expect_equal(params_mcp$mean[params_mcp$variable == "x_1"], unname(params_lmer[2]), tolerance = 0.05)

  # Compare variance components
  expect_equal(sigma_mcp, sigma_lmer, tolerance = 0.05)
  expect_equal(sd_mcp, sd_lmer, tolerance = 0.25)

  # Compare random effect deviations: ranef() now follows factor levels directly
  ran_mcp = ranef(fit_mcp)
  expect_equal(ran_mcp$variable, paste0("Intercept_1_id[", levels(df$id), "]"))
  ranef_lmer = lme4::ranef(fit_lmer)$id[, 1]
  expect_gt(stats::cor(ranef_lmer, ran_mcp$mean), 0.99)
  expect_lt(max(abs(ranef_lmer - ran_mcp$mean)), 0.1)

  # Compare conditional log-likelihood (mcp log_lik is conditional on random deviations)
  cond_loglik_lmer = sum(stats::dnorm(df$y, stats::predict(fit_lmer), sigma_lmer, log = TRUE))
  loglik_mcp = mean(rowSums(log_lik(fit_mcp)))
  expect_equal(loglik_mcp, cond_loglik_lmer, tolerance = 0.05)
})


test_that("Gaussian uncorrelated random intercept and slope model agrees with lme4::lmer()", {
  testthat::skip_if_not_installed("lme4")

  set.seed(42)
  N_groups = 15
  N_per_group = 30
  id = factor(rep(paste0("g", 1:N_groups), each = N_per_group), levels = paste0("g", 1:N_groups))
  x = rep(seq(-5, 5, length.out = N_per_group), N_groups)
  u_0 = rnorm(N_groups, mean = 0, sd = 1.0)
  u_1 = rnorm(N_groups, mean = 0, sd = 0.2)
  names(u_0) = paste0("g", 1:N_groups)
  names(u_1) = paste0("g", 1:N_groups)
  y = 3 + 0.4 * x + u_0[as.character(id)] + u_1[as.character(id)] * x + rnorm(length(x), sd = 0.8)
  df = data.frame(x = x, y = y, id = id)

  fit_lmer = lme4::lmer(y ~ x + (1 + x || id), data = df, REML = FALSE)

  model = list(y ~ 1 + x + (1 + x || id))
  fit_mcp = mcp(model, df, family = gaussian(), warmup = 500, iter = 2000, chains = 2, seed = 42, diagnostics = FALSE, quiet = TRUE)

  params_mcp = fixef(fit_mcp)
  capture.output({ summary_mcp = summary(fit_mcp) })
  sigma_mcp = summary_mcp$mean[summary_mcp$variable == "sigma_1"]
  sd_int_mcp = summary_mcp$mean[summary_mcp$variable == "Intercept_1_id_sd"]
  sd_x_mcp = summary_mcp$mean[summary_mcp$variable == "x_1_id_sd"]

  params_lmer = lme4::fixef(fit_lmer)
  sigma_lmer = stats::sigma(fit_lmer)
  vc = as.data.frame(lme4::VarCorr(fit_lmer))
  sd_int_lmer = vc$sdcor[which(vc$var1 == "(Intercept)")]
  sd_x_lmer = vc$sdcor[which(vc$var1 == "x")]

  expect_equal(params_mcp$mean[params_mcp$variable == "Intercept_1"], unname(params_lmer[1]), tolerance = 0.1)
  expect_equal(params_mcp$mean[params_mcp$variable == "x_1"], unname(params_lmer[2]), tolerance = 0.15)
  expect_equal(sigma_mcp, sigma_lmer, tolerance = 0.05)
  expect_equal(sd_int_mcp, sd_int_lmer, tolerance = 0.2)
  expect_equal(sd_x_mcp, sd_x_lmer, tolerance = 0.15)
})


test_that("Gaussian fixed change point with group intercepts agrees with lme4::lmer()", {
  testthat::skip_if_not_installed("lme4")

  set.seed(42)
  N_groups = 12
  N_per_group = 25
  id = factor(rep(paste0("g", 1:N_groups), each = N_per_group), levels = paste0("g", 1:N_groups))
  x = rep(seq(0, 20, length.out = N_per_group), N_groups)
  id_effects = rnorm(N_groups, mean = 0, sd = 1.2)
  names(id_effects) = paste0("g", 1:N_groups)
  y = 2 + 0.5 * x + 1.0 * pmax(0, x - 10) + id_effects[as.character(id)] + rnorm(length(x), sd = 0.9)
  df = data.frame(x = x, y = y, id = id)

  fit_lmer = lme4::lmer(y ~ x + I(pmax(0, x - 10)) + (1 | id), data = df, REML = FALSE)

  model = list(
    y ~ 1 + x + (1 | id),
    ~ 0 + x + same((1 | id))
  )
  fit_mcp = mcp(model, df, family = gaussian(), prior = list(cp_1 = 10),
                warmup = 500, iter = 2000, chains = 2, seed = 42, diagnostics = FALSE, quiet = TRUE)

  params_mcp = fixef(fit_mcp)
  capture.output({ summary_mcp = summary(fit_mcp) })
  sigma_mcp = summary_mcp$mean[summary_mcp$variable == "sigma_1"]
  sd_mcp = summary_mcp$mean[summary_mcp$variable == "Intercept_1_id_sd"]

  params_lmer = lme4::fixef(fit_lmer)
  sigma_lmer = stats::sigma(fit_lmer)
  sd_lmer = as.data.frame(lme4::VarCorr(fit_lmer))$sdcor[1]

  expect_equal(params_mcp$mean[params_mcp$variable == "Intercept_1"], unname(params_lmer[1]), tolerance = 0.1)
  expect_equal(params_mcp$mean[params_mcp$variable == "x_1"], unname(params_lmer[2]), tolerance = 0.05)
  expect_equal(params_mcp$mean[params_mcp$variable == "x_2"], unname(params_lmer[2] + params_lmer[3]), tolerance = 0.05)
  expect_equal(sigma_mcp, sigma_lmer, tolerance = 0.05)
  expect_equal(sd_mcp, sd_lmer, tolerance = 0.35)
})


test_that("Group-level simulation against lme4::lmer()", {
  testthat::skip_if_not_installed("lme4")

  set.seed(42)
  N_groups = 10
  N_per_group = 30
  df = data.frame(
    x = rep(seq(0, 20, length.out = N_per_group), N_groups),
    id = factor(rep(paste0("g", 1:N_groups), each = N_per_group)),
    y = 0
  )

  fit_empty = mcp(list(y ~ 1 + x + (1 | id)), df, family = gaussian(), sample = FALSE)
  df$y = fit_empty$simulate(
    fit_empty, df,
    Intercept_1 = 2, x_1 = 0.5, Intercept_1_id_sd = 1.5, sigma_1 = 1.0,
    .type = "predict"
  )

  fit_lmer_sim = lme4::lmer(y ~ x + (1 | id), data = df, REML = FALSE)
  expect_equal(unname(lme4::fixef(fit_lmer_sim)[2]), 0.5, tolerance = 0.05)
  expect_equal(stats::sigma(fit_lmer_sim), 1.0, tolerance = 0.1)
  expect_equal(as.data.frame(lme4::VarCorr(fit_lmer_sim))$sdcor[1], 1.5, tolerance = 0.35)
})
