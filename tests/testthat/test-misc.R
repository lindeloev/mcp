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


test_that("mcp rejects non-syntactic data column names", {
  data = data.frame(`x value` = 1:3, y = 0, check.names = FALSE)

  expect_error(
    mcp(list(y ~ 1), data, par_x = "x value", sample = FALSE),
    "`data` has non-syntactic column name(s): `x value`",
    fixed = TRUE
  )
})

test_that("mcp rejects data names that collide with generated JAGS nodes", {
  data = data.frame(mu_ = 1:5, x = 1:5)
  expect_error(
    mcp(list(mu_ ~ 1), data, par_x = "x", sample = FALSE),
    "Data column name(s) collide with mcp's generated JAGS namespace: 'mu_'",
    fixed = TRUE
  )
})

test_that("mcp rejects generated parameter names that collide with change points", {
  data = data.frame(x = 1:5, cp = c(0, 1, 0, 1, 0), y = 1:5)
  expect_error(
    mcp(list(y ~ cp, ~ 1), data, par_x = "x", sample = FALSE),
    "Generated parameter name(s) collide in the JAGS namespace: 'cp_1'",
    fixed = TRUE
  )
})


test_that("mcpfit model accessors follow standard R conventions", {
  data = data.frame(x = 1:5, y = c(1, 2, NA, 4, 5), unused = 6:10)
  expect_message(
    {
      fit = mcp(list(y ~ 1 + x), data, sample = FALSE)
    },
    "NA values detected in 'y'",
    fixed = TRUE
  )

  expect_identical(names(fit)[1:5], c("model", "data", "prior", "family", "call"))
  expect_type(fit$call, "language")
  expect_identical(fit$call[[1]], quote(mcp))
  expect_identical(family(fit), fit$family)
  expect_equal(nobs(fit), 4)
  expect_identical(model.frame(fit), fit$data)
  expect_identical(formula(fit), fit$model)
  expect_identical(formula(fit, segment = 1), fit$model[[1]])
  expect_equal(mcp_columns(fit)$par_x, "x")
  expect_equal(mcp_columns(fit)$response, "y")

  population = mcp_pars(fit, scope = "population")$name
  fixed = mcp_pars(fit, scope = "population", role = "fixed_effect")$name
  values = matrix(seq_len(3 * length(population)), nrow = 3)
  colnames(values) = population
  prior_values = values^2
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(values))
  fit$mcmc_prior = coda::mcmc.list(coda::mcmc(prior_values))
  expect_equal(vcov(fit), stats::cov(values[, fixed, drop = FALSE]))
  expect_equal(vcov(fit, prior = TRUE), stats::cov(prior_values[, fixed, drop = FALSE]))
  expect_equal(vcov(fit, correlation = TRUE), stats::cor(values[, fixed, drop = FALSE]))
  expect_equal(vcov(fit, pars = "all"), stats::cov(values))
  expect_error(vcov(fit, pars = "nonexistent"))
  expected_intervals = t(vapply(
    population,
    function(parameter) stats::quantile(values[, parameter], c(0.025, 0.975), names = FALSE),
    numeric(2)
  ))
  colnames(expected_intervals) = c("2.5 %", "97.5 %")
  expect_equal(confint(fit), expected_intervals)

  prior_intervals = t(vapply(
    population,
    function(parameter) stats::quantile(prior_values[, parameter], c(0.025, 0.975), names = FALSE),
    numeric(2)
  ))
  colnames(prior_intervals) = c("2.5 %", "97.5 %")
  expect_equal(confint(fit, prior = TRUE), prior_intervals)

  fit$mcmc_post = NULL
  expect_error(vcov(fit), "Posterior requested but the posterior was not drawn", fixed = TRUE)
  expect_error(confint(fit), "Posterior requested but the posterior was not drawn", fixed = TRUE)
  expect_equal(vcov(fit, prior = TRUE), stats::cov(prior_values[, fixed, drop = FALSE]))
})


test_that("fixef() and vcov() select distributional-parameter coefficients", {
  data = data.frame(x = 1:5, y = 0)
  fit = mcp(list(y ~ 1 + sigma(1 + x)), data, par_x = "x", sample = FALSE)
  pars = mcp_pars(fit, scope = "population")
  values = matrix(seq_len(3 * nrow(pars)), nrow = 3, dimnames = list(NULL, pars$name))
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(values))
  sigma_names = pars$name[pars$dpar == "sigma" & pars$role == "dpar_effect"]

  expect_setequal(fixef(fit, dpar = "sigma")$variable, sigma_names)
  expect_equal(vcov(fit, dpar = "sigma"), stats::cov(values[, sigma_names, drop = FALSE]))
})


test_that("model metadata uses aligned table names and group-effect selectors", {
  data = data.frame(
    x = 1:6,
    y = c(2, 4, 3, 8, 7, 9),
    id = rep(c("a", "b"), each = 3)
  )
  fit = mcp(list(y ~ 1, (1 | id) ~ 1), data, par_x = "x", sample = FALSE)
  tables = get_fit_model_tables(fit)

  expect_named(
    tables,
    c("data_columns", "segments", "cps", "predictors", "group_effects", "parameters", "design_specs")
  )
  expect_named(tables$data_columns, c("par_x", "response", "series", "weights"))
  expect_named(
    tables$group_effects,
    c(
      "population_name", "name", "part", "group_col", "segment", "dpar",
      "sd_name", "par_type", "matrix_name", "display_name", "order",
      "x_factor", "design_id", "design_col", "matrix_col", "matrix_data",
      "next_segment", "correlated"
    )
  )
  expect_named(
    tables$parameters,
    c("name", "part", "scope", "role", "segment", "dpar", "order", "group_col", "population_name")
  )
  expect_true(all(tables$parameters$part %in% c("cp", "predictor")))
  expect_true(all(tables$parameters$scope %in% c("population", "group")))
  expect_true(all(tables$parameters$role %in% c("change_point", "fixed_effect", "dpar_effect", "arma", "group_sd", "group_deviation")))
  expect_equal(mcp_pars(fit), tables$parameters)
  group = mcp_pars(fit, scope = "group")$name

  expect_equal(unpack_group_effects(fit, pars = "cp")$pars, group)
  expect_null(unpack_group_effects(fit, pars = "predictor")$pars)
  expect_equal(
    unpack_group_effects(fit, pars = group)$pars,
    group
  )
  expect_error(unpack_group_effects(fit, pars = "unknown"), "Unknown group-effect")
})


test_that("mcp_pars() and mcp_columns() are family-neutral metadata accessors", {
  data = data.frame(x = 1:8, y = c(1, 0, 3, 2, 4, 1, 5, 3))
  fit = mcp(
    list(y ~ 1 + x + shape(1 + x)),
    data,
    family = negbinomial(),
    sample = FALSE
  )

  parameters = mcp_pars(fit)
  expect_true(all(c("shape_1", "shape_x_1") %in% parameters$name))
  expect_equal(parameters$dpar[parameters$name == "shape_x_1"], "shape")
  expect_equal(mcp_columns(fit)$par_x, "x")
  expect_equal(mcp_columns(fit)$response, "y")
  expect_false("pars" %in% names(fit))
  expect_warning(expect_null(fit$pars), "deprecated")
  expect_warning(expect_null(fit[["pars"]]), "deprecated")

  binomial_fit = mcp(
    list(y | trials(N) ~ 1),
    data.frame(x = 1:4, y = 0:3, N = 3),
    family = binomial(), par_x = "x", sample = FALSE
  )
  expect_equal(mcp_columns(binomial_fit)$trials, "N")

  # Runtime methods read the authoritative tables.
  expect_equal(
    as.numeric(fit$simulate(
      fit, data.frame(x = 1:2),
      Intercept_1 = 0, x_1 = 0, shape_1 = 0, shape_x_1 = 0,
      .type = "fitted"
    )),
    c(1, 1)
  )
  expect_equal(interpolate_newdata(fit, x_values = 1:2)$x, 1:2)
})


test_that("summary, fixef, and ranef select parameter roles consistently", {
  data = data.frame(
    x = 1:8,
    y = c(1, 2, 4, 3, 5, 6, 8, 7),
    id = rep(c("a", "b"), each = 4)
  )
  fit = mcp(
    list(y ~ 1 + ar(1) + (1 | id), (1 | id) ~ 1),
    data, par_x = "x", sample = FALSE
  )
  pars = mcp_pars(fit)
  expect_equal(pars$dpar[pars$name == "ar1_1"], "ar")
  expect_equal(pars$order[pars$name == "ar1_1"], 1L)
  expect_equal(pars$group_col[pars$name == "cp_1_id"], "id")
  expect_equal(pars$population_name[pars$name == "cp_1_id"], "cp_1")
  group_draws = unlist(lapply(pars$name[pars$scope == "group"], function(name) {
    paste0(name, "[", c("a", "b"), "]")
  }))
  draw_names = c(pars$name[pars$scope == "population"], group_draws)
  draws = matrix(seq_len(3 * length(draw_names)), nrow = 3, dimnames = list(NULL, draw_names))
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  capture.output({ summary_result = summary(fit) })
  expect_setequal(summary_result$variable, pars$name[pars$scope == "population"])
  expect_setequal(
    fixef(fit)$variable,
    pars$name[pars$scope == "population" & pars$role == "fixed_effect"]
  )
  expect_setequal(ranef(fit)$variable, group_draws)
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
    "install_github",
    fixed = TRUE
  )

  missing_role = fit
  missing_role$.internal$model_tables$parameters$role = NULL
  expect_error(
    get_fit_model_tables(missing_role),
    "install_github",
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
    "install_github",
    fixed = FALSE
  )

  # v0.3.4 parameter names in common post-fit code
  expect_error(
    hypothesis(fit, "int_1 > 0"),
    "use `Intercept_i`",
    fixed = TRUE
  )
  expect_error(
    plot_pars(fit, pars = "int_1"),
    "install_github",
    fixed = TRUE
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
test_that("legacy prior scales are translated to JAGS parameters", {
  expect_equal(prior_to_jags("dnorm(5, 2)"), "dnorm(5, 1/(2)^2) ")
  expect_equal(prior_to_jags("dt(3, 2, 1)"), "dt(3, 1/(2)^2, 1) ")
  expect_equal(prior_to_jags("dlnorm(3, 2)"), "dlnorm(3, 1/(2)^2) ")
  expect_equal(prior_to_jags("ddexp(3, 2)"), "ddexp(3, 1/(2)) ")
  expect_equal(prior_to_jags("dlogis(3, 2) T(0, )"), "dlogis(3, 1/(2)) T(0,)")
  expect_error(
    prior_to_jags("dcauchy(3, 2)"),
    "use `dt(location, scale, 1)`",
    fixed = TRUE
  )
  expect_warning(
    sd_to_prec("dnorm(5, 2)"),
    "deprecated"
  )
})

test_that("priors are resolved without changing their parameterization", {
  data = data.frame(x = 1:6, y = c(2, 4, 3, 8, 7, 9))
  default_fit = mcp(list(y ~ 1 + x, ~ 1 + x), data, sample = FALSE)
  expect_equal(
    unclass(default_fit$prior),
    list(
      cp_1 = "dunif(1, 6)",
      Intercept_1 = "dt(5.5, 3.7, 3)",
      x_1 = "dt(0, 0.74, 3)",
      Intercept_2 = "dt(5.5, 3.7, 3)",
      x_2 = "dt(0, 0.74, 3)",
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
  expect_equal(fit$.internal$prior_format, "jags_string_v1")
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
demo_settings = mcp_example("demo", sample = "none", plot = FALSE)
demo_fit_iter = 50  # only niterations()/nchains() metadata is checked below, not recovery
expect_warning({
  demo_fit2 = mcp(demo_settings$model, demo_settings$data, warmup = 50, iter = demo_fit_iter, diagnostics = FALSE, quiet = TRUE)
}, "Adaptation incomplete", fixed = TRUE)

test_that("binomial example can be constructed without sampling", {
  fit = mcp_example("binomial", sample = "none", plot = FALSE)
  expect_s3_class(fit, "mcpfit")
  expect_true(all(fit$data$y <= fit$data$N))
})

test_that("group_mu example contains independent factor effects", {
  fit = mcp_example("group_mu", sample = "none", plot = FALSE)
  expect_s3_class(fit, "mcpfit")
  expect_type(fit$call, "language")
  expect_s3_class(fit$example_code, "mcptext")
  expect_equal(get_fit_model_tables(fit)$cps$name, "cp_1")
  expect_equal(
    mcp_pars(fit, scope = "group")$name,
    c("Intercept_1_id", "conditionB_1_id")
  )
  expect_match(fit$example_code, "condition || id", fixed = TRUE)
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

    if (requireNamespace("rstantools", quietly = TRUE)) {
      newdata = demo_fit$data[1:3, , drop = FALSE]
      expect_equal(dim(rstantools::posterior_epred(demo_fit, newdata, draws = 2)), c(2, 3))
      expect_equal(dim(rstantools::posterior_predict(demo_fit, newdata, draws = 2)), c(2, 3))
      expect_equal(dim(rstantools::posterior_linpred(demo_fit, newdata = newdata, draws = 2)), c(2, 3))

      predicted_seed_1 = rstantools::posterior_predict(demo_fit, newdata, draws = 5, seed = 1)
      predicted_seed_1_again = rstantools::posterior_predict(demo_fit, newdata, draws = 5, seed = 1)
      predicted_seed_2 = rstantools::posterior_predict(demo_fit, newdata, draws = 5, seed = 2)
      expect_identical(predicted_seed_1, predicted_seed_1_again)
      expect_false(identical(predicted_seed_1, predicted_seed_2))

      epred_draws = tidybayes::add_epred_draws(newdata, demo_fit, ndraws = 2)
      predicted_draws = tidybayes::add_predicted_draws(newdata, demo_fit, ndraws = 2)
      linpred_draws = tidybayes::add_linpred_draws(newdata, demo_fit, ndraws = 2)
      expect_s3_class(epred_draws, "tbl_df")
      expect_true(".epred" %in% names(epred_draws))
      expect_true(".prediction" %in% names(predicted_draws))
      expect_true(".linpred" %in% names(linpred_draws))
      expect_equal(nrow(epred_draws), 6)
      expect_equal(nrow(predicted_draws), 6)
      expect_equal(nrow(linpred_draws), 6)

      binomial_data = data.frame(x = 1:3, N = c(4, 10, 6), y = 0)
      binomial_fit = mcp(
        list(y | trials(N) ~ 1),
        binomial_data,
        family = binomial(),
        par_x = "x",
        sample = FALSE
      )
      binomial_fit$mcmc_post = coda::mcmc.list(coda::mcmc(
        matrix(0, nrow = 20, ncol = 1, dimnames = list(NULL, "Intercept_1"))
      ))

      epred = rstantools::posterior_epred(binomial_fit, binomial_data, draws = 2)
      expect_equal(unname(epred), matrix(rep(binomial_data$N / 2, each = 2), nrow = 2))

      set.seed(42)
      predicted = rstantools::posterior_predict(binomial_fit, binomial_data, draws = 20)
      expect_true(all(predicted == floor(predicted)))
      expect_true(any(predicted > 1))

      tidy_epred = tidybayes::add_epred_draws(binomial_data, binomial_fit, ndraws = 2)
      expect_equal(sort(unique(tidy_epred$.epred)), sort(binomial_data$N / 2))
    }
  }
})

test_that("summaries use central intervals and posterior diagnostics", {
  width = 0.8
  result = fixef(demo_fit, width = width)
  printed = capture.output(summary(demo_fit, width = width))
  parameter = result$name[[1]]
  raw = .subset2(demo_fit, "mcmc_post")
  values = unlist(lapply(raw, function(chain) chain[, parameter]))
  parameter_matrix = posterior::extract_variable_matrix(
    posterior_draws(demo_fit),
    variable = parameter
  )

  expect_equal(result$lower[[1]], unname(quantile(values, 0.1)))
  expect_equal(result$upper[[1]], unname(quantile(values, 0.9)))
  expect_equal(result$sd[[1]], stats::sd(values))
  expect_equal(result$rhat[[1]], posterior::rhat(parameter_matrix))
  expect_equal(result$ess_bulk[[1]], round(posterior::ess_bulk(parameter_matrix)))
  expect_equal(result$ess_tail[[1]], round(posterior::ess_tail(parameter_matrix)))
  expect_true(all(c("Change point parameters:", "Population-level parameters:") %in% printed))
  expect_true(any(grepl("\\b[0-9]+\\.[0-9]{2}\\b", printed)))
  expect_true(any(grepl(" 1\\.0[01] ", printed)))
  expect_true("sd" %in% names(result))
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(result)))
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
  fit$data[[mcp_columns(fit)$response]][2] = NA_real_
  fit$data$facet = factor(rep(1:2, length.out = nrow(fit$data)))

  expect_s3_class(pp_check(fit, ndraws = 5), "ggplot")
  loglik = log_lik(fit, ndraws = 10)
  expect_equal(dim(loglik), c(10, nrow(fit$data) - 1))
  expect_false(anyNA(loglik))

  loo_result = suppressWarnings(loo(fit, ndraws = 10, save_psis = TRUE))
  expect_equal(dim(loo_result$psis_object), dim(loglik))
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

test_that("legacy mcpfits explain how to reproduce v0.3.4", {
  legacy_fit = demo_fit
  legacy_fit$.internal$model_tables = NULL

  expect_true(is_legacy_mcpfit(legacy_fit))
  expect_error(mcp_pars(legacy_fit), "install_github", fixed = TRUE)
  expect_error(summary(legacy_fit), "v0.4.0 bug fixes", fixed = TRUE)
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
    "Directional Bayes factors require both prior and posterior draws",
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

  prior_directional = hypothesis(fit_asymmetric, directional, prior = TRUE)
  expect_equal(prior_directional$p, p_prior)
  expect_true(is.na(prior_directional$BF))

  fit_prior_only = fit_asymmetric
  fit_prior_only$mcmc_post = NULL
  expect_equal(hypothesis(fit_prior_only, directional, prior = TRUE)$p, p_prior)

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
  expect_warning(
    actual_equality <- hypothesis(fit_same, equality_expr),
    "Savage-Dickey Bayes factor was computed using default prior(s) for `cp_1`",
    fixed = TRUE
  )
  expect_s3_class(actual_equality, "data.frame")
  expect_equal(actual_equality$hypothesis, paste0("cp_1 - ", mid_val, " = 0"))
  expect_true(is.na(actual_equality$p))
  expect_false(is.na(actual_equality$BF))
  expect_equal(actual_equality$BF, 1, tolerance = 1e-3)

  # With user-specified prior, no default prior warning is emitted
  fit_user_prior = fit_same
  fit_user_prior$.internal$prior_table$source[fit_user_prior$.internal$prior_table$parameter == "cp_1"] = "user"
  expect_no_warning(hypothesis(fit_user_prior, equality_expr))

  prior_equality = hypothesis(fit_same, equality_expr, prior = TRUE)
  expect_true(is.na(prior_equality$BF))

  tail_val = format(max(cp_draws) + stats::sd(cp_draws), digits = 16)
  expect_warning(
    hypothesis(fit_user_prior, paste0("cp_1 = ", tail_val)),
    "tested value is in a sparse tail",
    fixed = TRUE
  )
})


test_that("Savage-Dickey hypotheses are affine", {
  parameters = c("x", "y", "group[x]")

  expect_true(validate_savage_dickey_expression("x = 1", parameters))
  expect_true(validate_savage_dickey_expression("2 * x - y = 3", parameters))
  expect_true(validate_savage_dickey_expression("`group[x]` - x = 0", parameters))

  nonlinear = c(
    "x / y = 1",
    "x * y = 1",
    "x^2 = 1",
    "exp(x) = 1",
    "exp(y) * (x - 1) = 0"
  )
  for (expression in nonlinear) {
    expect_error(
      validate_savage_dickey_expression(expression, parameters),
      "named scalar parameter or affine contrast",
      fixed = TRUE
    )
  }
})


test_that("Savage-Dickey density is evaluated directly", {
  x = seq(-3, 3, length.out = 101)
  bandwidth = stats::bw.SJ(x)

  expect_equal(
    get_density(x, 0),
    mean(stats::dnorm(0, mean = x, sd = bandwidth))
  )
  expect_gte(get_density(x, 10), 0)
  expect_false(is_sparse_tail(x, 0))
  expect_true(is_sparse_tail(x, 3))
})


test_that("diagnostic settings control fit warnings and summary footers", {
  defaults = resolve_diagnostics()
  expect_equal(defaults, list(rhat = 1.01, ess_bulk = 400, ess_tail = 400, ar = 0.1, ma = 0.1))
  custom = resolve_diagnostics(list(ess_bulk = 800, ar = NULL))
  expect_equal(custom$ess_bulk, 800)
  expect_null(custom$ar)
  expect_equal(custom$rhat, defaults$rhat)
  expect_true(all(vapply(resolve_diagnostics(FALSE), is.null, logical(1))))
  expect_error(resolve_diagnostics(list(unknown = 1)), "Unknown diagnostic")

  data = data.frame(x = 1:20, y = rnorm(20))
  model = list(y ~ 1, ~ 1)

  # diagnostics = FALSE suppresses fit warnings and the inherited footer.
  expect_warning({
    fit_nowarn = mcp(model, data, par_x = "x", iter = 50, warmup = 50, diagnostics = FALSE, quiet = TRUE)
  }, "Adaptation incomplete", fixed = TRUE)
  expect_true(all(vapply(fit_nowarn$.internal$diagnostics, is.null, logical(1))))

  # Force high rhat by modifying posterior draws.
  fit_bad = fit_nowarn
  raw_post = .subset2(fit_bad, "mcmc_post")
  raw_post[[1]][, "cp_1"] = 5
  raw_post[[2]][, "cp_1"] = 15
  raw_post[[3]][, "cp_1"] = 25
  fit_bad$mcmc_post = raw_post

  inherited = capture.output(summary(fit_bad))
  expect_false(any(grepl("poor convergence", inherited)))

  strict = list(rhat = 1.01, ess_bulk = NULL, ess_tail = NULL)
  strict_out = capture.output(summary(fit_bad, diagnostics = strict))
  expect_true(any(grepl("poor convergence.*rhat > 1.01", strict_out)))

  convergence_off = list(rhat = NULL, ess_bulk = NULL, ess_tail = NULL)
  diagnostics_off = capture.output(summary(fit_bad, diagnostics = convergence_off))
  expect_false(any(grepl("poor convergence", diagnostics_off)))

  expect_warning(warn_nonconvergence(raw_post, strict), "rhat > 1.01")
  expect_no_warning(warn_nonconvergence(raw_post, convergence_off))

  ex_fit = mcp_example("intercepts", sample = "none")
  expect_s3_class(ex_fit, "mcpfit")
})

test_that("bernoulli accepts logical TRUE/FALSE responses and converts to 0/1", {
  data_logical = data.frame(
    x = 1:6,
    y = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )
  fit = mcp(list(y ~ 1 + x), data = data_logical, family = bernoulli(), sample = FALSE)
  expect_s3_class(fit, "mcpfit")
  expect_equal(fit$data$y, c(1, 0, 1, 0, 1, 0))
})

test_that("posterior_linpred evaluates binomial models on probability scale when transform = TRUE", {
  data_bin = data.frame(
    x = c(1, 2, 3),
    y = c(2, 3, 4),
    N = c(5, 5, 5)
  )
  fit = mcp(list(y | trials(N) ~ 1 + x), data = data_bin, family = binomial(), sample = FALSE)
  
  # Mock 1 posterior draw where Intercept_1 = 0, x_1 = 0 -> logit(p) = 0 -> p = 0.5
  population = mcp_pars(fit, scope = "population")$name
  draws = matrix(rep(0, length(population)), nrow = 1)
  colnames(draws) = population
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  linpred_link = rstantools::posterior_linpred(fit, transform = FALSE)
  linpred_prob = rstantools::posterior_linpred(fit, transform = TRUE)
  epred_counts = rstantools::posterior_epred(fit)

  expect_equal(unname(linpred_link), matrix(c(0, 0, 0), nrow = 1))
  expect_equal(unname(linpred_prob), matrix(c(0.5, 0.5, 0.5), nrow = 1))
  expect_equal(unname(epred_counts), matrix(c(2.5, 2.5, 2.5), nrow = 1))
})

test_that("loo supports by_row and soft-deprecates pointwise", {
  data = data.frame(x = 1:5, y = c(2, 4, 6, 8, 10))
  fit = mcp(list(y ~ 1 + x), data = data, sample = FALSE)
  
  population = mcp_pars(fit, scope = "population")$name
  draws = matrix(rep(c(1, 2, 0.5), each = 5), nrow = 5)
  colnames(draws) = population
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  loo_by_row = suppressWarnings(loo(fit, by_row = TRUE))
  expect_s3_class(loo_by_row, "psis_loo")
  expect_warning(
    withCallingHandlers(
      loo(fit, pointwise = TRUE),
      warning = function(w) {
        if (!grepl("deprecated", conditionMessage(w)))
          invokeRestart("muffleWarning")
      }
    ),
    "deprecated"
  )
})


test_that("find_mixture_quantile handles large count means without integer overflow", {
  cdf_fn = function(q, dpars, data, rate = FALSE) stats::ppois(q, dpars$mu)
  dpars = list(mu = c(1e8, 1.5e9))
  res = mcp:::find_mixture_quantile(cdf_fn, dpars, data = NULL, p = 0.975, is_discrete = TRUE)
  expect_true(is.numeric(res) && is.finite(res) && res > 1.5e9)
})


test_that("format.mcpfamily and summary print all family links", {
  fam = mcpfamily(stats::gaussian())
  expect_output({ res = print(fam) }, "Family: gaussian\nLinks: mu = identity; sigma = identity", fixed = TRUE)
  expect_identical(res, fam)
  expect_error(print(fam, extra = 1), class = "rlang_error")

  data = data.frame(x = 1:10, y = 1:10)
  fit_gauss = mcp(list(y ~ 1 + x), data, sample = FALSE)
  expect_output(summary(fit_gauss), "Family: gaussian\nLinks: mu = identity; sigma = identity", fixed = TRUE)

  fit_sigma = mcp(list(y ~ 1 + x + sigma(1 + x)), data, sample = FALSE)
  expect_output(summary(fit_sigma), "Family: gaussian\nLinks: mu = identity; sigma = log", fixed = TRUE)
})

