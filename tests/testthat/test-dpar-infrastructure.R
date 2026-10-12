test_that("families define dpars independently of prior rows", {
  gaussian_family = as_mcpfamily(gaussian())
  poisson_family = as_mcpfamily(poisson())
  mean_only_families = list(
    binomial = as_mcpfamily(binomial()),
    bernoulli = bernoulli()
  )

  expect_equal(gaussian_family$dpar_specs$dpar, c("mu", "sigma"))
  expect_equal(gaussian_family$dpar_specs$link, c("identity", "identity"))
  expect_equal(gaussian_family$dpar_specs$link_constant, c("identity", "identity"))
  expect_equal(gaussian_family$dpar_specs$link_modeled, c("identity", "log"))
  expect_false(any(gaussian_family$dpar_specs$modeled))
  expect_equal(gaussian_family$dpar_specs$implicit, c(FALSE, TRUE))
  expect_equal(gaussian_family$dpar_specs$require_initial_predictor, c(FALSE, TRUE))
  expect_equal(gaussian_family$dpars, c("mu", "sigma"))
  expect_named(
    gaussian_family$default_prior,
    c("dpar", "par_type", "prior", "group_sd_prior", "description", "condition")
  )
  expect_equal(
    gaussian_family$default_prior$prior[
      gaussian_family$default_prior$dpar == "mu" &
        gaussian_family$default_prior$par_type == "slope"
    ],
    "dnorm(0, 2.5 * sd(.y) / predictor_scale())"
  )

  expect_equal(poisson_family$dpar_specs$dpar, "mu")
  expect_equal(poisson_family$dpar_specs$link, "log")
  expect_equal(poisson_family$dpars, "mu")

  # Compact prior tables are standardized in as_mcpfamily().
  standardized_prior = standardize_family_default_priors(data.frame(
    dpar = "mu",
    par_type = "Intercept",
    prior = "dnorm(0, 1)"
  ))
  expect_equal(standardized_prior$condition, "always")
  expect_true(is.na(standardized_prior$description))

  for (family in mean_only_families) {
    expect_equal(family$dpar_specs$dpar, "mu")
    expect_equal(family$dpars, "mu")
  }

  # Removing prior metadata must not redefine the mathematical family.
  gaussian_family$default_prior = dplyr::filter(
    gaussian_family$default_prior,
    .data$dpar != "sigma"
  )
  expect_equal(gaussian_family$dpar_specs$dpar, c("mu", "sigma"))
})


test_that("default coefficient priors are autoscaled by predictor SDs (rstanarm)", {
  x = 1:30
  data = data.frame(
    x = x,
    z = (x * 7) %% 31,
    w = exp((((x * 11) %% 37) + 1) / 10),
    b = rep(c(0, 10), 15),
    g = factor(rep(c("a", "b", "c"), 10)),
    y = sin(x / 3) + x / 10
  )
  fit = mcp(
    list(
      y ~ 1 + x + I(x^2) + z + log(w) + b + g + x:z + z:g +
        sigma(1 + z + x:z)
    ),
    data,
    family = gaussian(),
    par_x = "x",
    sample = FALSE
  )

  # Platforms can differ in the last bit of floating-point arithmetic, so allow ~1e-15
  expect_scaled = function(prior, dist, base, predictor_sd) {
    expect_true(startsWith(prior, paste0(dist, "(0, ")))
    expect_equal(as.numeric(strsplit(prior, "[(,)]")[[1]][3]), base / predictor_sd, tolerance = 1e-15)
  }
  mu_base = 2.5 * stats::sd(data$y)
  x_sd = stats::sd(data$x)
  z_sd = stats::sd(data$z)

  expect_scaled(fit$prior$x_1, "dnorm", mu_base, x_sd)
  expect_scaled(fit$prior$xE2_1, "dnorm", mu_base, stats::sd((data$x - min(data$x))^2))
  expect_scaled(fit$prior$z_1, "dnorm", mu_base, z_sd)
  expect_scaled(fit$prior$logw_1, "dnorm", mu_base, stats::sd(log(data$w)))
  expect_scaled(fit$prior$b_1, "dnorm", mu_base, stats::sd(data$b))
  expect_scaled(fit$prior$gb_1, "dnorm", mu_base, stats::sd(data$g == "b"))
  expect_scaled(fit$prior$x.z_1, "dnorm", mu_base, x_sd * z_sd)
  expect_scaled(fit$prior$z.gb_1, "dnorm", mu_base, stats::sd(data$z * (data$g == "b")))
  expect_scaled(fit$prior$sigma_z_1, "dt", 2.5, z_sd)
  expect_scaled(fit$prior$sigma_z.x_1, "dt", 2.5, x_sd * z_sd)
})


test_that("built-in families implement the shared internal contract", {
  families = list(
    as_mcpfamily(gaussian()),
    as_mcpfamily(binomial()),
    bernoulli(),
    as_mcpfamily(poisson()),
    negbinomial()
  )

  for (family in families) {
    expect_true(is_mcpfamily(family))
    expect_identical(assert_mcpfamily(family), family)
    expect_true(is.function(family$r$epred))
    expect_true(is.function(family$r$log_lik))
    expect_true(is.function(family$r$rng))
    expect_true(is.function(family$r$quantile))
    expect_true(all(family$dpars %in% known_dpar_wrappers()))  # Informative errors for wrappers used with other families
    expect_true(is.function(family$backends$jags$likelihood))
    expect_true(is.list(family$response$auxiliary))
    expect_true(is.function(family$response$observed))
    expect_true(is.function(family$response$probability))
  }

  expect_true(is.function(as_mcpfamily(gaussian())$garma$observed_r))
  expect_null(as_mcpfamily(gaussian(link = "log"))$garma)
})


test_that("family prior builders determine their supported links", {
  unsupported_gaussian = gaussian()
  unsupported_gaussian$link = "probit"
  unsupported_binomial = binomial()
  unsupported_binomial$link = "log"
  unsupported_poisson = poisson()
  unsupported_poisson$link = "probit"
  unsupported_negbinomial = negbinomial()
  unsupported_negbinomial$link = "identity"

  expect_error(as_mcpfamily(unsupported_gaussian), "no default priors for gaussian")
  expect_error(as_mcpfamily(unsupported_binomial), "no default priors for binomial")
  expect_error(as_mcpfamily(unsupported_poisson), "no default priors for poisson")
  expect_error(as_mcpfamily(unsupported_negbinomial), "require log links")
})


test_that("generated code separates link and distribution scales", {
  data = data.frame(x = 1:4, y = c(2, 3, 5, 7))
  constant_fit = mcp(
    list(y ~ 1 + x),
    data,
    family = gaussian(link = "log"),
    sample = FALSE
  )
  fit = mcp(
    list(y ~ 1 + x + sigma(1 + x)),
    data,
    family = gaussian(link = "log"),
    sample = FALSE
  )
  change_fit = mcp(
    list(y ~ 1, ~ 0 + sigma(1)),
    data,
    par_x = "x",
    sample = FALSE
  )

  expect_match(constant_fit$jags_code, "sigma_\\[i_\\] = max\\(1e-05, link_sigma_\\[i_\\]\\)")
  expect_false(grepl("max\\(1e-09, link_sigma_", constant_fit$jags_code))
  expect_equal(constant_fit$family$links["sigma"], c(sigma = "identity"))
  expect_false(constant_fit$family$dpar_specs$modeled[
    constant_fit$family$dpar_specs$dpar == "sigma"
  ])
  expect_equal(constant_fit$prior$sigma_1, "dt(0, 2.5, 3) T(0, )")

  expect_match(fit$.internal$formula_jags, "link_mu_\\[i_\\] =")
  expect_match(fit$.internal$formula_jags, "link_sigma_\\[i_\\] =")
  expect_match(fit$jags_code, "mu_\\[i_\\] = exp\\(link_mu_\\[i_\\]\\)")
  expect_match(fit$jags_code, "sigma_\\[i_\\] = max\\(1e-05, exp\\(link_sigma_\\[i_\\]\\)\\)")
  expect_equal(fit$family$links[c("mu", "sigma")], c(mu = "log", sigma = "log"))
  expect_true(fit$family$dpar_specs$modeled[fit$family$dpar_specs$dpar == "sigma"])
  expect_equal(fit$prior$sigma_Intercept_1, "dt(-sigma_x_1 * 1, 2.5, 3)")
  expect_equal(fit$prior$sigma_x_1, paste0("dt(0, ", format_prior_number(2.5 / stats::sd(data$x)), ", 3)"))
  expect_false(grepl("max\\(1e-09, link_sigma_", fit$jags_code))

  expect_equal(change_fit$family$links["sigma"], c(sigma = "log"))
  expect_equal(
    unname(unlist(change_fit$prior[c("sigma_Intercept_1", "sigma_Intercept_2")])),
    rep("dt(0, 2.5, 3)", 2)
  )
  expect_true(change_fit$family$dpar_specs$modeled[change_fit$family$dpar_specs$dpar == "sigma"])

  expect_false(any(grepl("^link_", mcp_pars(fit, scope = "population")$parameter)))
  expect_false(any(grepl("^(mu_|sigma_)\\[", mcp_pars(fit, scope = "population")$parameter)))
})


test_that("included dpar wrappers use the generic formula path", {
  data = data.frame(x = 1:4, y = c(2, 3, 5, 7))
  family = as_mcpfamily(gaussian())
  family$dpar_specs = dplyr::bind_rows(
    family$dpar_specs,
    new_dpar_spec("aux", "log", implicit = TRUE)
  )
  family = add_dpar_specs(family)

  predictors = get_predictors(
    list(y ~ 1 + x + aux(1 + x)),
    data,
    family,
    par_x = "x"
  )

  expect_setequal(unique(predictors$dpar), c("mu", "sigma", "aux"))
  expect_equal(
    predictors$name[predictors$dpar == "aux"],
    c("aux_Intercept_1", "aux_x_1")
  )
})

test_that("distributional population terms apply only where included", {
  data = data.frame(
    x = 1:8,
    z = c(0, 1, 0, 2, 1, 3, 2, 4),
    y = c(2, 3, 5, 4, 6, 8, 7, 9)
  )
  fit = mcp(
    list(
      y ~ 1 + sigma(1 + z),
      ~ 1,
      ~ 0 + sigma(0 + z)
    ),
    data,
    par_x = "x",
    sample = FALSE
  )
  expect_lifetimes(fit, c(sigma_z_1 = 2L))
  expect_match(
    fit$.internal$formula_jags,
    "(?s)x\\[i_\\] < cp_1.*c\\(sigma_z_1\\)",
    perl = TRUE
  )

  sigma = fit$simulate(
    fit,
    data,
    cp_1 = 3.5, cp_2 = 5.5,
    Intercept_1 = 0, Intercept_2 = 0,
    sigma_Intercept_1 = 0.5, sigma_z_1 = 1, sigma_z_3 = 2,
    .type = "fitted", .dpar = "sigma", .scale = "linear"
  )
  expected = 0.5 + ifelse(data$x < 3.5, data$z, ifelse(data$x < 5.5, 0, 2 * data$z))
  expect_equal(as.numeric(sigma), expected)

  shape = get_predictors(
    list(y ~ 1 + shape(1 + z), ~ 0 + shape(0 + z)),
    data,
    as_mcpfamily(negbinomial()),
    par_x = "x"
  )
  expect_lifetimes(shape, c(shape_z_1 = 2L))
})


test_that("required dpars have an initial predictor and implicit defaults retain their provenance", {
  data = data.frame(x = 1:6, z = 1:6, y = c(1, 2, 1, 3, 2, 4))

  expect_error(
    mcp(list(y ~ 1 + sigma(0)), data, par_x = "x", sample = FALSE),
    "`sigma\\(0\\)` cannot be the initial predictor"
  )
  expect_error(
    mcp(list(y ~ 1 + shape(0)), data, family = negbinomial(), par_x = "x", sample = FALSE),
    "`shape\\(0\\)` cannot be the initial predictor"
  )
  expect_silent(mcp(list(y ~ 1 + sigma(0 + x)), data, par_x = "x", sample = FALSE))
  expect_silent(mcp(list(y ~ 1 + shape(0 + z)), data, family = negbinomial(), par_x = "x", sample = FALSE))

  fit = mcp(list(y ~ 1, ~ 1 + sigma(1), ~ 1), data, par_x = "x", sample = FALSE)
  expect_false("sigma_Intercept_3" %in% mcp_pars(fit)$parameter)
  expect_true("sigma_Intercept_2" %in% mcp_pars(fit)$parameter)

  implicit = mcp(list(y ~ 1, ~ 1), data, par_x = "x", sample = FALSE)
  expect_equal(implicit$family$links["sigma"], c(sigma = "identity"))
  expect_equal(mcp_pars(implicit)$parameter[grepl("^sigma", mcp_pars(implicit)$parameter)], "sigma_1")
})


test_that("each component is allowed once per segment", {
  data = data.frame(x = 1:4, y = c(2, 3, 5, 7))
  expect_error(mcp(list(y ~ 1 + sigma(1) + sigma(1)), data, par_x = "x", sample = FALSE), "Only one `sigma()` allowed", fixed = TRUE)
  expect_error(mcp(list(y ~ mu(1) + mu(x)), data, par_x = "x", sample = FALSE), "Only one `mu()` allowed", fixed = TRUE)
  expect_error(mcp(list(y ~ 1 + sigma(1, x)), data, par_x = "x", sample = FALSE), "`sigma()` takes exactly one formula", fixed = TRUE)
})


test_that("unsupported dpar wrappers give a family-specific error", {
  data = data.frame(x = 1:4, y = c(2, 3, 5, 7))

  expect_error(
    mcp(list(y ~ 1 + x + sigma(1)), data, family = poisson(), sample = FALSE),
    paste0(
      "`sigma()` is not a distributional parameter for family = poisson(). ",
      "See available parameters with `as_mcpfamily(poisson())$dpars`."
    ),
    fixed = TRUE
  )

  expect_error(
    mcp(list(y ~ 1 + x + shape(1)), data, family = gaussian(), sample = FALSE),
    "`shape()` is not a distributional parameter for family = gaussian()",
    fixed = TRUE
  )
})


test_that("AR and MA component names are reserved", {
  family = as_mcpfamily(gaussian())
  family$dpar_specs = dplyr::bind_rows(
    family$dpar_specs,
    new_dpar_spec("ma", "identity")
  )

  expect_error(
    add_dpar_specs(family),
    "'epred', 'ar', and 'ma' are reserved and cannot be family distributional parameters.",
    fixed = TRUE
  )
})


test_that("R-side dpar evaluation supports link and response scales", {
  data = data.frame(x = 1:4, y = c(2, 3, 5, 7))
  fit = mcp(
    list(y ~ 1 + x + sigma(1)),
    data,
    family = gaussian(link = "log"),
    sample = FALSE
  )

  args = list(
    fit,
    data,
    Intercept_1 = log(2),
    x_1 = 0,
    sigma_Intercept_1 = 1,
    .type = "fitted",
    .dpar = "mu"
  )
  response = rlang::exec(fit$simulate, !!!args, .scale = "response")
  linear = rlang::exec(fit$simulate, !!!args, .scale = "linear")

  expect_equal(as.numeric(response), rep(2, nrow(data)))
  expect_equal(as.numeric(linear), rep(log(2), nrow(data)))

  args$Intercept_1 = 0
  args$sigma_Intercept_1 = log(3)
  args$.dpar = "sigma"
  sigma_response = rlang::exec(fit$simulate, !!!args, .scale = "response")
  sigma_linear = rlang::exec(fit$simulate, !!!args, .scale = "linear")
  expect_equal(as.numeric(sigma_response), rep(3, nrow(data)))
  expect_equal(as.numeric(sigma_linear), rep(log(3), nrow(data)))
})


test_that("evaluated-draw helpers enforce draw and data-row identity", {
  samples = expand.grid(
    .data_row = 1:2,
    .draw = 1:3
  )
  samples$.chain = 1L
  samples$.iteration = samples$.draw
  samples$x = 5
  samples$group = factor(c("A", "B")[samples$.data_row])
  samples$`.epred` = 10 * samples$.data_row + samples$.draw
  samples = samples[c(6, 1, 4, 2, 5, 3), ]

  expect_silent(validate_eval_draws(samples, ".epred"))

  # Equal x values do not define the evaluation unit. Quantiles stay separate
  # by .data_row and receive plotting metadata only after summarisation.
  quantiles = get_quantiles(
    samples,
    c(0.25, 0.75),
    ".epred",
    keep = c("x", "group")
  )
  expect_equal(nrow(quantiles), 4)
  expect_equal(unique(quantiles$.epred[quantiles$.data_row == 1]), c(11.5, 12.5))
  expect_equal(unique(quantiles$.epred[quantiles$.data_row == 2]), c(21.5, 22.5))

  fitted_matrix = tidy_to_matrix(samples, type = ".epred", data_rows = c(2, 1))
  expect_equal(colnames(fitted_matrix), c("2", "1"))
  expect_equal(
    unname(fitted_matrix),
    rbind(c(21, 11), c(22, 12), c(23, 13))
  )

  expect_error(
    validate_eval_draws(rbind(samples, samples[1, ]), ".epred"),
    "one `.epred` value per `.draw` and `.data_row`",
    fixed = TRUE
  )
  expect_error(
    validate_eval_draws(samples[-1, ], ".epred"),
    "same complete set of posterior draws",
    fixed = TRUE
  )

  inconsistent_grid = samples
  inconsistent_grid$x[1] = 6
  expect_error(
    get_quantiles(
      inconsistent_grid,
      0.5,
      ".epred",
      keep = "x"
    ),
    "metadata differs across draws",
    fixed = TRUE
  )
})


test_that("distributional parameters use the same evaluation identity", {
  data = data.frame(
    x = rep(1:3, each = 2),
    group = factor(rep(c("A", "B"), 3)),
    y = 0
  )
  fit = mcp(
    list(y ~ 1 + group + sigma(1 + group)),
    data,
    par_x = "x",
    sample = FALSE
  )

  draws = cbind(
    Intercept_1 = rep(0, 4),
    groupB_1 = rep(10, 4),
    sigma_Intercept_1 = rep(0, 4),
    sigma_groupB_1 = rep(log(3), 4)
  )
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))
  newdata = data.frame(
    x = c(1, 1),
    group = factor(c("A", "B"), levels = levels(data$group))
  )

  mu = fitted(
    fit,
    newdata = newdata,
    dpar = "mu",
    summary = FALSE,
    width = FALSE
  )
  sigma = fitted(
    fit,
    newdata = newdata,
    dpar = "sigma",
    summary = FALSE,
    width = FALSE
  )
  keys = c(".chain", ".iteration", ".draw", ".data_row")

  expect_equal(mu[, keys], sigma[, keys])
  expect_silent(validate_eval_draws(mu, ".mu"))
  expect_silent(validate_eval_draws(sigma, ".sigma"))

  sigma_summary = fitted(
    fit,
    newdata = newdata,
    dpar = "sigma",
    width = 0.5
  )
  expect_equal(sigma_summary$.mean, c(1, 3))
  expect_equal(sigma_summary$.lower, c(1, 3))
  expect_equal(sigma_summary$.upper, c(1, 3))

  predicted = pp_eval(
    fit,
    newdata = newdata,
    summary = FALSE,
    type = "predict",
    .return_fitted = TRUE
  )
  expect_equal(predicted$.epred, mu$.mu)
  expect_null(attr(predicted$.prediction, ".epred"))
})


test_that("per-draw value columns are named after the dpar and scale", {
  value_col = function(...) utils::tail(names(fitted(demo_fit, summary = FALSE, ndraws = 2, ...)), 1)
  expect_equal(value_col(), ".epred")
  expect_equal(value_col(dpar = "mu"), ".mu")
  expect_equal(value_col(dpar = "sigma"), ".sigma")
  expect_equal(value_col(scale = "linear"), ".linpred")
  expect_equal(value_col(dpar = "sigma", scale = "linear"), ".linpred")
  expect_error(
    fitted(demo_fit, newdata = data.frame(time = 1:2, .sigma = 1:2), dpar = "sigma"),
    "reserved output namespace"
  )
})


test_that("quantiles stay separate for rows and categorical curves sharing x", {
  data = data.frame(
    x = 1:4,
    group = factor(c("A", "B", "A", "B")),
    y = c(0, 10, 0, 10)
  )
  fit = mcp(list(y ~ 1 + group), data, par_x = "x", sample = FALSE)
  draws = cbind(
    Intercept_1 = c(0, 0),
    groupB_1 = c(10, 10),
    sigma_1 = c(1, 1)
  )
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))
  newdata = data.frame(x = c(2, 2), group = c("A", "B"))

  summary = fitted(fit, newdata = newdata, width = 0.5)
  expect_equal(summary$.mean, c(0, 10))
  expect_equal(summary$.lower, c(0, 10))
  expect_equal(summary$.upper, c(0, 10))

  samples = fitted(fit, newdata = newdata, summary = FALSE, width = FALSE) %>% add_plot_groups()
  quantiles = get_quantiles(samples, c(0.25, 0.75), ".epred", c("x", ".group", ".color"))
  quantile_layer = geom_quantiles(
    quantiles,
    yvar = rlang::sym(".epred"),
    use_color = TRUE
  )
  expect_equal(nrow(quantile_layer$data), 4)
  expect_equal(sort(unname(quantile_layer$data$.epred)), c(0, 0, 10, 10))
  expect_equal(length(unique(quantile_layer$data$.group)), 2)

  automatic = plot(fit, lines = 0, q_fit = c(0.25, 0.75))
  automatic_quantiles = ggplot2::ggplot_build(automatic)$data[[2]]
  expect_equal(automatic$labels$colour, "group")
  expect_equal(length(unique(automatic_quantiles$colour)), 2)

  plot = ggplot2::ggplot(
    samples,
    ggplot2::aes(x = .data$x, color = .data$.group)
  ) + quantile_layer
  built_quantiles = ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(length(unique(built_quantiles$group)), 4)
})


test_that("color_by controls color without pooling categorical curves", {
  expect_equal(formals(plot.mcpfit)$ndraws, 500)
  expect_equal(formals(plot_dpar)$ndraws, 500)

  plot_colors = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  data = expand.grid(
    x = 1:4,
    group = factor(c("A", "B")),
    state = factor(c("control", "treatment"))
  )
  data$y = 0
  fit = mcp(
    list(y ~ 1 + group + state),
    data,
    par_x = "x",
    sample = FALSE
  )

  population = mcp_pars(fit, scope = "population")$parameter
  draws = matrix(
    0,
    nrow = 20,
    ncol = length(population),
    dimnames = list(NULL, population)
  )
  draws[, "groupB_1"] = seq(9, 11, length.out = 20)
  draws[, "statetreatment_1"] = seq(99, 101, length.out = 20)
  draws[, "sigma_1"] = 1
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  expect_error(
    plot(fit, lines = 0, q_fit = c(0.25, 0.75)),
    "Unmapped categorical predictors: 'group', 'state'.",
    fixed = TRUE
  )
  expect_error(
    plot(fit, color_by = "group", lines = 0, q_fit = c(0.25, 0.75)),
    "Unmapped categorical predictors: 'state'.",
    fixed = TRUE
  )

  # NULL disables color, but the two quantiles for all four categorical
  # curves must still remain separate.
  no_color = plot(fit, color_by = NULL, lines = 0, q_fit = c(0.25, 0.75), q_predict = c(0.25, 0.75))
  no_color_layers = ggplot2::ggplot_build(no_color)$data
  no_color_quantiles = no_color_layers[[2]]
  expect_equal(unique(no_color_quantiles$colour), "red")
  expect_equal(unique(no_color_quantiles$linetype), "dashed")
  expect_equal(unique(no_color_quantiles$linewidth), 0.85)
  expect_equal(unique(no_color_layers[[3]]$colour), "#009E73")
  expect_equal(unique(no_color_layers[[3]]$linetype), "dashed")
  expect_equal(length(unique(no_color_quantiles$group)), 8)
  expect_null(no_color$labels$colour)

  # All categorical curves are identified by either color or facets.
  one_color = plot(
    fit,
    color_by = "group",
    facet_by = "state",
    lines = 0,
    q_fit = c(0.25, 0.75),
    q_predict = c(0.25, 0.75)
  )
  one_color_layers = ggplot2::ggplot_build(one_color)$data
  one_color_quantiles = one_color_layers[[2]]
  expect_setequal(unique(one_color_quantiles$colour), plot_colors[1:2])
  expect_equal(unique(one_color_quantiles$linetype), "dashed")
  expect_equal(length(unique(one_color_quantiles$group)), 8)
  expect_setequal(unique(one_color_layers[[3]]$colour), plot_colors[1:2])
  expect_equal(unique(one_color_layers[[3]]$linetype), "dashed")
  expect_equal(length(unique(one_color_layers[[3]]$group)), 8)
  expect_equal(one_color$labels$colour, "group")

  interaction_color = plot(
    fit,
    color_by = c("group", "state"),
    lines = 0,
    q_fit = c(0.25, 0.75)
  )
  interaction_quantiles = ggplot2::ggplot_build(interaction_color)$data[[2]]
  expect_setequal(unique(interaction_quantiles$colour), plot_colors)
  expect_equal(interaction_color$labels$colour, "group:state")

  faceted = plot(
    fit,
    facet_by = "state",
    color_by = "group",
    lines = 0,
    q_fit = c(0.25, 0.75)
  )
  faceted_quantiles = ggplot2::ggplot_build(faceted)$data[[2]]
  expect_equal(length(unique(faceted_quantiles$PANEL)), 2)
  expect_equal(length(unique(faceted_quantiles$colour)), 2)

  expect_error(
    plot(fit, color_by = "does_not_exist"),
    "Invalid: 'does_not_exist'. Valid columns are 'group', 'state'.",
    fixed = TRUE
  )
  expect_error(
    plot(fit, color_by = "x"),
    "Invalid: 'x'. Valid columns are 'group', 'state'.",
    fixed = TRUE
  )

  dpar_plot = plot_dpar(
    fit,
    dpar = "sigma",
    color_by = "group",
    facet_by = "state",
    lines = 0,
    q_fit = c(0.25, 0.75)
  )
  dpar_quantiles = ggplot2::ggplot_build(dpar_plot)$data[[1]]
  expect_setequal(unique(dpar_quantiles$colour), plot_colors[1:2])
  expect_equal(length(unique(dpar_quantiles$group)), 8)

  # A posterior draw is joint across curves, so every curve should use the
  # same sampled draw IDs rather than sampling independently by group.
  plotted = plot(fit, color_by = "group", lines = 3)
  line_data = plotted$layers[[2]]$data
  draws_by_curve = split(line_data$.draw, line_data$.group)
  draws_by_curve = lapply(draws_by_curve, unique)

  expect_length(unique(line_data$.draw), 3)
  expect_true(all(vapply(
    draws_by_curve,
    identical,
    logical(1),
    draws_by_curve[[1]]
  )))
})


test_that("change point density colors adapt to fitted colors", {
  plots = list(plot(demo_fit, lines = 1), plot_dpar(demo_fit, dpar = "sigma", lines = 1))
  for (plot in plots) {
    layer = which(vapply(plot$layers, function(x) inherits(x$geom, "GeomPolygon"), logical(1)))
    density = ggplot2::ggplot_build(plot)$data[[layer]]
    expect_equal(unique(density$fill), "#6BAED6")
    expect_equal(unique(density$colour), "#3182BD")
    expect_equal(unique(density$linewidth), 0.25)
    expect_equal(unique(density$alpha), 0.4 / nchains(demo_fit))
  }

  gray_density = geom_cp_density(demo_fit, NULL, FALSE, c(0, 1), use_color = TRUE)
  expect_equal(gray_density$aes_params$fill, "#A6A6A6")
  expect_equal(gray_density$aes_params$colour, "#666666")
})


test_that("integer grouping factors can control plot color", {
  data = expand.grid(x = 1:4, id = 1:2, state = factor(c("A", "B")))
  data$y = 0
  fit = mcp(
    list(y ~ 1 + state, ~ cp(1 + (1 | id)) + 1),
    data,
    par_x = "x",
    sample = FALSE
  )

  sample_names = c(mcp_pars(fit, scope = "population")$parameter, "cp_1_id[1]", "cp_1_id[2]")
  draws = matrix(
    0,
    nrow = 10,
    ncol = length(sample_names),
    dimnames = list(NULL, sample_names)
  )
  draws[, "cp_1"] = 2.5
  draws[, "cp_1_id_sd"] = 1
  draws[, "cp_1_id[1]"] = -0.5
  draws[, "cp_1_id[2]"] = 0.5
  draws[, "sigma_1"] = 1
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  interpolated = interpolate_newdata(fit, by = "id", x_values = 1:4)
  expect_equal(sort(unique(interpolated$id)), 1:2)
  expect_false("id" %in% names(interpolate_newdata(fit, x_values = 1:4)))

  automatic = plot(fit, lines = 0, q_fit = c(0.25, 0.75), cp_dens = FALSE)
  automatic_quantiles = ggplot2::ggplot_build(automatic)$data[[2]]
  expect_equal(automatic$labels$colour, "state")
  expect_equal(length(unique(automatic_quantiles$colour)), 2)

  # Without an explicit facet or color mapping, show the population curve
  # once rather than duplicating it for every integer group ID.
  population_plot = plot(
    fit,
    color_by = NULL,
    lines = 0,
    q_fit = c(0.25, 0.75),
    cp_dens = FALSE
  )
  population_quantiles = ggplot2::ggplot_build(population_plot)$data[[2]]
  expect_equal(length(unique(population_quantiles$group)), 4)

  plotted = plot(
    fit,
    color_by = "id",
    lines = 0,
    q_fit = c(0.25, 0.75),
    cp_dens = FALSE
  )
  quantiles = ggplot2::ggplot_build(plotted)$data[[2]]

  expect_equal(length(unique(quantiles$colour)), 2)
  expect_equal(length(unique(quantiles$group)), 8)

  grouped = plot(
    fit,
    color_by = "state",
    facet_by = "id",
    lines = 0,
    q_fit = c(0.25, 0.75),
    cp_dens = FALSE
  )
  grouped_quantiles = ggplot2::ggplot_build(grouped)$data[[2]]
  expect_equal(length(unique(grouped_quantiles$PANEL)), 2)
  expect_equal(length(unique(grouped_quantiles$colour)), 2)
  expect_equal(length(unique(grouped_quantiles$group)), 8)

  # color_by = "auto" still colors the sole categorical predictor when faceting by another column
  auto_grouped = plot(fit, facet_by = "id", lines = 0, q_fit = c(0.25, 0.75), cp_dens = FALSE)
  expect_equal(auto_grouped$labels$colour, "state")
  expect_equal(ggplot2::ggplot_build(auto_grouped)$data[[2]], grouped_quantiles)
  expect_equal(plot_dpar(fit, dpar = "sigma", facet_by = "id", lines = 0, cp_dens = FALSE)$labels$colour, "state")

  # ... but not when it is already faceted
  auto_faceted = plot(fit, facet_by = "state", lines = 0, cp_dens = FALSE)
  expect_null(auto_faceted$labels$colour)
})


test_that("interpolate_newdata uses the ungrouped adaptive grid by default", {
  expect_equal(
    interpolate_newdata(demo_fit)[[mcp_columns(demo_fit)$par_x]],
    get_x_values(demo_fit)
  )

  # The prior grid is refined around the prior change points
  prior_x = interpolate_newdata(demo_fit, prior = TRUE)[[mcp_columns(demo_fit)$par_x]]
  expect_equal(prior_x, get_x_values(demo_fit, prior = TRUE))
  prior_cps = as_draws_df(demo_fit, prior = TRUE)$cp_1
  expect_true(all(stats::quantile(prior_cps, c(0, 1)) %in% prior_x))
})


test_that("fitted() intervals differ across group levels and contain the fitted value", {
  data = data.frame(
    x = rep(1:4, 2),
    id = rep(1:2, each = 4),
    y = 0
  )
  fit = mcp(
    list(y ~ 1, ~ cp(1 + (1 | id)) + 1),
    data,
    par_x = "x",
    sample = FALSE
  )

  # cp_1_id[1] spans the query point x = 2.5, so the group-level change point for
  # id 1 is sometimes before and sometimes after it (mixed segment
  # membership across draws). cp_1_id[2] stays below x = 2.5 for every draw,
  # so id 2 is deterministically in segment 2.
  sample_names = c(mcp_pars(fit, scope = "population")$parameter, "cp_1_id[1]", "cp_1_id[2]")
  n_draws = 20
  draws = matrix(
    0,
    nrow = n_draws,
    ncol = length(sample_names),
    dimnames = list(NULL, sample_names)
  )
  draws[, "cp_1"] = 2.5
  draws[, "cp_1_id_sd"] = 1
  draws[, "cp_1_id[1]"] = seq(-1, 1, length.out = n_draws)
  draws[, "cp_1_id[2]"] = seq(-0.9, -0.6, length.out = n_draws)
  draws[, "Intercept_1"] = 0
  draws[, "Intercept_2"] = seq(9, 11, length.out = n_draws)
  draws[, "sigma_1"] = 1
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  newdata = data.frame(x = c(2.5, 2.5), id = c(1, 2))
  result = fitted(fit, newdata = newdata, width = 0.95)

  # The interval must reflect each row's own group level rather than being
  # identical (or pooled) across levels.
  expect_false(isTRUE(all.equal(result$.lower[1], result$.lower[2])))
  expect_false(isTRUE(all.equal(result$.upper[1], result$.upper[2])))

  # The fitted point estimate must fall within its own reported interval.
  expect_true(all(result$.mean >= result$.lower))
  expect_true(all(result$.mean <= result$.upper))
})


test_that("sigma() and ar() support categorical predictors and interactions", {
  data = data.frame(
    x = rep(1:4, 3),
    group = rep(c("A", "B", "C"), each = 4),
    y = 0
  )

  fit_sigma = mcp(
    list(y ~ 1 + sigma(1 + group + x:group)),
    data,
    par_x = "x",
    sample = FALSE
  )
  expect_true(all(
    c("sigma_groupB_1", "sigma_groupC_1", "sigma_groupA.x_1", "sigma_groupB.x_1", "sigma_groupC.x_1") %in%
      mcp_pars(fit_sigma, scope = "population")$parameter
  ))

  # x is not sorted across groups here, which is irrelevant for parsing but
  # triggers an unrelated ar()/ma() row-order notice.
  fit_ar = suppressMessages(mcp(
    list(y ~ ar(1, 1 + group + x:group)),
    data,
    par_x = "x",
    sample = FALSE
  ))
  expect_true(all(
    c("ar1_groupB_1", "ar1_groupC_1", "ar1_groupA.x_1", "ar1_groupB.x_1", "ar1_groupC.x_1") %in%
      mcp_pars(fit_ar, scope = "population")$parameter
  ))
})
