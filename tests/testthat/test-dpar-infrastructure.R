test_that("families declare dpars independently of prior rows", {
  gaussian_family = mcpfamily(gaussian())
  poisson_family = mcpfamily(poisson())
  mean_only_families = list(
    binomial = mcpfamily(binomial()),
    bernoulli = bernoulli()
  )

  expect_equal(gaussian_family$dpar_specs$dpar, c("mu", "sigma"))
  expect_equal(gaussian_family$dpar_specs$link, c("identity", "identity"))
  expect_equal(gaussian_family$dpar_specs$link_constant, c("identity", "identity"))
  expect_equal(gaussian_family$dpar_specs$link_modeled, c("identity", "log"))
  expect_false(any(gaussian_family$dpar_specs$modeled))
  expect_equal(gaussian_family$dpar_specs$implicit, c(FALSE, TRUE))
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
    "dt(0, max(2.5, round(mad(.y), 1)) / predictor_scale(), 3)"
  )

  expect_equal(poisson_family$dpar_specs$dpar, "mu")
  expect_equal(poisson_family$dpar_specs$link, "log")
  expect_equal(poisson_family$dpars, "mu")

  # Compact prior tables are normalized at the family boundary.
  normalized_prior = normalize_family_default_priors(data.frame(
    dpar = "mu",
    par_type = "Intercept",
    prior = "dnorm(0, 1)"
  ))
  expect_equal(normalized_prior$condition, "always")
  expect_true(is.na(normalized_prior$description))

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


test_that("default coefficient priors scale representative predictor changes", {
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

  rounded = function(x) as.numeric(format_prior_number(x))
  scaled_t = function(base, reference_change) {
    paste0(
      "dt(0, ",
      format_prior_number(base / rounded(reference_change)),
      ", 3)"
    )
  }
  mu_scale = max(2.5, round(stats::mad(data$y), 1))
  segment_width = diff(range(data$x))
  z_scale = rounded(2 * stats::sd(data$z))

  expect_equal(fit$prior$x_1, scaled_t(mu_scale, segment_width))
  expect_equal(fit$prior$xE2_1, scaled_t(mu_scale, segment_width^2))
  expect_equal(fit$prior$z_1, scaled_t(mu_scale, z_scale))
  expect_equal(
    fit$prior$logw_1,
    scaled_t(mu_scale, 2 * stats::sd(log(data$w)))
  )
  expect_equal(fit$prior$b_1, scaled_t(mu_scale, 10))
  expect_equal(fit$prior$gb_1, "dt(0, 2.5, 3)")
  expect_equal(fit$prior$xz_1, scaled_t(mu_scale, segment_width * z_scale))
  expect_equal(
    fit$prior$zgb_1,
    scaled_t(mu_scale, 2 * stats::sd(data$z * (data$g == "b")))
  )
  expect_equal(fit$prior$sigma_z_1, scaled_t(2.5, z_scale))
  expect_equal(
    fit$prior$sigma_zx_1,
    scaled_t(2.5, segment_width * z_scale)
  )
})


test_that("built-in families implement the shared internal contract", {
  families = list(
    mcpfamily(gaussian()),
    mcpfamily(binomial()),
    bernoulli(),
    mcpfamily(poisson()),
    negbinomial()
  )

  for (family in families) {
    expect_true(is.mcpfamily(family))
    expect_true(is.function(family$r$epred))
    expect_true(is.function(family$r$log_lik))
    expect_true(is.function(family$r$rng))
    expect_true(is.function(family$backends$jags$likelihood))
    expect_true(is.list(family$response$auxiliary))
    expect_true(is.function(family$response$observed))
    expect_true(is.function(family$response$probability))
  }

  expect_true(is.function(mcpfamily(gaussian())$garma$observed_r))
  expect_null(mcpfamily(gaussian(link = "log"))$garma)
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

  expect_error(mcpfamily(unsupported_gaussian), "no default priors for gaussian")
  expect_error(mcpfamily(unsupported_binomial), "no default priors for binomial")
  expect_error(mcpfamily(unsupported_poisson), "no default priors for poisson")
  expect_error(mcpfamily(unsupported_negbinomial), "require log links")
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

  expect_match(constant_fit$jags_code, "sigma_\\[i_\\] = link_sigma_\\[i_\\]")
  expect_false(grepl("max\\(1e-09, link_sigma_", constant_fit$jags_code))
  expect_equal(constant_fit$family$links["sigma"], c(sigma = "identity"))
  expect_false(constant_fit$family$dpar_specs$modeled[
    constant_fit$family$dpar_specs$dpar == "sigma"
  ])
  expect_equal(constant_fit$prior$sigma_1, "dt(0, 2.5, 3) T(0, )")

  expect_match(fit$.internal$formula_jags, "link_mu_\\[i_\\] =")
  expect_match(fit$.internal$formula_jags, "link_sigma_\\[i_\\] =")
  expect_match(fit$jags_code, "mu_\\[i_\\] = exp\\(link_mu_\\[i_\\]\\)")
  expect_match(fit$jags_code, "sigma_\\[i_\\] = exp\\(link_sigma_\\[i_\\]\\)")
  expect_equal(fit$family$links[c("mu", "sigma")], c(mu = "log", sigma = "log"))
  expect_true(fit$family$dpar_specs$modeled[fit$family$dpar_specs$dpar == "sigma"])
  expect_equal(fit$prior$sigma_1, "dt(0, 2.5, 3)")
  expect_equal(fit$prior$sigma_x_1, "dt(0, 0.8333333, 3)")
  expect_false(grepl("max\\(1e-09, link_sigma_", fit$jags_code))

  expect_equal(change_fit$family$links["sigma"], c(sigma = "log"))
  expect_equal(
    unname(unlist(change_fit$prior[c("sigma_1", "sigma_2")])),
    rep("dt(0, 2.5, 3)", 2)
  )
  predictors = get_fit_model_tables(change_fit)$predictors
  expect_false(predictors$explicit[
    predictors$code_name == "sigma_1"
  ])
  expect_true(predictors$explicit[
    predictors$code_name == "sigma_2"
  ])

  expect_false(any(grepl("^link_", fit$pars$population)))
  expect_false(any(grepl("^(mu_|sigma_)\\[", fit$pars$population)))
})


test_that("declared dpar wrappers use the generic formula path", {
  data = data.frame(x = 1:4, y = c(2, 3, 5, 7))
  family = mcpfamily(gaussian())
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
    predictors$code_name[predictors$dpar == "aux"],
    c("aux_1", "aux_x_1")
  )
})


test_that("unsupported dpar wrappers give a family-specific error", {
  data = data.frame(x = 1:4, y = c(2, 3, 5, 7))

  expect_error(
    mcp(list(y ~ 1 + x + sigma(1)), data, family = poisson(), sample = FALSE),
    paste0(
      "`sigma()` is not a distributional parameter for family = poisson(). ",
      "See available parameters with `mcpfamily(poisson())$dpars`."
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
  family = mcpfamily(gaussian())
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
    sigma_1 = 1,
    .type = "fitted",
    .dpar = "mu"
  )
  response = rlang::exec(fit$simulate, !!!args, .scale = "response")
  linear = rlang::exec(fit$simulate, !!!args, .scale = "linear")

  expect_equal(as.numeric(response), rep(2, nrow(data)))
  expect_equal(as.numeric(linear), rep(log(2), nrow(data)))

  args$Intercept_1 = 0
  args$sigma_1 = log(3)
  args$.dpar = "sigma"
  sigma_response = rlang::exec(fit$simulate, !!!args, .scale = "response")
  sigma_linear = rlang::exec(fit$simulate, !!!args, .scale = "linear")
  expect_equal(as.numeric(sigma_response), rep(3, nrow(data)))
  expect_equal(as.numeric(sigma_linear), rep(log(3), nrow(data)))
})


test_that("evaluated-draw helpers enforce draw and data-row identity", {
  samples = expand.grid(
    data_row = 1:2,
    .draw = 1:3
  )
  samples$.chain = 1L
  samples$.iteration = samples$.draw
  samples$x = 5
  samples$group = factor(c("A", "B")[samples$data_row])
  samples$fitted = 10 * samples$data_row + samples$.draw
  samples = samples[c(6, 1, 4, 2, 5, 3), ]

  expect_silent(validate_eval_draws(samples, "fitted"))

  # Equal x values do not define the evaluation unit. Quantiles stay separate
  # by data_row and receive plotting metadata only after summarisation.
  quantiles = get_quantiles(
    samples,
    c(0.25, 0.75),
    "fitted",
    keep = c("x", "group")
  )
  expect_equal(nrow(quantiles), 4)
  expect_equal(unique(quantiles$fitted[quantiles$data_row == 1]), c(11.5, 12.5))
  expect_equal(unique(quantiles$fitted[quantiles$data_row == 2]), c(21.5, 22.5))

  fitted_matrix = tidy_to_matrix(samples, type = "fitted", data_rows = c(2, 1))
  expect_equal(colnames(fitted_matrix), c("2", "1"))
  expect_equal(
    unname(fitted_matrix),
    rbind(c(21, 11), c(22, 12), c(23, 13))
  )

  expect_error(
    validate_eval_draws(rbind(samples, samples[1, ]), "fitted"),
    "one `fitted` value per `.draw` and `data_row`",
    fixed = TRUE
  )
  expect_error(
    validate_eval_draws(samples[-1, ], "fitted"),
    "same complete set of posterior draws",
    fixed = TRUE
  )

  inconsistent_grid = samples
  inconsistent_grid$x[1] = 6
  expect_error(
    get_quantiles(
      inconsistent_grid,
      0.5,
      "fitted",
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
    sigma_1 = rep(0, 4),
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
    probs = FALSE
  )
  sigma = fitted(
    fit,
    newdata = newdata,
    dpar = "sigma",
    summary = FALSE,
    probs = FALSE
  )
  keys = c(".chain", ".iteration", ".draw", "data_row")

  expect_equal(mu[, keys], sigma[, keys])
  expect_silent(validate_eval_draws(mu, "fitted"))
  expect_silent(validate_eval_draws(sigma, "fitted"))

  sigma_summary = fitted(
    fit,
    newdata = newdata,
    dpar = "sigma",
    probs = c(0.25, 0.75)
  )
  expect_equal(sigma_summary$fitted, c(1, 3))
  expect_equal(sigma_summary$Q25, c(1, 3))
  expect_equal(sigma_summary$Q75, c(1, 3))

  predicted = pp_eval(
    fit,
    newdata = newdata,
    summary = FALSE,
    type = "predict",
    .include_fitted = TRUE
  )
  expect_equal(predicted$fitted, mu$fitted)
  expect_null(attr(predicted$predict, "fitted"))
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

  summary = fitted(fit, newdata = newdata, probs = c(0.25, 0.75))
  expect_equal(summary$fitted, c(0, 10))
  expect_equal(summary$Q25, c(0, 10))
  expect_equal(summary$Q75, c(0, 10))

  samples = fitted(fit, newdata = newdata, summary = FALSE, probs = FALSE) %>% add_plot_groups()
  quantile_layer = geom_quantiles(
    samples,
    quantiles = c(0.25, 0.75),
    xvar = rlang::sym("x"),
    yvar = rlang::sym("fitted"),
    facet_by = NULL
  )
  expect_equal(nrow(quantile_layer$data), 4)
  expect_equal(sort(unname(quantile_layer$data$fitted)), c(0, 0, 10, 10))
  expect_equal(length(unique(quantile_layer$data$.group)), 2)

  plot = ggplot2::ggplot(
    samples,
    ggplot2::aes(x = .data$x, color = .data$.group)
  ) + quantile_layer
  built_quantiles = ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(length(unique(built_quantiles$group)), 4)
})


test_that("color_by controls color without pooling categorical curves", {
  plot_colors = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  data = expand.grid(
    x = 1:4,
    group = factor(c("A", "B")),
    condition = factor(c("control", "treatment"))
  )
  data$y = 0
  fit = mcp(
    list(y ~ 1 + group + condition),
    data,
    par_x = "x",
    sample = FALSE
  )

  draws = matrix(
    0,
    nrow = 20,
    ncol = length(fit$pars$population),
    dimnames = list(NULL, fit$pars$population)
  )
  draws[, "groupB_1"] = seq(9, 11, length.out = 20)
  draws[, "conditiontreatment_1"] = seq(99, 101, length.out = 20)
  draws[, "sigma_1"] = 1
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  # NULL disables color, but the two quantiles for all four categorical
  # curves must still remain separate.
  no_color = plot(fit, color_by = NULL, lines = 0, q_fit = c(0.25, 0.75))
  no_color_quantiles = ggplot2::ggplot_build(no_color)$data[[2]]
  expect_equal(length(unique(no_color_quantiles$colour)), 1)
  expect_equal(length(unique(no_color_quantiles$group)), 8)
  expect_null(no_color$labels$colour)

  # Selecting one column changes only color: condition still defines separate
  # curves and is not pooled into the group-specific quantiles.
  one_color = plot(
    fit,
    color_by = "group",
    lines = 0,
    q_fit = c(0.25, 0.75),
    q_predict = c(0.25, 0.75)
  )
  one_color_layers = ggplot2::ggplot_build(one_color)$data
  one_color_quantiles = one_color_layers[[2]]
  expect_setequal(unique(one_color_quantiles$colour), plot_colors[1:2])
  expect_equal(length(unique(one_color_quantiles$group)), 8)
  expect_setequal(unique(one_color_layers[[3]]$colour), plot_colors[1:2])
  expect_equal(length(unique(one_color_layers[[3]]$group)), 8)
  expect_equal(one_color$labels$colour, "group")

  interaction_color = plot(
    fit,
    color_by = c("group", "condition"),
    lines = 0,
    q_fit = c(0.25, 0.75)
  )
  interaction_quantiles = ggplot2::ggplot_build(interaction_color)$data[[2]]
  expect_setequal(unique(interaction_quantiles$colour), plot_colors)
  expect_equal(interaction_color$labels$colour, "group:condition")

  faceted = plot(
    fit,
    facet_by = "condition",
    color_by = "group",
    lines = 0,
    q_fit = c(0.25, 0.75)
  )
  faceted_quantiles = ggplot2::ggplot_build(faceted)$data[[2]]
  expect_equal(length(unique(faceted_quantiles$PANEL)), 2)
  expect_equal(length(unique(faceted_quantiles$colour)), 2)

  expect_error(
    plot(fit, color_by = "does_not_exist"),
    "Invalid: 'does_not_exist'. Valid columns are 'group', 'condition'.",
    fixed = TRUE
  )
  expect_error(
    plot(fit, color_by = "x"),
    "Invalid: 'x'. Valid columns are 'group', 'condition'.",
    fixed = TRUE
  )

  dpar_plot = plot_dpar(
    fit,
    dpar = "sigma",
    color_by = "group",
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


test_that("integer varying groups can control plot color", {
  data = data.frame(
    x = rep(1:4, 2),
    id = rep(1:2, each = 4),
    y = 0
  )
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1),
    data,
    par_x = "x",
    sample = FALSE
  )

  sample_names = c(fit$pars$population, "cp_1_id[1]", "cp_1_id[2]")
  draws = matrix(
    0,
    nrow = 10,
    ncol = length(sample_names),
    dimnames = list(NULL, sample_names)
  )
  draws[, "cp_1"] = 2.5
  draws[, "cp_1_sd"] = 1
  draws[, "cp_1_id[1]"] = -0.5
  draws[, "cp_1_id[2]"] = 0.5
  draws[, "sigma_1"] = 1
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  interpolated = interpolate_newdata(fit, by = "id", x_values = 1:4)
  expect_equal(sort(unique(interpolated$id)), 1:2)

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
  expect_equal(length(unique(population_quantiles$group)), 2)

  plotted = plot(
    fit,
    color_by = "id",
    lines = 0,
    q_fit = c(0.25, 0.75),
    cp_dens = FALSE
  )
  quantiles = ggplot2::ggplot_build(plotted)$data[[2]]

  expect_equal(length(unique(quantiles$colour)), 2)
  expect_equal(length(unique(quantiles$group)), 4)
})


test_that("interpolate_newdata uses the ungrouped adaptive grid by default", {
  expect_equal(
    interpolate_newdata(demo_fit)[[demo_fit$pars$x]],
    get_x_values(demo_fit)
  )
})


test_that("fitted() intervals differ across varying levels and contain the fitted value", {
  data = data.frame(
    x = rep(1:4, 2),
    id = rep(1:2, each = 4),
    y = 0
  )
  fit = mcp(
    list(y ~ 1, 1 + (1 | id) ~ 1),
    data,
    par_x = "x",
    sample = FALSE
  )

  # cp_1_id[1] spans the query point x = 2.5, so the varying change point for
  # id 1 is sometimes before and sometimes after it (mixed segment
  # membership across draws). cp_1_id[2] stays below x = 2.5 for every draw,
  # so id 2 is deterministically in segment 2.
  sample_names = c(fit$pars$population, "cp_1_id[1]", "cp_1_id[2]")
  n_draws = 20
  draws = matrix(
    0,
    nrow = n_draws,
    ncol = length(sample_names),
    dimnames = list(NULL, sample_names)
  )
  draws[, "cp_1"] = 2.5
  draws[, "cp_1_sd"] = 1
  draws[, "cp_1_id[1]"] = seq(-1, 1, length.out = n_draws)
  draws[, "cp_1_id[2]"] = seq(-0.9, -0.6, length.out = n_draws)
  draws[, "Intercept_1"] = 0
  draws[, "Intercept_2"] = seq(9, 11, length.out = n_draws)
  draws[, "sigma_1"] = 1
  fit$mcmc_post = coda::mcmc.list(coda::mcmc(draws))

  newdata = data.frame(x = c(2.5, 2.5), id = c(1, 2))
  result = fitted(fit, newdata = newdata, probs = c(0.025, 0.975))

  # The interval must reflect each row's own varying level rather than being
  # identical (or pooled) across levels.
  expect_false(isTRUE(all.equal(result$Q2.5[1], result$Q2.5[2])))
  expect_false(isTRUE(all.equal(result$Q97.5[1], result$Q97.5[2])))

  # The fitted point estimate must fall within its own reported interval.
  expect_true(all(result$fitted >= result$Q2.5))
  expect_true(all(result$fitted <= result$Q97.5))
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
    c("sigma_groupB_1", "sigma_groupC_1", "sigma_groupAx_1", "sigma_groupBx_1", "sigma_groupCx_1") %in%
      fit_sigma$pars$population
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
    c("ar1_groupB_1", "ar1_groupC_1", "ar1_groupAx_1", "ar1_groupBx_1", "ar1_groupCx_1") %in%
      fit_ar$pars$population
  ))
})
