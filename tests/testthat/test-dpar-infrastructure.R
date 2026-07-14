test_that("families declare dpars independently of prior rows", {
  gaussian_family = mcpfamily(gaussian())
  poisson_family = mcpfamily(poisson())
  mean_only_families = list(
    binomial = mcpfamily(binomial()),
    bernoulli = bernoulli()
  )

  expect_equal(gaussian_family$dpar_specs$dpar, c("mu", "sigma"))
  expect_equal(gaussian_family$dpar_specs$link, c("identity", "identity"))
  expect_equal(gaussian_family$dpar_specs$implicit, c(FALSE, TRUE))
  expect_equal(gaussian_family$dpars, c("mu", "sigma", "ar", "ma"))
  expect_named(
    gaussian_family$default_prior,
    c("dpar", "par_type", "prior", "description", "condition")
  )
  expect_equal(
    gaussian_family$default_prior$prior[
      gaussian_family$default_prior$dpar == "mu" &
        gaussian_family$default_prior$par_type == "slope"
    ],
    "dt(0, max(2.5, round(mad(.y), 1)) / segment_width(.x), 3)"
  )

  expect_equal(poisson_family$dpar_specs$dpar, "mu")
  expect_equal(poisson_family$dpar_specs$link, "log")
  expect_equal(poisson_family$dpars, c("mu", "ar", "ma"))

  # Existing custom-family tables need only the original three columns.
  custom_family = gaussian()
  custom_family$default_prior = data.frame(
    dpar = "mu",
    par_type = "Intercept",
    prior = "dnorm(0, 1)"
  )
  custom_family = mcpfamily(custom_family)
  expect_equal(custom_family$default_prior$condition, "always")
  expect_true(is.na(custom_family$default_prior$description))

  for (family in mean_only_families) {
    expect_equal(family$dpar_specs$dpar, "mu")
    expect_equal(family$dpars, c("mu", "ar", "ma"))
  }

  # Removing prior metadata must not redefine the mathematical family.
  gaussian_family$default_prior = dplyr::filter(
    gaussian_family$default_prior,
    .data$dpar != "sigma"
  )
  expect_equal(gaussian_family$dpar_specs$dpar, c("mu", "sigma"))
})


test_that("family prior builders determine their supported links", {
  unsupported = gaussian()
  unsupported$link = "probit"

  expect_null(default_priors_gaussian("probit"))
  expect_null(default_priors_binomial("log"))
  expect_null(default_priors_poisson("probit"))
  expect_null(default_priors_negbinomial("identity"))
  expect_error(mcpfamily(unsupported), "no default priors for gaussian")
})


test_that("generated code separates link and distribution scales", {
  data = data.frame(x = 1:4, y = c(2, 3, 5, 7))
  fit = mcp(
    list(y ~ 1 + x + sigma(1)),
    data,
    family = gaussian(link = "log"),
    sample = FALSE
  )

  expect_match(fit$.internal$formula_jags, "link_mu_\\[i_\\] =")
  expect_match(fit$.internal$formula_jags, "link_sigma_\\[i_\\] =")
  expect_match(fit$jags_code, "mu_\\[i_\\] = exp\\(link_mu_\\[i_\\]\\)")
  expect_match(fit$jags_code, "sigma_\\[i_\\] = max\\(1e-09, link_sigma_\\[i_\\]\\)")

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

  rhs_table = get_rhs_table(
    list(y ~ 1 + x + aux(1 + x)),
    data,
    family,
    par_x = "x"
  )

  expect_setequal(unique(rhs_table$dpar), c("mu", "sigma", "aux"))
  expect_equal(
    rhs_table$code_name[rhs_table$dpar == "aux"],
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

  samples = fitted(fit, newdata = newdata, summary = FALSE, probs = FALSE) %>% add_group()
  quantile_layer = geom_quantiles(
    samples,
    quantiles = c(0.25, 0.75),
    xvar = rlang::sym("x"),
    yvar = rlang::sym("fitted"),
    facet_by = NULL
  )
  expect_equal(nrow(quantile_layer$data), 4)
  expect_equal(sort(unname(quantile_layer$data$y)), c(0, 0, 10, 10))
  expect_equal(length(unique(quantile_layer$data$.group)), 2)

  plot = ggplot2::ggplot(
    samples,
    ggplot2::aes(x = .data$x, color = .data$.group)
  ) + quantile_layer
  built_quantiles = ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(length(unique(built_quantiles$group)), 4)
})
