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

  expect_equal(poisson_family$dpar_specs$dpar, "mu")
  expect_equal(poisson_family$dpar_specs$link, "log")
  expect_equal(poisson_family$dpars, c("mu", "ar", "ma"))

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
