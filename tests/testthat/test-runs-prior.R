###############
# TEST PRIORS #
###############
prior_model = list(
  y ~ 1 + x,
  1 + (1|id) ~ 1 + x,
  ~ 0
)

bad_prior = list(
  list(
    cp_1 = "dirichlet(1)",  # Has to be all-dirichlet
    cp_2 = "dnorm(3, 10)"
  ),
  list(
    cp_1 = "dirichlet(1)",
    cp_2 = "dirichlet(0)"  # alpha has to be > 0
  )
)

for (prior in bad_prior) {
  test_name = paste0("Bad priors: ", paste0(prior, collapse=", "))
  testthat::test_that(test_name, {
    testthat::expect_error(test_runs(prior_model, sample = FALSE, prior = prior))
  })
}


testthat::test_that("Prior entries must all have nonempty names", {
  testthat::expect_error(
    test_runs(prior_model, sample = FALSE, prior = list("dnorm(999, 1)")),
    "completely named list"
  )
  testthat::expect_error(
    test_runs(prior_model, sample = FALSE, prior = structure(list("dnorm(999, 1)"), names = NA_character_)),
    "completely named list"
  )
  testthat::expect_error(
    test_runs(prior_model, sample = FALSE, prior = structure(list("dnorm(999, 1)"), names = "")),
    "completely named list"
  )
})


testthat::test_that("Prior entries must have unique names", {
  prior = structure(list("dnorm(999, 1)", "dnorm(999, 1)"), names = c("cp_1", "cp_1"))
  testthat::expect_error(
    test_runs(prior_model, sample = FALSE, prior = prior),
    "duplicated entries"
  )
})


testthat::test_that("Prior entries must be scalar numbers or strings", {
  testthat::expect_error(
    validate_prior_v1(list(cp_1 = c(1, 2))),
    "finite numeric scalar or one nonempty character string"
  )
  testthat::expect_error(
    validate_prior_v1(list(cp_1 = NA_real_)),
    "finite numeric scalar or one nonempty character string"
  )
  testthat::expect_error(
    validate_prior_v1(list(cp_1 = " ")),
    "finite numeric scalar or one nonempty character string"
  )
  testthat::expect_invisible(validate_prior_v1(list(cp_1 = 1, Intercept_1 = "dnorm(0, 1)")))
})


good_prior = list(
  list(  # Fixed values and non-default change point
    Intercept_2 = "Intercept_1",
    cp_1 = "dnorm(3, 10)",
    x_2 = "-0.5"
  ),
  list(  # Outside the observed range is allowed
    cp_1 = "dunif(-100, -90)",
    cp_2 = "dnorm(100, 20) T(100, 110)"
  ),
  list(
    cp_1 = "dirichlet(1)",  # Dirichlet prior on change points
    cp_2 = "dirichlet(1)"
  ),
  list(
    cp_1 = "dirichlet(10)",  # Dirichlet prior on change points
    cp_2 = "dirichlet(10)"
  )
)

for (prior in good_prior) {
  test_name = paste0("Good priors: ", paste0(prior, collapse=", "))
  testthat::test_that(test_name, {
    test_runs(prior_model, prior = prior)
  })
}


testthat::test_that("Dirichlet change point priors use a common alpha", {
  fit = mcp(
    prior_model,
    data = data_gauss,
    prior = list(cp_1 = "dirichlet(0.5)", cp_2 = "dirichlet(0.5)"),
    par_x = "x",
    sample = FALSE,
    quiet = TRUE
  )
  testthat::expect_match(fit$jags_code, "ddirch(c(0.5, 0.5, 0.5))", fixed = TRUE)
  fit_10 = mcp(
    prior_model,
    data = data_gauss,
    prior = list(cp_1 = "dirichlet(10)", cp_2 = "dirichlet(10)"),
    par_x = "x",
    sample = FALSE,
    quiet = TRUE
  )
  testthat::expect_match(fit_10$jags_code, "ddirch(c(10, 10, 10))", fixed = TRUE)
  testthat::expect_no_error(suppressWarnings(mcp(
    prior_model,
    data = data_gauss,
    prior = list(cp_1 = "dirichlet(0.5)", cp_2 = "dirichlet(0.5)"),
    par_x = "x",
    sample = "prior",
    chains = 1,
    iter = 4,
    warmup = 4,
    quiet = TRUE
  )))
  testthat::expect_error(
    mcp(prior_model, data = data_gauss, prior = list(cp_1 = "dirichlet(2)", cp_2 = "dirichlet(3)"), par_x = "x", sample = FALSE, quiet = TRUE),
    "same alpha"
  )
  testthat::expect_error(
    mcp(prior_model, data = data_gauss, prior = list(cp_1 = "dirichlet(0)", cp_2 = "dirichlet(0)"), par_x = "x", sample = FALSE, quiet = TRUE),
    "finite alpha > 0"
  )
})


testthat::test_that("Dirichlet change point prior matches direct R simulation", {
  alpha = 2.5
  n_draws = 2000L
  model = list(y ~ 1, ~ 1, ~ 1, ~ 1, ~ 1)
  prior = stats::setNames(
    as.list(rep(paste0("dirichlet(", alpha, ")"), 4)),
    paste0("cp_", 1:4)
  )
  fit = suppressWarnings(mcp(
    model,
    data = data.frame(x = seq(0, 1, length.out = 20), y = 0),
    prior = prior,
    par_x = "x",
    sample = "prior",
    chains = 1,
    iter = n_draws,
    warmup = 100,
    quiet = TRUE
  ))
  jags_draws = as.matrix(.subset2(fit, "mcmc_prior")[[1]])[, names(prior), drop = FALSE]

  set.seed(123)
  r_spacings = matrix(stats::rgamma(n_draws * 5, shape = alpha), ncol = 5)
  r_spacings = r_spacings / rowSums(r_spacings)
  r_draws = t(apply(r_spacings, 1, cumsum))[, 1:4, drop = FALSE]

  testthat::expect_lt(max(abs(colMeans(jags_draws) - colMeans(r_draws))), 0.025)
})
