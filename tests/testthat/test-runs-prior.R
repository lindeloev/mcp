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
