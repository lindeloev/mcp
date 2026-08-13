seed_model = list(y ~ 1)
seed_data = data.frame(
  x = 1:12,
  y = c(-1, 0, 1, 2, 0, 1, -1, 2, 1, 0, 2, -1)
)


test_that("get_jags_inits creates separate reproducible streams", {
  post = get_jags_inits(list(Intercept_1 = 0), 42, 3, "post")
  post_again = get_jags_inits(list(Intercept_1 = 0), 42, 3, "post")
  prior = get_jags_inits(list(Intercept_1 = 0), 42, 3, "prior")

  expect_identical(post, post_again)
  expect_length(post, 3)
  expect_length(unique(vapply(post, `[[`, integer(1), ".RNG.seed")), 3)
  expect_false(any(
    vapply(post, `[[`, integer(1), ".RNG.seed") %in%
      vapply(prior, `[[`, integer(1), ".RNG.seed")
  ))
  expect_true(all(vapply(post, `[[`, numeric(1), "Intercept_1") == 0))
})


test_that("get_jags_inits uses shared initial values and replaces RNG values", {
  inits = list(
    Intercept_1 = -1,
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = 999
  )
  seeded = get_jags_inits(inits, 42, 2, "post")

  expect_true(all(vapply(seeded, `[[`, numeric(1), "Intercept_1") == -1))
  expect_true(all(
    vapply(seeded, `[[`, character(1), ".RNG.name") ==
      "base::Wichmann-Hill"
  ))
  expect_equal(vapply(seeded, `[[`, integer(1), ".RNG.seed"), c(42L, 43L))
  expect_error(
    get_jags_inits(
      list(list(Intercept_1 = -1), list(Intercept_1 = 1)),
      42, 2, "post"
    ),
    "`inits` must be a single named list",
    fixed = TRUE
  )
})


test_that("mcp seed is validated", {
  expect_error(
    mcp(seed_model, seed_data, par_x = "x", sample = FALSE, seed = 0),
    "Element 1 is not >= 1",
    fixed = TRUE
  )
  expect_error(
    mcp(seed_model, seed_data, par_x = "x", sample = FALSE, seed = 1.5),
    "single integerish value",
    fixed = TRUE
  )
})


test_that("mcp sampling iterations are validated", {
  expect_error(
    mcp(seed_model, seed_data, par_x = "x", sample = FALSE, iter = 0),
    "Element 1 is not >= 1",
    fixed = TRUE
  )
  expect_error(
    mcp(seed_model, seed_data, par_x = "x", sample = FALSE, adapt = 1.5),
    "single integerish value",
    fixed = TRUE
  )
})


test_that("mcp seed reproduces prior and posterior samples", {
  old_plan = future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  fit_seeded = function(seed) {
    suppressWarnings(mcp(
      seed_model,
      seed_data,
      par_x = "x",
      sample = "both",
      chains = 2,
      adapt = 20,
      iter = 30,
      inits = list(Intercept_1 = 0),
      seed = seed,
      diagnostics = FALSE,
      quiet = TRUE
    ))
  }

  fit_1 = fit_seeded(123)
  fit_2 = fit_seeded(123)
  fit_3 = fit_seeded(124)

  expect_identical(.subset2(fit_1, "mcmc_post"), .subset2(fit_2, "mcmc_post"))
  expect_identical(.subset2(fit_1, "mcmc_prior"), .subset2(fit_2, "mcmc_prior"))
  expect_false(identical(.subset2(fit_1, "mcmc_post"), .subset2(fit_3, "mcmc_post")))
  expect_false(identical(.subset2(fit_1, "mcmc_prior"), .subset2(fit_3, "mcmc_prior")))
  mcmc1 = .subset2(fit_1, "mcmc_post")
  expect_false(identical(mcmc1[[1]], mcmc1[[2]]))
})
