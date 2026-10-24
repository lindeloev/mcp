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

  # Nearby seeds share no chains
  chain_seeds = lapply(1:20, function(seed) vapply(get_jags_inits(NULL, seed, 4, "post"), `[[`, integer(1), ".RNG.seed"))
  expect_length(unique(unlist(chain_seeds)), 80)

  # The user's random number stream is unchanged
  set.seed(1)
  expected = stats::runif(1)
  set.seed(1)
  get_jags_inits(NULL, 42, 3, "post")
  expect_identical(stats::runif(1), expected)
})


test_that("get_median_inits uses numeric priors and truncation", {
  prior = list(a = "dnorm(1, 2)", b = "dt(0, 3, 3) T(0, )", c = "dbeta(1, 3)", d = "dnorm(0, b)", e_inverse = "dgamma(2, 3)")
  expect_equal(
    get_median_inits(lapply(prior, parse_prior)),
    list(a = 1, b = 3 * stats::qt(0.75, 3), c = stats::qbeta(0.5, 1, 3), e_inverse = stats::qgamma(0.5, 2, 3))
  )
})


test_that("get_jags_inits handles shared and chain-specific initial values", {
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
  default_seeds = vapply(get_jags_inits(NULL, 42, 2, "post"), `[[`, integer(1), ".RNG.seed")
  expect_equal(vapply(seeded, `[[`, integer(1), ".RNG.seed"), default_seeds)

  # Chain-specific inits attach RNG settings to each chain
  chain_inits = list(
    list(Intercept_1 = -1, .RNG.seed = 100),
    list(Intercept_1 = 1, .RNG.seed = 200)
  )
  chain_seeded = get_jags_inits(chain_inits, 42, 2, "post")
  expect_equal(vapply(chain_seeded, `[[`, numeric(1), "Intercept_1"), c(-1, 1))
  expect_true(all(
    vapply(chain_seeded, `[[`, character(1), ".RNG.name") ==
      "base::Wichmann-Hill"
  ))
  expect_equal(vapply(chain_seeded, `[[`, integer(1), ".RNG.seed"), default_seeds)

  # Mismatched length throws informative error
  expect_error(
    get_jags_inits(
      list(list(Intercept_1 = -1)),
      42, 2, "post"
    ),
    "`inits` must have length equal to `chains` (2)",
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
    mcp(seed_model, seed_data, par_x = "x", sample = FALSE, warmup = 1.5),
    "single integerish value",
    fixed = TRUE
  )
})


test_that("mcp seed reproduces prior and posterior draws", {
  old_plan = future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  fit_seeded = function(seed) {
    expect_warning({
      fit = mcp(
        seed_model,
        seed_data,
        par_x = "x",
        sample = "both",
        chains = 2,
        warmup = 20,
        iter = 30,
        inits = list(Intercept_1 = 0),
        seed = seed,
        diagnostics = FALSE,
        quiet = TRUE
      )
    }, "Adaptation incomplete", fixed = TRUE)
    fit
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


test_that("set.seed() reproduces sequential sampling without `seed`", {
  old_plan = future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  fit_set_seed = function() {
    set.seed(123)
    suppressWarnings(mcp(seed_model, seed_data, par_x = "x", chains = 2, warmup = 20, iter = 30, diagnostics = FALSE, quiet = TRUE))
  }
  expect_identical(.subset2(fit_set_seed(), "mcmc_post"), .subset2(fit_set_seed(), "mcmc_post"))
})


test_that("chain-specific inits succeed when sequential change points are present", {
  data = data.frame(x = 1:20, y = c(rnorm(10, 1), rnorm(10, 5)))
  inits_chains = list(
    list(Intercept_1 = 0.5),
    list(Intercept_1 = 1.5)
  )

  # Sequential execution
  expect_no_error(
    suppressWarnings(
      mcp(
        list(y ~ 1, ~ 1),
        data,
        par_x = "x",
        inits = inits_chains,
        chains = 2,
        iter = 5,
        warmup = 5,
        quiet = TRUE,
        diagnostics = FALSE
      )
    )
  )

  # Parallel execution (verifies internal seed generation accepts chain-specific inits)
  old_plan = future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  expect_no_error(
    suppressWarnings(
      mcp(
        list(y ~ 1, ~ 1),
        data,
        par_x = "x",
        inits = inits_chains,
        chains = 2,
        iter = 5,
        warmup = 5,
        quiet = TRUE,
        diagnostics = FALSE
      )
    )
  )
})


test_that("inits only accept directly sampled parameters", {
  model = list(y ~ 1 + x, ~ cp(1 + (1 | id)) + 0 + x, ~ 1)
  data = data.frame(x = 1:40, y = rnorm(40), id = rep(1:4, 10))
  expect_no_error(mcp(model, data, inits = list(x_1 = 0, cp_1_id_sd = 1, Intercept_3 = 0), sample = FALSE))
  expect_no_error(mcp(model, data, inits = list(list(cp_2 = 30), list(cp_2 = 35)), chains = 2, sample = FALSE))
  expect_error(mcp(model, data, inits = list(Intercept_1 = 0), sample = FALSE), "Got `Intercept_1`")  # Reparameterized
  expect_error(mcp(model, data, inits = list(cp_1_id = 0), sample = FALSE), "Got `cp_1_id`")  # Deterministic
  expect_error(mcp(model, data, inits = list(list(cp_2 = 30), list(Intercept_1_mid_ = 0)), chains = 2, sample = FALSE), "Got `Intercept_1_mid_`")  # Sampler node
})
