######################
# TEST CHANGE POINTS #
######################
bad_cps = list(
  list(y ~ x,
       0 ~ 1),  # Needs changepoint stuff
  list(y ~ x,
       q ~ 1),  # Slope not allowed for changepoint
  list(y ~ 1,
       (goat|id) ~ 1),  # No varying slope allowed
  list(y ~ 1,
       y ~ ~ 1),  # Needs to be explicit if y is defined
  list(y ~ 1,
       1 + (1|bad_id) ~ 1),  # decimal group
  list(y ~ 1,
       1 + (1|id) ~ 1,
       1 + (1|ok_id_integer) ~ 1),  # CP grouping factors must agree
  list(y ~ 1,
       1 + (1|id) ~ 1,
       1 + (1|ok_id_factor) ~ 1)  # CP grouping factors must agree
)

test_bad(bad_cps)




good_cps = list(
  list(y ~ 0 + x,  # Regular cp
       1 ~ 1),
  list(y ~ 1,  # Implicit cp
       ~ 1,
       ~ 0),
  list(y ~ 0,  # Varying
       1 + (1|id) ~ 1),
  list(y ~ 0,  # Chained varying cps
       y ~ 1 ~ 1,
       1 + (1|id) ~ 0,
       1 + (1|id) ~ 0,
       ~ x),
  list(y ~ 1,
       (1|id) ~ 0),  # Intercept is implicit. I don't like it, but OK.
  list(y ~ 1,
       1 + (1|id) ~ 1,
       1 + (1|id) ~ 1,
       1 + (1|id) ~ 1)
)

test_good(good_cps)



# Test detection of par_x
bad_par_x = list(
  list(y ~ 0),  # no par_x
  list(y ~ 1),  # no par_x
  list(y ~ 1 + bad_x_char,
       ~ 0 + bad_x_factor),  # only invalid par_x
  list(y ~ bad_x_char),  # Has to be continuous
  list(y ~ bad_x_factor),  # Has to be continuous
  list(y ~ x + ok_x),  # multiple in one segment
  list(y ~ x,
       ~ ok_x)  # Multiple in two segments
)

test_bad(bad_par_x, par_x = NULL)


good_par_x = list(
  list(y ~ 1 + x),
  list(y ~ 0,
       ~ 1,
       ~ 0 + ok_x)
)

test_good(good_par_x, par_x = NULL)



################
# TEST WEIGHTS #
################
bad_weights = list(
  list(y + weights(weights_ok) ~ 1),  # weights added
  list(weights(y) ~ 1),  # just wrong :-)
  list(y | weights_ok ~ 1),  # Has to be weights(weights_ok)
  list(y | weights(weights_bad) ~ 1),  # Bad weights
  list(y | weights(weights_ok) ~ 1,
       y | weights(weights_bad) ~ 1 ~ 1)  # Different weights
)

test_bad(bad_weights)

good_weights = list(
  list(y | weights(weights_ok) ~ 1),  # Regular
  list(y | weights(weights_ok) ~ 1,
       ~ 1 + x + I(x^2),
       1 + (1|id) ~ 1)  # With multiple segments and functions and varying
)

test_good(good_weights)


test_that("Gaussian JAGS weights implement a likelihood power", {
  weighted_data = data.frame(x = 1:4, y = c(0, NA, 2, 3), w = c(0.5, 2, 3, 1))
  expect_message({
    fit = mcp(
      list(y | weights(w) ~ 1), weighted_data,
      par_x = "x", sample = FALSE
    )
  }, "NA values detected in 'y'")
  tables = get_fit_model_tables(fit)
  jags_data = get_jags_data(
    fit$data, fit$family, tables$segments, tables$predictors,
    tables$group_effects, fit$jags_code
  )

  expect_equal(jags_data$response_observed_, c(1, 0, 1, 1))
  expect_equal(jags_data$likelihood_zero_, numeric(4))
  expect_equal(jags_data$w, weighted_data$w)
  expect_match(fit$jags_code, "likelihood_weight_[i_] = 1 + response_observed_[i_] * (w[i_] - 1)", fixed = TRUE)
  expect_match(fit$jags_code, "likelihood_phi_[i_] = response_observed_[i_] * max(0, (likelihood_weight_[i_] - 1) * log(sigma_[i_]) + abs(likelihood_weight_[i_] - 1) * 20)", fixed = TRUE)
  expect_match(fit$jags_code, "likelihood_zero_[i_] ~ dpois(likelihood_phi_[i_])", fixed = TRUE)

  # The JAGS normal density plus its Poisson zeros correction is the desired
  # powered normal density up to a constant that depends only on the weight.
  grid = expand.grid(mu = c(-1, 0.5), sigma = c(0.4, 2), w = c(0.5, 2, 3))
  y = 1.2
  phi = (grid$w - 1) * log(grid$sigma) + abs(grid$w - 1) * 20
  jags_kernel = stats::dnorm(y, grid$mu, grid$sigma / sqrt(grid$w), log = TRUE) +
    stats::dpois(0, lambda = phi, log = TRUE)
  power_kernel = grid$w * stats::dnorm(y, grid$mu, grid$sigma, log = TRUE)
  kernel_difference = jags_kernel - power_kernel
  expect_equal(
    kernel_difference,
    0.5 * log(grid$w) + 0.5 * (grid$w - 1) * log(2 * pi) - abs(grid$w - 1) * 20
  )
})


test_that("Weighted Gaussian samples on ordinary large weights without numerical underflow", {
  df_gauss = data.frame(x = 1:4, y = c(1.2, 1.4, 0.9, 1.1), w = 1000)
  # Previously failed during compile/initialization with:
  # 'Error in node likelihood_zero_[1]: Invalid parent values'
  fit = mcp(
    list(y | weights(w) ~ 1),
    data = df_gauss,
    family = gaussian(),
    par_x = "x",
    prior = list(Intercept_1 = "dnorm(0, 1)", sigma_1 = "dunif(0.001, 10)"),
    diagnostics = FALSE,
    quiet = TRUE
  )
  expect_true(is.list(fit$model))
  draws = as.matrix(coda::as.mcmc(fit))
  expect_equal(mean(draws[, "Intercept_1"]), 1.15, tolerance = 0.05)
  expect_equal(mean(draws[, "sigma_1"]), 0.18, tolerance = 0.05)
})
