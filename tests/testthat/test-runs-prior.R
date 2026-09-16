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


testthat::test_that("Change-point priors outside the observed range are allowed", {
  test_runs(prior_model, sample = FALSE, prior = list(
    cp_1 = "dunif(-100, -90)",
    cp_2 = "dnorm(100, 20) T(100, 110)"
  ))

  outside_range_model = list(
    y ~ 1 + x,
    ~ 1 + x,
    ~ 0
  )
  test_runs(outside_range_model, prior = list(
    cp_1 = "dunif(-100, -90)",
    cp_2 = "dnorm(100, 20) T(100, 110)"
  ))
})


testthat::test_that("Dirichlet change point priors use a common alpha", {
  fit = mcp(
    prior_model,
    data = data_gauss,
    prior = list(cp_1 = "dirichlet(0.5)", cp_2 = "dirichlet(0.5)"),
    par_x = "x",
    sample = FALSE,
    quiet = TRUE
  )
  testthat::expect_match(fit$jags_code, "dbeta(0.5, 1)", fixed = TRUE)
  fit_10 = mcp(
    prior_model,
    data = data_gauss,
    prior = list(cp_1 = "dirichlet(10)", cp_2 = "dirichlet(10)"),
    par_x = "x",
    sample = FALSE,
    quiet = TRUE
  )
  testthat::expect_match(fit_10$jags_code, "dbeta(10, 20)", fixed = TRUE)
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


testthat::test_that("parse_prior_call parses prior calls, arguments, and truncation", {
  # Ordinary distributions and nested calls
  testthat::expect_equal(parse_prior_call("dnorm(0, 1)"), list(name = "dnorm", args = c("0", "1")))
  testthat::expect_equal(parse_prior_call("dt(0, 2.5, 3)"), list(name = "dt", args = c("0", "2.5", "3")))
  testthat::expect_equal(parse_prior_call("dunif(min(x), max(x))"), list(name = "dunif", args = c("min(x)", "max(x)")))
  testthat::expect_equal(parse_prior_call("dirichlet(1)"), list(name = "dirichlet", args = "1"))
  testthat::expect_equal(parse_prior_call("dnorm()"), list(name = "dnorm", args = character()))

  # Truncation syntax and missing bounds
  testthat::expect_equal(parse_prior_call("T(, 3)"), list(name = "T", args = c("", "3")))
  testthat::expect_equal(parse_prior_call("T(0, )"), list(name = "T", args = c("0", "")))
  testthat::expect_equal(parse_prior_call("T(,)"), list(name = "T", args = c("", "")))

  # Expressions with commas or parentheses inside symbols/brackets
  testthat::expect_equal(parse_prior_call("dnorm(`a,b`, 1)"), list(name = "dnorm", args = c("`a,b`", "1")))
  testthat::expect_equal(parse_prior_call("dnorm(`a)b`, 1)"), list(name = "dnorm", args = c("`a)b`", "1")))
  testthat::expect_equal(parse_prior_call("dnorm(a[1, 2], 1)"), list(name = "dnorm", args = c("a[1, 2]", "1")))

  # Non-calls and malformed inputs
  testthat::expect_null(parse_prior_call(NULL))
  testthat::expect_null(parse_prior_call(NA_character_))
  testthat::expect_null(parse_prior_call(""))
  testthat::expect_null(parse_prior_call("   "))
  testthat::expect_null(parse_prior_call("5"))
  testthat::expect_null(parse_prior_call("-0.5"))
  testthat::expect_null(parse_prior_call("Intercept_1"))
  testthat::expect_null(parse_prior_call("x_1 + x_2"))
  testthat::expect_null(parse_prior_call("(dnorm(1, 2))"))
  testthat::expect_null(parse_prior_call("dnorm(1, "))
  testthat::expect_null(parse_prior_call("dnorm(1, 2) T(0, )"))

  # Named arguments are explicitly rejected
  testthat::expect_error(parse_prior_call("dnorm(sd = 1, mean = 10)"), "Named arguments are not supported")
  testthat::expect_error(parse_prior_call("T(lower = 0, )"), "Named arguments are not supported")
})


testthat::test_that("default log-count priors use normal population slopes/contrasts and half-normal group SDs", {
  d = data.frame(
    x = 1:10,
    cat = factor(rep(c("A", "B"), 5)),
    id = factor(rep(1:5, each = 2)),
    y = 1:10
  )
  for (fam in list(poisson(), negbinomial())) {
    fit = mcp(
      list(y ~ 1 + x + cat + (1 + x + cat || id)),
      data = d,
      family = fam,
      sample = FALSE
    )
    # Population priors
    testthat::expect_match(fit$prior$Intercept_1, "^dnorm\\(")
    testthat::expect_equal(fit$prior$x_1, "dnorm(0, 0.2777778)")
    testthat::expect_equal(fit$prior$catB_1, "dnorm(0, 2.5)")

    # Group SD priors
    testthat::expect_equal(fit$prior$Intercept_1_id_sd, "dnorm(0, 2.5) T(0, )")
    testthat::expect_equal(fit$prior$x_1_id_sd, "dnorm(0, 0.2777778) T(0, )")
    testthat::expect_equal(fit$prior$catB_1_id_sd, "dnorm(0, 2.5) T(0, )")
  }

  # Modeled shape in negative-binomial
  fit_shape = mcp(
    list(y ~ 1 + x + shape(1 + x + cat + (1 + x + cat || id))),
    data = d,
    family = negbinomial(),
    sample = FALSE
  )
  testthat::expect_equal(fit_shape$prior$shape_1, "dnorm(0, 2.5)")
  testthat::expect_equal(fit_shape$prior$shape_x_1, "dnorm(0, 0.2777778)")
  testthat::expect_equal(fit_shape$prior$shape_catB_1, "dnorm(0, 2.5)")
  testthat::expect_equal(fit_shape$prior$shape_1_id_sd, "dnorm(0, 2.5) T(0, )")
  testthat::expect_equal(fit_shape$prior$shape_x_1_id_sd, "dnorm(0, 0.2777778) T(0, )")
  testthat::expect_equal(fit_shape$prior$shape_catB_1_id_sd, "dnorm(0, 2.5) T(0, )")
})


testthat::test_that("offset adjusts default intercept priors to log-rate", {
  d = data.frame(
    x = 1:30,
    y = c(rep(5, 10), rep(10, 10), rep(20, 10)),
    exposure = rep(10, 30),
    exposure2 = rep(100, 30)
  )

  for (fam in list(poisson(), negbinomial(), gaussian(link = "log"))) {
    # Single offset
    fit = mcp(
      list(y ~ 1 + x + offset(log(exposure))),
      data = d,
      family = fam,
      sample = FALSE
    )
    expected_rate = log(pmax(d$y, 0.1)) - log(d$exposure)
    expected_loc = round(median(expected_rate), 1)
    expected_scale = max(2.5, round(mad(expected_rate), 1))
    testthat::expect_equal(fit$prior$Intercept_1, paste0("dnorm(", expected_loc, ", ", expected_scale, ")"))

    ps = prior_summary(fit, verbose = TRUE)
    int_row = ps[ps$parameter == "Intercept_1", ]
    testthat::expect_match(int_row$rule, "- offset", fixed = TRUE)
    testthat::expect_equal(int_row$description, "Robustly centered log-rate intercept with a minimum scale of 2.5")

    # A later segment needs its own offset declaration.
    fit_multi = mcp(
      list(
        y ~ 1 + x + offset(log(exposure)),
        1 ~ 1 + x,
        1 ~ 1 + offset(0)
      ),
      data = d,
      family = fam,
      sample = FALSE
    )
    expected_count = log(pmax(d$y, 0.1))
    expected_count_loc = round(median(expected_count), 1)
    expected_count_scale = max(2.5, round(mad(expected_count), 1))
    testthat::expect_equal(fit_multi$prior$Intercept_1, paste0("dnorm(", expected_loc, ", ", expected_scale, ")"))
    testthat::expect_equal(fit_multi$prior$Intercept_2, paste0("dnorm(", expected_count_loc, ", ", expected_count_scale, ")"))

    # Segment 3 reverted to log-count / log-mean
    testthat::expect_equal(fit_multi$prior$Intercept_3, paste0("dnorm(", expected_count_loc, ", ", expected_count_scale, ")"))

    ps_multi = prior_summary(fit_multi, verbose = TRUE)
    testthat::expect_equal(ps_multi$description[ps_multi$parameter == "Intercept_1"], "Robustly centered log-rate intercept with a minimum scale of 2.5")
    testthat::expect_match(ps_multi$description[ps_multi$parameter == "Intercept_2"], "Robustly centered log-(count|mean) intercept with a minimum scale of 2.5")
    testthat::expect_match(ps_multi$description[ps_multi$parameter == "Intercept_3"], "Robustly centered log-(count|mean) intercept with a minimum scale of 2.5")

    # Multiple distinct offsets in different segments
    fit_diff = mcp(
      list(
        y ~ 1 + x + offset(log(exposure)),
        1 ~ 1 + x + offset(log(exposure2))
      ),
      data = d,
      family = fam,
      sample = FALSE
    )
    expected_rate2 = log(pmax(d$y, 0.1)) - log(d$exposure2)
    expected_loc2 = round(median(expected_rate2), 1)
    expected_scale2 = max(2.5, round(mad(expected_rate2), 1))
    testthat::expect_equal(fit_diff$prior$Intercept_1, paste0("dnorm(", expected_loc, ", ", expected_scale, ")"))
    testthat::expect_equal(fit_diff$prior$Intercept_2, paste0("dnorm(", expected_loc2, ", ", expected_scale2, ")"))

    ps_diff = prior_summary(fit_diff, verbose = TRUE)
    testthat::expect_match(ps_diff$rule[ps_diff$parameter == "Intercept_1"], "- offset_1", fixed = TRUE)
    testthat::expect_match(ps_diff$rule[ps_diff$parameter == "Intercept_2"], "- offset_2", fixed = TRUE)
  }
})


testthat::test_that("offsets on other distributional parameters do not overwrite mu intercept priors", {
  d = data.frame(
    x = 1:10,
    y = 10,
    exposure = 10,
    exposure2 = 100
  )

  fit1 = mcp(
    list(y ~ 1 + offset(log(exposure))),
    data = d,
    family = negbinomial(),
    par_x = "x",
    sample = FALSE
  )
  testthat::expect_equal(fit1$prior$Intercept_1, "dnorm(0, 2.5)")

  fit2 = mcp(
    list(y ~ 1 + offset(log(exposure)) + shape(1 + offset(log(exposure2)))),
    data = d,
    family = negbinomial(),
    par_x = "x",
    sample = FALSE
  )
  testthat::expect_equal(fit2$prior$Intercept_1, "dnorm(0, 2.5)")
  testthat::expect_equal(fit2$prior$shape_1, "dnorm(0, 2.5)")
})


testthat::test_that("gaussian(link = 'log') aligns with log-link model default priors", {
  d = data.frame(
    time = 1:20,
    y = exp(seq(1, 4, length.out = 20)),
    group = rep(c("A", "B"), 10),
    id = rep(1:5, each = 4)
  )

  fit = mcp(
    list(y ~ 1 + time + group + (1 + time || id)),
    data = d,
    family = gaussian(link = "log"),
    sample = FALSE
  )

  ps = prior_summary(fit, verbose = TRUE)
  int_row = ps[ps$parameter == "Intercept_1", ]
  testthat::expect_match(int_row$prior, "^normal\\(mean = [0-9.]+, sd = [0-9.]+\\)$")
  testthat::expect_match(int_row$rule, "log\\(pmax\\(y, 0.1\\)\\)")

  # Categorical contrast on log link uses dnorm(0, 2.5)
  group_row = ps[ps$parameter == "groupB_1", ]
  testthat::expect_equal(group_row$prior, "normal(mean = 0, sd = 2.5)")
  testthat::expect_equal(group_row$rule, "normal(mean = 0, sd = 2.5)")
  testthat::expect_equal(fit$prior$groupB_1, "dnorm(0, 2.5)")

  # Slope on log link uses dnorm(0, 2.5 / predictor_scale())
  time_row = ps[ps$parameter == "time_1", ]
  time_span = diff(range(d$time))
  expected_slope_sd = format_prior_number(2.5 / time_span)
  testthat::expect_equal(time_row$prior, paste0("normal(mean = 0, sd = ", expected_slope_sd, ")"))
  testthat::expect_equal(fit$prior$time_1, paste0("dnorm(0, ", expected_slope_sd, ")"))

  # Group SDs are half-normals
  id_int_row = ps[ps$parameter == "Intercept_1_id_sd", ]
  testthat::expect_equal(id_int_row$prior, "normal(mean = 0, sd = 2.5)")
  testthat::expect_equal(id_int_row$bounds, "[0, Inf]")

  id_time_row = ps[ps$parameter == "time_1_id_sd", ]
  testthat::expect_equal(id_time_row$prior, paste0("normal(mean = 0, sd = ", expected_slope_sd, ")"))
  testthat::expect_equal(id_time_row$bounds, "[0, Inf]")

  # Sigma_1 uses response-scale half-Student-t
  sig_row = ps[ps$parameter == "sigma_1", ]
  testthat::expect_match(sig_row$prior, "^student_t\\(df = 3, location = 0, scale = [0-9.]+\\)$")
  testthat::expect_equal(sig_row$bounds, "[0.001, Inf]")
})
