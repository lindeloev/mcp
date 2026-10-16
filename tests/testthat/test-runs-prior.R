###############
# TEST PRIORS #
###############
prior_model = list(
  y ~ 1 + x,
  ~ cp(1 + (1|id)) + 1 + x,
  ~ 0
)

bad_prior = list(
  list(
    cp_1 = "dirichlet(1)",  # Has to be all-dirichlet
    cp_2 = "dnorm(3, 10)"
  ),
  list(
    cp_1 = "dirichlet(0.5)",
    cp_2 = "dirichlet(0.5)"  # alpha has to be >= 1
  ),
  list(
    cp_1 = "dunif(min(x), cp_2)"  # Limits can only refer to the previous change point
  ),
  list(
    cp_2 = "dnorm(1, 1) T(cp_1 - 0.5, )"  # ... and only as a minimum segment length
  ),
  list(
    cp_2 = "dunif(0, 3) T(0, 2)"  # JAGS cannot truncate dunif()
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
  ),
  list(
    cp_1 = "dnorm(0, 1)",
    cp_2 = "cp_1 + 1"  # Fixed change point relative to another
  ),
  list(
    cp_1 = "dunif(-1, 3)",  # Identical priors
    cp_2 = "dunif(-1, 3)"
  ),
  list(
    cp_2 = "dunif(cp_1 + 0.5, max(x))"  # Minimum segment length
  ),
  list(
    cp_2 = "dnorm(2, 1) T(cp_1 + 0.5, )"  # Minimum segment length with truncation
  ),
  list(
    cp_2 = 2  # Fixed change point
  ),
  list(
    cp_2 = "dnorm(cp_1 + 1, 1)"  # Conditional prior
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
    prior = list(cp_1 = "dirichlet(1)", cp_2 = "dirichlet(1)"),
    par_x = "x",
    sample = FALSE,
    quiet = TRUE
  )
  testthat::expect_match(fit$jags_code, "cp_2 ~ dunif(cp_0, cp_3)", fixed = TRUE)
  testthat::expect_false(grepl("cp_weight_", fit$jags_code, fixed = TRUE))
  fit_10 = mcp(
    prior_model,
    data = data_gauss,
    prior = list(cp_1 = "dirichlet(10)", cp_2 = "dirichlet(10)"),
    par_x = "x",
    sample = FALSE,
    quiet = TRUE
  )
  testthat::expect_match(fit_10$jags_code, "cp_weight_ ~ dpois(-9 * (log((cp_1 - cp_0) / (cp_3 - cp_0)) + ", fixed = TRUE)
  testthat::expect_error(
    mcp(prior_model, data = data_gauss, prior = list(cp_1 = "dirichlet(2)", cp_2 = "dirichlet(3)"), par_x = "x", sample = FALSE, quiet = TRUE),
    "same alpha"
  )
  testthat::expect_error(
    mcp(prior_model, data = data_gauss, prior = list(cp_1 = "dirichlet(0.5)", cp_2 = "dirichlet(0.5)"), par_x = "x", sample = FALSE, quiet = TRUE),
    "finite alpha >= 1"
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


testthat::test_that("Minimum segment lengths go in cp_order_", {
  fit = mcp(list(y ~ 1, ~ 1, ~ 1), data_gauss, par_x = "x", prior = list(cp_2 = "dnorm(2, 1) T(cp_1 + 0.5, )"), sample = FALSE)
  testthat::expect_equal(fit$prior$cp_2, "dnorm(2, 1) T(cp_1 + 0.5, max(x))")
  testthat::expect_match(fit$jags_code, "cp_2 ~ dnorm(2, 1/(1)^2) T(cp_0,cp_3)", fixed = TRUE)
  testthat::expect_match(fit$jags_code, "cp_order_ ~ dbern(step(cp_2 - cp_1 - 0.5))", fixed = TRUE)

  # A lower bound at the previous change point is a minimum segment length of 0, i.e., just the ordering
  testthat::expect_warning(
    fit <- mcp(list(y ~ 1, ~ 1, ~ 1), data_gauss, par_x = "x", prior = list(cp_2 = "dunif(cp_1, max(x))"), sample = FALSE),
    "as the lower bound is already embodied in the mcp model"
  )
  testthat::expect_match(fit$jags_code, "cp_2 ~ dunif(cp_0, cp_3)", fixed = TRUE)
  testthat::expect_match(fit$jags_code, "cp_order_ ~ dbern(step(cp_2 - cp_1))", fixed = TRUE)
})


testthat::test_that("Change points keep their specified priors, restricted to ordered change points", {
  data = data.frame(x = seq(0, 1, length.out = 20), y = 0)
  set.seed(123)
  n = 1e5
  rtnorm = function(mean, sd) stats::qnorm(stats::runif(n, stats::pnorm(0, mean, sd), stats::pnorm(1, mean, sd)), mean, sd)
  restrict = function(cp_1, cp_2, min_length = 0) cbind(cp_1, cp_2)[cp_2 - cp_1 >= min_length, ]
  cases = list(
    list(prior = list(cp_1 = "dnorm(0.6, 0.2)", cp_2 = "dnorm(0.4, 0.2)"), r_draws = restrict(rtnorm(0.6, 0.2), rtnorm(0.4, 0.2))),
    list(prior = list(cp_2 = "dunif(cp_1 + 0.3, max(x))"), r_draws = restrict(stats::runif(n), stats::runif(n), 0.3))
  )
  for (case in cases) {
    fit = suppressWarnings(mcp(list(y ~ 1, ~ 1, ~ 1), data, par_x = "x", prior = case$prior, sample = "prior", chains = 1, iter = 6000, warmup = 100, seed = 1, quiet = TRUE))
    jags_draws = as.matrix(.subset2(fit, "mcmc_prior")[[1]])[, c("cp_1", "cp_2")]
    testthat::expect_lt(max(abs(colMeans(jags_draws) - colMeans(case$r_draws))), 0.025)
  }
})


testthat::test_that("Change point inits are ordered and within the priors, around fixed and user-supplied values", {
  fit = mcp(list(y ~ 1, ~ 1, ~ 1, ~ 1), data.frame(x = 0:100, y = 0), par_x = "x", sample = FALSE,
    prior = list(cp_2 = "dunif(80, 90)", cp_3 = "dnorm(50, 5)"))
  jags_prior = function(fit) as.list(stats::setNames(fit$.internal$prior_table$value, fit$.internal$prior_table$parameter))  # With cp_0 and cp_K
  testthat::expect_equal(get_cp_inits(jags_prior(fit), c(0, 100)), list(cp_1 = 25, cp_2 = 81, cp_3 = 82))
  testthat::expect_equal(get_cp_inits(jags_prior(fit), c(0, 100), list(cp_2 = 85)), list(cp_1 = 25, cp_3 = 86))
  testthat::expect_error(get_cp_inits(jags_prior(fit), c(0, 100), list(cp_3 = 70)), "no room for ordered change points")

  fit = mcp(list(y ~ 1, ~ 1, ~ 1), data.frame(x = 0:100, y = 0), par_x = "x", sample = FALSE, prior = list(cp_2 = 20))
  testthat::expect_equal(get_cp_inits(jags_prior(fit), c(0, 100)), list(cp_1 = 19))
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
  x_scale = format_prior_number(2.5 / stats::sd(d$x))
  cat_scale = format_prior_number(2.5 / stats::sd(d$cat == "B"))
  for (fam in list(poisson(), negbinomial())) {
    fit = mcp(
      list(y ~ 1 + x + cat + (1 + x + cat || id)),
      data = d,
      family = fam,
      sample = FALSE
    )
    # Population priors
    testthat::expect_match(fit$prior$Intercept_1, "^dnorm\\(")
    expect_prior(fit$prior$x_1, paste0("dnorm(0, ", x_scale, ")"))
    expect_prior(fit$prior$catB_1, paste0("dnorm(0, ", cat_scale, ")"))

    # Group SD priors
    testthat::expect_equal(fit$prior$Intercept_1_id_sd, "dnorm(0, 2.5) T(0, )")
    testthat::expect_equal(fit$prior$x_1_id_sd, paste0("dnorm(0, ", x_scale, ") T(0, )"))
    testthat::expect_equal(fit$prior$catB_1_id_sd, "dnorm(0, 2.5) T(0, )")
  }

  # Modeled shape in negative-binomial
  fit_shape = mcp(
    list(y ~ 1 + x + shape(1 + x + cat + (1 + x + cat || id))),
    data = d,
    family = negbinomial(),
    sample = FALSE
  )
  testthat::expect_equal(fit_shape$prior$shape_Intercept_1, "dnorm(-shape_x_1 * 1, 2.5)")
  testthat::expect_equal(fit_shape$prior$shape_x_1, paste0("dnorm(0, ", x_scale, ")"))
  testthat::expect_equal(fit_shape$prior$shape_catB_1, "dnorm(0, 2.5)")
  testthat::expect_equal(fit_shape$prior$shape_Intercept_1_id_sd, "dnorm(0, 2.5) T(0, )")
  testthat::expect_equal(fit_shape$prior$shape_x_1_id_sd, paste0("dnorm(0, ", x_scale, ") T(0, )"))
  testthat::expect_equal(fit_shape$prior$shape_catB_1_id_sd, "dnorm(0, 2.5) T(0, )")
})


testthat::test_that("offset adjusts default log-link intercept priors to log(y) - offset", {
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
    testthat::expect_equal(fit$prior$Intercept_1, paste0("dnorm(", expected_loc, " - x_1 * 1, ", expected_scale, ")"))

    ps = prior_summary(fit, verbose = TRUE)
    int_row = ps[ps$parameter == "Intercept_1", ]
    testthat::expect_match(int_row$rule, "- offset", fixed = TRUE)
    testthat::expect_match(int_row$description, "^Robustly centered log-(count|mean) intercept with a minimum scale of 2.5; adjusted for the offset; on the value at min\\(x\\)$")

    # A later segment must include its own offset.
    fit_multi = mcp(
      list(
        y ~ 1 + x + offset(log(exposure)),
        ~ 1 + x + offset(0),
        ~ 1 + offset(0)
      ),
      data = d,
      family = fam,
      sample = FALSE
    )
    expected_count = log(pmax(d$y, 0.1))
    expected_count_loc = round(median(expected_count), 1)
    expected_count_scale = max(2.5, round(mad(expected_count), 1))
    testthat::expect_equal(fit_multi$prior$Intercept_1, paste0("dnorm(", expected_loc, " - x_1 * 1, ", expected_scale, ")"))
    testthat::expect_equal(fit_multi$prior$Intercept_2, paste0("dnorm(", expected_count_loc, ", ", expected_count_scale, ")"))

    # Segment 3 reverted to log-count / log-mean
    testthat::expect_equal(fit_multi$prior$Intercept_3, paste0("dnorm(", expected_count_loc, ", ", expected_count_scale, ")"))

    ps_multi = prior_summary(fit_multi, verbose = TRUE)
    testthat::expect_match(ps_multi$description[ps_multi$parameter == "Intercept_1"], "^Robustly centered log-(count|mean) intercept with a minimum scale of 2.5; adjusted for the offset; on the value at min\\(x\\)$")
    testthat::expect_match(ps_multi$description[ps_multi$parameter == "Intercept_2"], "Robustly centered log-(count|mean) intercept with a minimum scale of 2.5")
    testthat::expect_match(ps_multi$description[ps_multi$parameter == "Intercept_3"], "Robustly centered log-(count|mean) intercept with a minimum scale of 2.5")

    # Multiple distinct offsets in different segments
    fit_diff = mcp(
      list(
        y ~ 1 + x + offset(log(exposure)),
        ~ 1 + x + offset(log(exposure2))
      ),
      data = d,
      family = fam,
      sample = FALSE
    )
    expected_rate2 = log(pmax(d$y, 0.1)) - log(d$exposure2)
    expected_loc2 = round(median(expected_rate2), 1)
    expected_scale2 = max(2.5, round(mad(expected_rate2), 1))
    testthat::expect_equal(fit_diff$prior$Intercept_1, paste0("dnorm(", expected_loc, " - x_1 * 1, ", expected_scale, ")"))
    testthat::expect_equal(fit_diff$prior$Intercept_2, paste0("dnorm(", expected_loc2, ", ", expected_scale2, ")"))

    ps_diff = prior_summary(fit_diff, verbose = TRUE)
    testthat::expect_match(ps_diff$rule[ps_diff$parameter == "Intercept_1"], "- offset_1", fixed = TRUE)
    testthat::expect_match(ps_diff$rule[ps_diff$parameter == "Intercept_2"], "- offset_2", fixed = TRUE)
  }
})


testthat::test_that("offset adjusts identity-link default intercept priors to y - offset", {
  d = data.frame(
    x = 1:30,
    y = 1000 + c(rep(5, 10), rep(10, 10), rep(20, 10)),
    baseline = 1000 + rep(0:2, 10),
    baseline2 = 990
  )

  for (fam in list(gaussian(), poisson(link = "identity"))) {
    fit = mcp(
      list(y ~ 1 + offset(baseline), ~ 1 + offset(baseline2)),
      data = d,
      par_x = "x",
      family = fam,
      sample = FALSE
    )
    for (i in 1:2) {
      adjusted = d$y - d[[c("baseline", "baseline2")[i]]]
      expect_prior(split_prior_truncation(fit$prior[[paste0("Intercept_", i)]])$distribution, paste0("dnorm(", format_prior_number(mean(adjusted)), ", ", format_prior_number(2.5 * sd(adjusted)), ")"))
    }

    ps = prior_summary(fit, verbose = TRUE)
    testthat::expect_match(ps$rule[ps$parameter == "Intercept_1"], "mean(y - offset_1)", fixed = TRUE)

    # Positive means (Poisson) truncate the intercept at -offset rather than 0
    if (fam$family == "poisson") {
      testthat::expect_match(fit$prior$Intercept_1, paste0("T(-", min(d$baseline), ", )"), fixed = TRUE)
      testthat::expect_match(fit$prior$Intercept_2, paste0("T(-", min(d$baseline2), ", )"), fixed = TRUE)
    }
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
  testthat::expect_equal(fit2$prior$shape_Intercept_1, "dnorm(0, 2.5)")
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
  testthat::expect_match(int_row$prior, "^dnorm\\([0-9.]+ - time_1 \\* 1, [0-9.]+\\)$")
  testthat::expect_match(int_row$rule, "log(ifelse(y <= 0, 0.1, y))", fixed = TRUE)

  # Categorical contrast on log link uses dnorm(0, 2.5 / sd(dummy))
  group_row = ps[ps$parameter == "groupB_1", ]
  expected_group_sd = format_prior_number(2.5 / stats::sd(d$group == "B"))
  testthat::expect_equal(group_row$prior, format_prior_display(paste0("dnorm(0, ", expected_group_sd, ")")))
  expect_prior(fit$prior$groupB_1, paste0("dnorm(0, ", expected_group_sd, ")"))

  # Slope on log link uses dnorm(0, 2.5 / sd(x))
  time_row = ps[ps$parameter == "time_1", ]
  expected_slope_sd = format_prior_number(2.5 / stats::sd(d$time))
  testthat::expect_equal(time_row$prior, format_prior_display(paste0("dnorm(0, ", expected_slope_sd, ")")))
  expect_prior(fit$prior$time_1, paste0("dnorm(0, ", expected_slope_sd, ")"))

  # Group SDs are half-normals
  id_int_row = ps[ps$parameter == "Intercept_1_id_sd", ]
  testthat::expect_equal(id_int_row$prior, "dnorm(0, 2.5) T(0, )")
  testthat::expect_equal(id_int_row$bounds, "[0, Inf]")

  id_time_row = ps[ps$parameter == "time_1_id_sd", ]
  testthat::expect_equal(id_time_row$prior, format_prior_display(paste0("dnorm(0, ", expected_slope_sd, ") T(0, )")))
  testthat::expect_equal(id_time_row$bounds, "[0, Inf]")

  # Sigma_1 uses response-scale half-Student-t
  sig_row = ps[ps$parameter == "sigma_1", ]
  testthat::expect_match(sig_row$prior, "^dt\\(0, [0-9.]+, 3\\) T\\(0, \\)$")
  testthat::expect_equal(sig_row$bounds, "[0, Inf]")
})


test_that("legacy prior scales are translated to JAGS parameters", {
  expect_equal(prior_to_jags("dnorm(5, 2)"), "dnorm(5, 1/(2)^2) ")
  expect_equal(prior_to_jags("dt(3, 2, 1)"), "dt(3, 1/(2)^2, 1) ")
  expect_equal(prior_to_jags("dlnorm(3, 2)"), "dlnorm(3, 1/(2)^2) ")
  expect_equal(prior_to_jags("ddexp(3, 2)"), "ddexp(3, 1/(2)) ")
  expect_equal(prior_to_jags("dlogis(3, 2) T(0, )"), "dlogis(3, 1/(2)) T(0,)")
  expect_error(
    prior_to_jags("dcauchy(3, 2)"),
    "use `dt(location, scale, 1)`",
    fixed = TRUE
  )
  expect_warning(
    sd_to_prec("dnorm(5, 2)"),
    "deprecated"
  )
})


test_that("priors are resolved without changing their parameterization", {
  data = data.frame(x = 1:6, y = c(2, 4, 3, 8, 7, 9))
  default_fit = mcp(list(y ~ 1 + x, ~ 1 + x), data, sample = FALSE)
  slope = paste0("dnorm(0, ", format_prior_number(2.5 * stats::sd(data$y) / stats::sd(data$x)), ")")
  intercept = paste0("dnorm(", format_prior_number(mean(data$y)), ", ", format_prior_number(2.5 * stats::sd(data$y)), ")")
  expect_prior(
    unclass(default_fit$prior),
    list(
      cp_1 = "dunif(min(x), max(x))",
      Intercept_1 = sub(",", " - x_1 * 1,", intercept),
      x_1 = slope,
      Intercept_2 = intercept,
      x_2 = slope,
      sigma_1 = "dt(0, 3.7, 3) T(0, )"
    )
  )

  fit = mcp(
    list(y ~ 1 + x, ~ 1 + x),
    data,
    prior = list(
      Intercept_1 = "dt(median(y), mad(y), 3)T(, max(y))",
      x_1 = 5,
      Intercept_2 = "Intercept_1",
      x_2 = "x_1 / (max(x) - min(x))"
    ),
    sample = FALSE
  )

  expect_prior(fit$prior$Intercept_1, paste0("dt(5.5, ", format_prior_number(stats::mad(data$y)), ", 3) T(, 9)"))
  expect_equal(fit$prior$x_2, "x_1/5")
  expect_equal(fit$.internal$prior_format, "jags_string_v1")
  expect_false(grepl("MINX|MAXX|N_CP|LINKY", fit$jags_code))
  expect_match(fit$jags_code, "# User-specified prior", fixed = TRUE)

  compact = prior_summary(fit)
  verbose = prior_summary(fit, verbose = TRUE)
  expect_named(compact, c("parameter", "segment", "dpar", "prior", "bounds"))
  expect_named(
    verbose,
    c("parameter", "segment", "dpar", "prior", "bounds", "rule", "description", "source", "kind")
  )
  expect_equal(
    verbose$kind[match(c("Intercept_1", "x_1", "Intercept_2", "x_2"), verbose$parameter)],
    c("distribution", "constant", "alias", "expression")
  )
  expect_equal(verbose$source[verbose$parameter == "Intercept_1"], "user")

  legacy_fit = NULL
  expect_warning(
    {
      legacy_fit = mcp(
        list(y ~ 1, ~ 1), data,
        par_x = "x", prior = list(cp_1 = "dunif(MINX, MAXX)"), sample = FALSE
      )
    },
    "Deprecated prior data constant"
  )
  expect_equal(legacy_fit$prior$cp_1, "dunif(1, 6)")
})


test_that("resolved priors keep full precision regardless of options(digits)", {
  withr::local_options(digits = 3)
  data = data.frame(x = c(1, 2, 4, 7.3, 8, 9), y = 1:6)
  fit = mcp(list(y ~ 1 + x), data, sample = FALSE, prior = list(Intercept_1 = "dnorm(0, 1/3)"))
  slope_sd = 2.5 * stats::sd(data$y) / stats::sd(data$x)
  expect_equal(fit$prior$Intercept_1, "dnorm(0, 0.33333333333333331)")
  expect_identical(as.numeric(sub("dnorm\\(0, (.*)\\)", "\\1", fit$prior$x_1)), slope_sd)
  expect_match(fit$jags_code, paste0("x_1 ~ dnorm(0, 1/(", format_prior_number(slope_sd), ")^2)"), fixed = TRUE)

  # Only prior_summary() is rounded for display
  ps = prior_summary(fit)
  expect_equal(ps$prior[ps$parameter == "Intercept_1"], "dnorm(0, 0.33333)")
  expect_equal(ps$prior[ps$parameter == "x_1"], paste0("dnorm(0, ", signif(slope_sd, 5), ")"))
})


test_that("change point priors refer to the exact ends of the observed range, cp_0 and cp_K in JAGS", {
  data = data.frame(x = c(0.618632254190743, 2:9, 9.87654321), y = 1:10)
  fit = mcp(
    list(y ~ 1, ~ 1, ~ 1), data, par_x = "x", sample = FALSE,
    prior = list(cp_1 = "dunif(min(x), 5)", cp_2 = "dnorm(7, 1)", Intercept_1 = "dnorm(min(x), 1)")
  )
  expect_equal(fit$prior$cp_1, "dunif(min(x), 5)")
  expect_equal(fit$prior$cp_2, "dnorm(7, 1) T(min(x), max(x))")
  expect_equal(fit$prior$Intercept_1, "dnorm(0.618632254190743, 1)")
  expect_match(fit$jags_code, "cp_0 = min(x)\n  cp_1 ~ dunif(cp_0, 5)  # User-specified prior, restricted to ordered change points\n  cp_2 ~ dnorm(7, 1/(1)^2) T(cp_0,cp_3)  # User-specified prior, truncated to the observed range of x, restricted to ordered change points\n  cp_3 = max(x)\n  cp_order_ ~ dbern(step(cp_2 - cp_1))  # Ordered change points\n", fixed = TRUE)

  # cp_0 and cp_K are JAGS nodes, not parameters
  expect_error(mcp(list(y ~ 1, ~ 1, ~ 1), data, par_x = "x", sample = FALSE, prior = list(cp_1 = "dunif(cp_0, 5)")), "Use `min(x)` and `max(x)`", fixed = TRUE)
  expect_error(mcp(list(y ~ 1, ~ 1, ~ 1), data, par_x = "x", sample = FALSE, prior = list(cp_2 = "dnorm(7, 1) T(, cp_3)")), "`cp_3`", fixed = TRUE)

  # Passing fit$prior back gives the same model
  refit = mcp(list(y ~ 1, ~ 1, ~ 1), data, par_x = "x", sample = FALSE, prior = fit$prior)
  expect_equal(refit$prior, fit$prior)

  # prior_summary() shows the ends of the observed range as numbers
  ps = prior_summary(fit)
  expect_equal(ps$prior[1:2], c("dunif(0.61863, 5)", "dnorm(7, 1) T(0.61863, 9.8765)"))
  expect_equal(ps$bounds[1:2], c("[0.61863, 5], ordered", "[0.61863, 9.8765], ordered"))
})


test_that("partial user-specified Dirichlet priors propagate source and description", {
  data_dir = data.frame(x = 1:10, y = 1:10)
  fit_dir = mcp(list(y ~ 1, ~ 1, ~ 1), data = data_dir, par_x = "x", prior = list(cp_1 = "dirichlet(2)"), sample = FALSE)
  verbose_dir = prior_summary(fit_dir, verbose = TRUE)

  expect_equal(verbose_dir$prior[verbose_dir$parameter == "cp_1"], "dirichlet(2)")
  expect_equal(verbose_dir$source[verbose_dir$parameter == "cp_1"], "user")
  expect_equal(verbose_dir$description[verbose_dir$parameter == "cp_1"], "User-specified prior, restricted to ordered change points")

  expect_equal(verbose_dir$prior[verbose_dir$parameter == "cp_2"], "dirichlet(2)")
  expect_equal(verbose_dir$source[verbose_dir$parameter == "cp_2"], "user")
  expect_equal(verbose_dir$description[verbose_dir$parameter == "cp_2"], "Inherited from user-specified `cp_1`, restricted to ordered change points")
})


test_that("Gaussian defaults use coherent response and link scales", {
  log_data = data.frame(
    x = 1:6,
    y = c(-10, 80, 90, 100, 110, 120)
  )
  log_fit = mcp(
    list(y ~ 1 + x), log_data,
    family = gaussian(link = "log"), sample = FALSE
  )

  # Non-positive responses are valid: the log link applies to mu, not y.
  expect_equal(log_fit$prior$Intercept_1, "dnorm(4.6 - x_1 * 1, 2.5)")
  expect_prior(log_fit$prior$x_1, paste0("dnorm(0, ", format_prior_number(2.5 / stats::sd(log_data$x)), ")"))
  expect_equal(log_fit$prior$sigma_1, "dt(0, 22.2, 3) T(0, )")
  expect_false(grepl("log\\(y\\)", log_fit$jags_code))

  rules = prior_summary(log_fit, verbose = TRUE)
  sigma_rule = rules$rule[rules$parameter == "sigma_1"]
  expect_match(sigma_rule, "mad(y)", fixed = TRUE)
  expect_false(grepl("mad(log(y))", sigma_rule, fixed = TRUE))

  wide_data = data.frame(x = 1:6, y = exp(c(-10, -5, 0, 5, 10, 15)))
  wide_fit = mcp(
    list(y ~ 1 + x), wide_data,
    family = gaussian(link = "log"), sample = FALSE
  )
  expect_equal(wide_fit$prior$Intercept_1, "dnorm(2.5 - x_1 * 1, 11.1)")  # mad(log(y)), with no values floored
  expect_prior(wide_fit$prior$x_1, paste0("dnorm(0, ", format_prior_number(2.5 / stats::sd(wide_data$x)), ")"))
  expect_equal(wide_fit$prior$sigma_1, "dt(0, 110.8, 3) T(0, )")

  zero_data = data.frame(x = 1:6, y = 0:5)
  zero_fit = mcp(
    list(y ~ 1 + x), zero_data,
    family = gaussian(link = "log"), sample = FALSE
  )
  expect_equal(zero_fit$prior$Intercept_1, "dnorm(0.9 - x_1 * 1, 2.5)")

  # Positive responses below 0.1 are not floored, so the location is scale-free (as in brms)
  tiny_data = data.frame(x = 1:6, y = c(1, 2, 3, 4, 5, 6) / 1000)
  tiny_fit = mcp(
    list(y ~ 1 + x), tiny_data,
    family = gaussian(link = "log"), sample = FALSE
  )
  tiny_prior = prior_summary(tiny_fit)$prior[prior_summary(tiny_fit)$parameter == "Intercept_1"]
  expect_equal(tiny_prior, paste0("dnorm(", round(median(log(tiny_data$y)), 1), " - x_1 * 1, 2.5)"))

  small_data = data.frame(x = 1:4, y = c(0.01, 0.02, 0.03, 0.04))
  small_fit = mcp(list(y ~ 1 + x), small_data, sample = FALSE)
  expect_prior(small_fit$prior$Intercept_1, paste0("dnorm(", format_prior_number(mean(small_data$y)), " - x_1 * 1, ", format_prior_number(2.5 * stats::sd(small_data$y)), ")"))
  expect_equal(small_fit$prior$sigma_1, "dt(0, 2.5, 3) T(0, )")
})


test_that("x-term powers are autoscaled by the sd of the power measured from min(x)", {
  data = data.frame(x = 1:10, y = 1:10)
  fit = mcp(list(y ~ 1 + x + I(x^2) + I(x^3), ~ 1), data, sample = FALSE)
  base = 2.5 * stats::sd(data$y)
  expect_prior(fit$prior$x_1, paste0("dnorm(0, ", format_prior_number(base / stats::sd(data$x)), ")"))
  expect_prior(fit$prior$xE2_1, paste0("dnorm(0, ", format_prior_number(base / stats::sd((data$x - 1)^2)), ")"))
  expect_prior(fit$prior$xE3_1, paste0("dnorm(0, ", format_prior_number(base / stats::sd((data$x - 1)^3)), ")"))
})


test_that("x-interactions are autoscaled by two-valued numeric variables but not by indicators", {
  data = data.frame(x = 1:10, y = 1:10, z = rep(c(100, 200), 5), b = rep(0:1, 5))
  fit = mcp(list(y ~ 1 + x:z + x:b, ~ 1), data, sample = FALSE)
  prior_sd = function(prior) as.numeric(sub("dnorm\\(0, (.*)\\)", "\\1", prior))
  base = 2.5 * stats::sd(data$y) / stats::sd(data$x)
  expect_equal(prior_sd(fit$prior$x.z_1), base / stats::sd(data$z))
  expect_equal(prior_sd(fit$prior$x.b_1), base)
})


test_that("default segment-1 intercept priors are placed at min(x)", {
  data = data.frame(x = 101:130, y = 1:30 + 0.5 * (1:30 > 15))
  fit = withr::with_options(list(mcp_reparam = FALSE), mcp(list(y ~ 1 + x + I(x^2), ~ 0 + x), data, par_x = "x", sample = FALSE))
  location = paste0(format_prior_number(mean(data$y)), " - x_1 * 101 - xE2_1 * 10201")
  expect_prior(fit$prior$Intercept_1, paste0("dnorm(", location, ", ", format_prior_number(2.5 * stats::sd(data$y)), ")"))
  expect_match(fit$jags_code, paste0("Intercept_1 ~ dnorm(", location, ", "), fixed = TRUE)
  ps = prior_summary(fit, verbose = TRUE)
  expect_equal(ps$rule[ps$parameter == "Intercept_1"], "dnorm(mean(y) - x_1 * min(x) - xE2_1 * min(x)^2, 2.5 * sd(y))")
  expect_match(ps$description[ps$parameter == "Intercept_1"], "on the value at min(x)", fixed = TRUE)

  # Bounded priors shift their bounds too
  data$y = rpois(30, 10)
  poisson_fit = mcp(list(y ~ 1 + x), data, family = poisson(link = "identity"), sample = FALSE)
  expect_match(poisson_fit$prior$Intercept_1, " - x_1 * 101, ", fixed = TRUE)
  expect_match(poisson_fit$prior$Intercept_1, ") T(-x_1 * 101, )", fixed = TRUE)
  expect_match(poisson_fit$jags_code, "Intercept_1_mid_ ~ dnorm(", fixed = TRUE)
  data$y = rbinom(30, 1, 0.5)
  bernoulli_fit = mcp(list(y ~ 1 + x), data, family = bernoulli(link = "identity"), sample = FALSE)
  expect_equal(bernoulli_fit$prior$Intercept_1, "dunif(-x_1 * 101, 1 - x_1 * 101)")
  expect_match(bernoulli_fit$jags_code, "Intercept_1_mid_ ~ dunif(", fixed = TRUE)

  # Factor levels and contrasts with their own x-terms are placed at min(x) too
  data$state = rep(c("a", "b"), 15)
  level_fit = mcp(list(y ~ 0 + state + state:x), data, sample = FALSE)
  expect_match(level_fit$prior$stateb_1, " - stateb.x_1 * 101, ", fixed = TRUE)
  expect_match(level_fit$jags_code, "stateb_1_mid_ ~ dnorm(", fixed = TRUE)
  contrast_fit = mcp(list(y ~ 1 + state * x), data, sample = FALSE)
  expect_match(contrast_fit$prior$stateb_1, "dnorm(-stateb.x_1 * 101, ", fixed = TRUE)
  expect_match(contrast_fit$jags_code, "stateb_1_mid_ ~ dnorm(", fixed = TRUE)
  ps = prior_summary(contrast_fit, verbose = TRUE)
  expect_match(ps$description[ps$parameter == "stateb_1"], "on the difference at min(x)", fixed = TRUE)

  # User-specified intercept priors stay on the intercept at x = 0
  user_fit = mcp(list(y ~ 1 + x), data, prior = list(Intercept_1 = "dnorm(0, 10)"), sample = FALSE)
  expect_match(user_fit$jags_code, "Intercept_1 ~ dnorm(", fixed = TRUE)
  expect_false(grepl("Intercept_1_", user_fit$jags_code))

  # Nothing to move when x starts at zero or segment 1 has no x-terms
  expect_false(grepl("x_", mcp(list(y ~ 1 + x), transform(data, x = x - 101), sample = FALSE)$prior$Intercept_1))
  expect_false(grepl("x_", mcp(list(y ~ 1, ~ 0 + x), data, par_x = "x", sample = FALSE)$prior$Intercept_1))
})


test_that("the prior on the value at min(x) matches the default intercept prior", {
  data = data.frame(x = 101:130, y = 1:30)
  get_draws = function(prior = list()) {
    fit = suppressWarnings(mcp(
      list(y ~ 1 + x), data, prior = prior, sample = "prior",
      chains = 1, iter = 4000, warmup = 100, seed = 42, quiet = TRUE
    ))
    list(fit = fit, draws = as_draws_df(fit, prior = TRUE))
  }
  default = get_draws()
  level = default$draws$Intercept_1 + default$draws$x_1 * 101
  expect_equal(mean(level), mean(data$y), tolerance = 0.1)
  expect_equal(stats::sd(level), 2.5 * stats::sd(data$y), tolerance = 0.1)

  # Reusing fit$prior gives the same model
  copy = get_draws(as.list(default$fit$prior))
  expect_equal(stats::sd(copy$draws$Intercept_1), stats::sd(default$draws$Intercept_1), tolerance = 0.1)
  expect_equal(stats::sd(copy$draws$Intercept_1 + copy$draws$x_1 * 101), stats::sd(level), tolerance = 0.1)
})


test_that("factors without an intercept get intercept priors where a segment starts at zero", {
  data = data.frame(x = 101:200, state = factor(rep(c("A", "B"), 50)), y = rnorm(100, 50))
  intercept = mcp(list(y ~ 1, ~ 1), data, par_x = "x", sample = FALSE)$prior$Intercept_1
  get_prior = function(model, family = gaussian())
    mcp(model, data, family = family, par_x = "x", sample = FALSE)$prior

  # Segments that start at a fixed zero: each factor level gets the intercept prior
  expect_equal(get_prior(list(y ~ 0 + state, ~ 1))[c("stateA_1", "stateB_1")], list(stateA_1 = intercept, stateB_1 = intercept))
  expect_equal(get_prior(list(y ~ 0, ~ 0 + state))$stateA_2, intercept)
  data$y = rbinom(100, 1, 0.5)
  expect_equal(get_prior(list(y ~ 0 + state, ~ 1), bernoulli(link = "identity"))$stateA_1, "dunif(0, 1)")

  # Otherwise contrast and slope priors, e.g., jumps from an estimated value
  expect_match(get_prior(list(y ~ 0 + x, ~ 0 + state))$stateA_2, "^dnorm\\(0, ")
  expect_match(get_prior(list(y ~ 1, ~ 1 + state))$stateB_2, "^dnorm\\(0, ")
  expect_match(get_prior(list(y ~ 0 + state:x, ~ 1))$stateA.x_1, "^dnorm\\(0, ")

  # Each factor level's prior is on its value at min(x), as for Intercept_1
  data$y = rnorm(100, 50)
  fit = suppressWarnings(mcp(
    list(y ~ 0 + state + x), data, par_x = "x", sample = "prior",
    chains = 1, iter = 4000, warmup = 100, seed = 42, quiet = TRUE
  ))
  draws = as_draws_df(fit, prior = TRUE)
  for (name in c("stateA_1", "stateB_1")) {
    level = draws[[name]] + draws$x_1 * 101
    expect_equal(mean(level), mean(data$y), tolerance = 0.1)
    expect_equal(stats::sd(level), 2.5 * stats::sd(data$y), tolerance = 0.1)
  }
})


test_that("a constant response uses unit scale for sd(y)-based default priors", {
  data = data.frame(x = 1:10, y = 2)
  fit = mcp(list(y ~ 1 + x), data, sample = FALSE)
  expect_equal(fit$prior$Intercept_1, "dnorm(2 - x_1 * 1, 2.5)")
  expect_prior(fit$prior$x_1, paste0("dnorm(0, ", format_prior_number(2.5 / stats::sd(data$x)), ")"))
})


test_that("reparameterizations keep segment lines in place when change points move", {
  data = data.frame(x = 1:100, y = rnorm(100))
  model = list(y ~ 1 + x, ~ 1 + x, ~ 0 + x + I(x^2))
  fit = mcp(model, data, par_x = "x", sample = FALSE)
  # Segment 1: intercept sampled as the value at the expected midpoint of segment 1, 1 + 99 / 6
  expect_match(fit$jags_code, "Intercept_1 = Intercept_1_mid_ - (x_1 * (17.5))", fixed = TRUE)
  # Disjoined segment: intercept sampled as the value at the segment end
  expect_match(fit$jags_code, "Intercept_2 = Intercept_2_end_ - (x_2 * (cp_2 - cp_1))", fixed = TRUE)
  # Joined segment: first x-term sampled as the value at the segment end, others as the rise
  expect_match(fit$jags_code, "x_3 = (x_3_end_ - (Intercept_2 + x_2 * (cp_2 - cp_1))) / (cp_3 - cp_2)", fixed = TRUE)
  expect_match(fit$jags_code, "xE2_3 = xE2_3_rise_ / (cp_3 - cp_2)^2", fixed = TRUE)
  # Distributional parameters follow the same rules
  sigma_fit = mcp(list(y ~ 1, ~ 0 + sigma(1 + x)), data, par_x = "x", sample = FALSE)
  expect_match(sigma_fit$jags_code, "sigma_Intercept_2 = sigma_Intercept_2_end_ - (sigma_x_2 * (cp_2 - cp_1))", fixed = TRUE)

  # The plain model
  plain_fit = withr::with_options(list(mcp_reparam = FALSE), mcp(model, data, par_x = "x", sample = FALSE))
  expect_false(grepl("_mid_|_end_|_rise_", plain_fit$jags_code))
  expect_match(plain_fit$jags_code, "\n  x_3 ~ dnorm")
  expect_equal(plain_fit$prior, fit$prior)

  # Not for user-specified priors or segments bounded by group-level change points
  user_fit = mcp(list(y ~ 1, ~ 0 + x), data, par_x = "x", prior = list(x_2 = "dnorm(0, 1)"), sample = FALSE)
  expect_false(grepl("_end_|_rise_", user_fit$jags_code))
  data$id = rep(1:5, 20)
  group_fit = mcp(list(y ~ 1, ~ cp(1 + (1 | id)) + 0 + x), data, par_x = "x", sample = FALSE)
  expect_false(grepl("_end_|_rise_", group_fit$jags_code))
})


test_that("reparameterizations imply the same priors", {
  data = data.frame(x = 101:200, y = rnorm(100))
  model = list(y ~ 1 + x, ~ 0 + x, ~ 1 + x + sigma(1 + x))
  get_quantiles = function(reparam) {
    withr::local_options(mcp_reparam = reparam)
    fit = suppressWarnings(mcp(
      model, data, par_x = "x", sample = "prior",
      chains = 2, iter = 5000, warmup = 200, seed = 42, quiet = TRUE
    ))
    draws = as_draws_df(fit, prior = TRUE)
    sapply(mcp_pars(fit)$parameter, function(name) stats::quantile(draws[[name]], c(0.1, 0.5, 0.9)))
  }
  q_reparam = get_quantiles(TRUE)
  q_plain = get_quantiles(FALSE)
  spread = q_plain[3, ] - q_plain[1, ]
  for (name in colnames(q_plain))
    expect_true(all(abs(q_reparam[, name] - q_plain[, name]) < 0.1 * spread[name]), label = name)
})
