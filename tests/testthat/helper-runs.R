get_mcp_test_level = function() {
  if (identical(Sys.getenv("MCP_TEST_LEVEL"), "release")) return("release")
  "default"
}

# testthat deletes any vdiffr/snapshot file it didn't see referenced during the
# current run (testthat:::SnapshotReporter$end_reporter), *unless* it thinks
# it's running on CI (Sys.getenv("CI") == "true"). Since test-fits-example-*.R
# files (which own tests/testthat/_snaps/fits-examples/) only run
# at MCP_TEST_LEVEL = "release", every default-level run would otherwise wipe
# those reference snapshots from disk. Spoof CI so cleanup is skipped whenever
# we're not doing a release-level run; leave it alone at release level so
# snapshot review/cleanup keeps working normally.
if (get_mcp_test_level() != "release" && !isTRUE(as.logical(Sys.getenv("CI", "false")))) {
  Sys.setenv(CI = "true")
}


test_runs = function(model,
                     data = data_gauss,
                     prior = list(),
                     family = gaussian(),
                     par_x = "x",
                     sample = TRUE,
                     test_s3 = FALSE) {

  # Without sampling, on a data.frame.
  empty = mcp(
    model = model,
    data = data,
    prior = prior,
    family = family,
    par_x = par_x,
    sample = FALSE,
    quiet = TRUE
  )

  testthat::expect_true(is.list(empty$model), model)
  testthat::expect_true(is.data.frame(empty$data), model)
  testthat::expect_true(is.list(empty$prior), model)
  testthat::expect_true(all(class(empty$family) == c("mcpfamily", "family")), model)
  testthat::expect_true(is.null(empty$samples), model)
  testthat::expect_true(is.null(empty$loglik), model)
  testthat::expect_true(is.null(empty$loo), model)
  testthat::expect_true(is.null(empty$waic), model)
  testthat::expect_true(is.list(empty$pars), model)
  testthat::expect_true(is.character(empty$pars$population), model)
  testthat::expect_true((is.character(empty$pars$group) | is.null(empty$pars$group)), model)
  testthat::expect_true(is.character(empty$pars$x), model)
  testthat::expect_true(is.character(empty$pars$y), model)
  testthat::expect_true(is.character(empty$jags_code), model)
  testthat::expect_true(is.function(empty$simulate), model)
  testthat::expect_true(is.list(empty$.internal), model)

  level = get_mcp_test_level()

  if (sample) {
    if (rbinom(1, 1, 0.5) == 1)
      data = tibble::as_tibble(data)

    warnings = messages = character()
    fit = withCallingHandlers(
      mcp(
        model = model,
        data = data,
        family = family,
        sample = "both",  # prior and posterior to check hypotheses
        par_x = par_x,
        adapt = 6,
        iter = 18,  # Keep enough draws for the minimum supported loo version.
        chains = 2,  # run sequentially under the default future plan
        quiet = TRUE
      ),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      },
      message = function(condition) {
        messages <<- c(messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      }
    )

    # Allow for known messages and wornings that does not signify errors
    accepted_warnings = c(
      "Adaptation incomplete",  # due to very small test datasets
      "Some parameters may not have converged well",  # ditto - too few iterations to mix
      "Posterior AR/MA root smoke test"
    )
    accepted_messages = c(
      "Autoregression currently assumes homoskedasticity",
      "You are using ar\\(\\) together"
    )

    for (warn in warnings) {
      if (!any(stringr::str_starts(warn, accepted_warnings))) {
        testthat::fail("Got an unknown warning: ", warn)
      }
    }
    for (msg in messages) {
      if (!any(stringr::str_starts(msg, accepted_messages))) {
        testthat::fail("Got an unknown message: ", msg)
      }
    }

    # Assign globally so errors can be inspected upon hard fail

    if (is_arma(fit))
      test_arma_simulation(fit)

    if (test_s3 || level == "release") {
      test_s3_methods(fit)
    }
  }
}


test_s3_methods = function(fit) {
  # Test criterions. Will warn about very few samples
  if (!is.null(.subset2(fit, "mcmc_post"))) {
    fit$loo = suppressMessages(suppressWarnings(loo(fit)))
    fit$waic = suppressMessages(suppressWarnings(waic(fit)))
    testthat::expect_true(loo::is.psis_loo(fit$loo))
    testthat::expect_true(loo::is.waic(fit$waic))

    # Test pointwise
    fit$loo_pointwise = suppressMessages(suppressWarnings(loo(fit, pointwise = TRUE)))
    rownames(fit$loo$pointwise) = NULL
    testthat::expect_equal(fit$loo$pointwise, fit$loo_pointwise$pointwise)
  }

  for (col in c("mcmc_post", "mcmc_prior")) {
    use_prior = col == "mcmc_prior"
    fit_to_test = fit

    # Test the informative fallback once. All intentional prior tests below
    # request prior samples explicitly to avoid repeated messages.
    if (use_prior) {
      fit_to_test$mcmc_post = NULL
      testthat::expect_message(
        mcmclist_samples(fit_to_test),
        "Posterior was not sampled. Using prior samples"
      )
    }

    # Check that samples are the correct format
    samples_col = .subset2(fit_to_test, col)
    testthat::expect_true(is.list(samples_col))
    testthat::expect_true(coda::is.mcmc(samples_col[[1]]))
    testthat::expect_true(all(fit_to_test$pars$population %in% colnames(samples_col[[1]])))

    # Test mcpfit functions
    varying_cols = na.omit(get_fit_model_tables(fit_to_test)$group_effects$group_col)
    test_summary(fit_to_test, varying_cols, prior = use_prior)
    test_plot_pars(fit_to_test, prior = use_prior)
    test_pp_eval(fit_to_test, prior = use_prior)
    test_hypothesis(fit, prior = use_prior)
  }
}


# Exercise fresh-series simulation for every accepted AR/MA model. Detailed
# recursion and R/JAGS agreement are tested separately in test-validate-garma.R.
test_arma_simulation = function(fit) {
  control_args = c("fit", "newdata", ".type", ".rate", ".dpar", ".arma", ".scale")
  component_args = setdiff(names(formals(fit$simulate)), control_args)
  component_values = stats::setNames(as.list(rep(0, length(component_args))), component_args)
  call = c(
    list(fit = fit, newdata = fit$data),
    component_values,
    list(.type = "predict")
  )
  simulated = suppressWarnings(suppressMessages(do.call(fit$simulate, call)))

  testthat::expect_length(simulated, nrow(fit$data))
  testthat::expect_true(is.numeric(simulated))
  testthat::expect_false(anyNA(simulated))
  testthat::expect_true(all(is.finite(simulated)))

  if (fit$family$family == "binomial") {
    trials = fit$data[[fit$pars$trials]]
    testthat::expect_true(all(simulated >= 0 & simulated <= trials))
    testthat::expect_true(all(simulated == floor(simulated)))
  } else if (fit$family$family %in% c("poisson", "negbinomial")) {
    testthat::expect_true(all(simulated >= 0))
    testthat::expect_true(all(simulated == floor(simulated)))
  }
}


# Tests if summary(fit) and ranef(fit) work as expected
test_summary = function(fit, varying_cols, prior = FALSE) {
  summary_cols = c('name','mean','lower','upper','rhat','ess_bulk','ess_tail')
  verbose_summary_cols = c('name','segment','dpar','mean','lower','upper','rhat','ess_bulk','ess_tail')
  if (!is.null(attr(fit$data[, fit$pars$y], "simulated"))) {
    summary_cols = append(summary_cols, c("match", "sim"), after = 1)
    verbose_summary_cols = append(verbose_summary_cols, c("match", "sim"), after = 3)
  }
  output = capture.output(result <- summary(fit, prior = prior))
  testthat::expect_named(result, summary_cols)
  testthat::expect_true(all(result$name %in% fit$pars$population))  # All parameters
  capture.output(verbose_result <- summary(fit, prior = prior, verbose = TRUE))
  testthat::expect_named(verbose_result, verbose_summary_cols)
  testthat::expect_named(fixef(fit, prior = prior), summary_cols)
  testthat::expect_named(fixef(fit, prior = prior, verbose = TRUE), verbose_summary_cols)

  # If there are group-level effects
  if (length(varying_cols) > 0) {
    testthat::expect_true(any(grepl("ranef\\(", output)))  # notice about group-level effects
    varying = ranef(fit, prior = prior)
    testthat::expect_true(is.character(varying$name))
    testthat::expect_true(is.numeric(varying$mean))
    testthat::expect_named(ranef(fit, prior = prior, verbose = TRUE), verbose_summary_cols)

    effects = get_fit_model_tables(fit)$group_effects
    group_level_counts = lapply(
      effects$group_col,
      function(col) length(unique(fit$data[[col]]))
    )
    n_varying_levels = sum(unlist(group_level_counts))
    testthat::expect_equal(nrow(varying), n_varying_levels)
  }
}


# Test the regular plot, including faceting
test_plot = function(fit, varying_cols) {
  q_fit = rbinom(1, 1, 0.5) == 1  # add quantiles sometimes
  q_predict = rbinom(1, 1, 0.5) == 1  # add quantiles sometimes
  # To facet or not to facet
  if (length(varying_cols) > 0) {
    gg = try(plot(fit, facet_by = varying_cols[1], color_by = NULL, q_fit = q_fit, q_predict = q_predict, lines = 3, ndraws = NULL), silent = TRUE)  # just take the first
  } else {
    gg = try(plot(fit, color_by = NULL, q_fit = q_fit, q_predict = q_predict, lines = 3, ndraws = NULL), silent = TRUE)
  }
  # Is it a ggplot or a known error?
  if (inherits(gg, "try-error")) {
    # (the error is an artifact of very small test data --> wide posteriors.)
    if (fit$family$family == "poisson") {
      expected_error = "predict"  # column "predict" OK: a side-effect of the small data and short sampling.
    } else {
      expected_error = ">>>>do_not_expect_any_errors<<<<<"
    }
    error_message = attr(gg, "condition")$message
    is_expected = any(stringr::str_detect(error_message, expected_error))
    testthat::expect_true(is_expected)
  } else {
    testthat::expect_s3_class(gg, "ggplot")
  }
}

# Test plot() calls to bayesplot
test_plot_pars = function(fit, prior = FALSE) {
  gg = plot_pars(
    fit,
    type = "dens_overlay",
    prior = prior,
    nvariables = length(fit$pars$population)
  )
  testthat::expect_s3_class(gg, "ggplot")
}



test_hypothesis = function(fit, prior) {
  # Function to test directional, interval, and point equality hypotheses
  run_test_hypothesis = function(fit, base, prior = prior) {
    hypotheses = paste0(base, " > 1")  # Directional
    if (prior == FALSE) {
      hypotheses = c(
        hypotheses,
        paste0(base, " = -1"),  # Savage-Dickey (point); only works if both prior and posterior is present.
        paste0(base, " > -10 & ", base, " < 10")  # Interval hypothesis requiring prior
      )
    }

    result = suppressWarnings(hypothesis(fit, hypotheses, prior = prior))
    testthat::expect_true(is.data.frame(result) & nrow(result) == length(hypotheses))
  }

  # Test single pop effect
  run_test_hypothesis(fit, fit$pars$population[1], prior = prior)

  # Test multiple pop effect
  if (length(fit$pars$population) > 1)
    run_test_hypothesis(fit, paste0(fit$pars$population[1] , " + ", fit$pars$population[2]), prior = prior)

  # Varying
  if (!is.null(fit$pars$group)) {
    mcmc_vars = colnames(mcmclist_samples(fit)[[1]])
    varying_starts = paste0("^", fit$pars$group[1], "\\[")
    varying_col_ids = stringr::str_detect(mcmc_vars, varying_starts)
    varying_cols = paste0("`", mcmc_vars[varying_col_ids], "`")  # Add these for varying

    # Test one group-level deviation
    run_test_hypothesis(fit, varying_cols[1], prior = prior)

    # Test multiple group-level deviations
    if (length(varying_cols) > 1)
      run_test_hypothesis(fit, paste0(varying_cols[1], " + ", varying_cols[2]), prior = prior)
  }
}


test_pp_eval_func = function(fit, func, colname, prior = FALSE) {
  # Settings
  varying_cols = na.omit(unique(
    get_fit_model_tables(fit)$group_effects$group_col
  ))
  rhs_cols = intersect(
    setdiff(get_rhs_vars(fit$model), varying_cols),
    colnames(fit$data)
  )
  expected_colnames = c(
    fit$pars$x,
    fit$pars$trials,
    rhs_cols,
    varying_cols,
    colname, "error", "Q2.5", "Q97.5"  # substitute-stuff just gets the func name as string
  )
  if (length(fit$pars$arma) > 0 || colname %in% c("loglik", "residuals"))
    expected_colnames = c(expected_colnames, fit$pars$y)

  # `log_lik()` follows the conventional draws-by-observation matrix default.
  # Its former observation-level data-frame output remains available explicitly.
  if (colname == "loglik") {
    default_result = try(func(fit, prior = prior), silent = TRUE)
    if (!inherits(default_result, "try-error")) {
      testthat::expect_true(is.matrix(default_result))
      testthat::expect_equal(
        ncol(default_result),
        sum(!is.na(fit$data[[fit$pars$y]]))
      )
      testthat::expect_equal(
        default_result,
        func(fit, prior = prior, summary = FALSE, samples_format = "matrix")
      )
    }
  }

  # Run and test
  result = if (colname == "loglik") {
    try(func(fit, prior = prior, summary = TRUE), silent = TRUE)
  } else {
    try(func(fit, prior = prior), silent = TRUE)
  }
  if (inherits(result, "try-error")) {
    error_message = as.character(result)
    is_count_family = fit$family$family %in% c("poisson", "negbinomial")
    expected_error_count = "Modelled extremely large count mean"  # A test-specific side-effect of the small data and short sampling.
    testthat::expect_true(!is_count_family | stringr::str_detect(error_message, stringr::fixed(expected_error_count)))
  } else {
    testthat::expect_true(is.data.frame(result))
    testthat::expect_equal(nrow(result), nrow(fit$data))  # Returns same number of rows as data
    if (fit$family$family == "poisson") {
      # Extremely diffuse draws from these deliberately tiny test fits can
      # overflow derived Poisson summaries. Data keys must remain intact and
      # the estimate itself must not be entirely missing.
      data_cols = intersect(
        c(fit$pars$x, fit$pars$y, fit$pars$trials, rhs_cols, varying_cols),
        colnames(result)
      )
      testthat::expect_false(anyNA(result[, data_cols, drop = FALSE]))
      testthat::expect_false(all(is.na(result[, colname])))
    } else {
      testthat::expect_false(anyNA(result))
    }
    testthat::expect_true(dplyr::setequal(colnames(result), expected_colnames))  # Exactly these columns regardless of order
    testthat::expect_true(all(result[, fit$pars$x] == fit$data[, fit$pars$x]))  # Output should have same order as input

    if (colname == "residuals") {
      fitted_result = fitted(fit, rate = FALSE, prior = prior)  # residuals are on observed-data scale
      observed = fit$data[[fit$pars$y]]
      testthat::expect_equal(result$residuals, observed - fitted_result$fitted)
      testthat::expect_equal(result$error, fitted_result$error)
      testthat::expect_equal(result$Q2.5, observed - fitted_result$Q97.5)
      testthat::expect_equal(result$Q97.5, observed - fitted_result$Q2.5)
    }
  }
}


# Weighted Gaussian evaluation must use the same observation-level SD as JAGS:
# precision = weight / sigma^2, or equivalently SD = sigma / sqrt(weight).
test_pp_eval_weights = function(fit, prior = FALSE) {
  if (fit$family$family != "gaussian" || length(fit$pars$weights) == 0)
    return(invisible(NULL))

  weight_col = fit$pars$weights
  keys = c(".chain", ".iteration", ".draw", "data_row")
  mu = fitted(fit, summary = FALSE, probs = FALSE, prior = prior, dpar = "mu")
  sigma = fitted(fit, summary = FALSE, probs = FALSE, prior = prior, dpar = "sigma")
  loglik = log_lik(fit, summary = FALSE, probs = FALSE, prior = prior)
  weights = fit$data[[weight_col]][loglik$data_row]
  observation_sd = sigma$.epred / sqrt(weights)
  observed = fit$data[[fit$pars$y]][loglik$data_row]

  testthat::expect_equal(loglik[, keys], mu[, keys])
  testthat::expect_equal(loglik[, keys], sigma[, keys])
  testthat::expect_equal(
    loglik$.loglik,
    stats::dnorm(observed, mu$.epred, observation_sd, log = TRUE)
  )

  had_random_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed)
    random_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(123)
  prediction = predict(fit, summary = FALSE, probs = FALSE, prior = prior)
  set.seed(123)
  expected = stats::rnorm(nrow(prediction), mu$.epred, observation_sd)
  testthat::expect_equal(prediction[, keys], mu[, keys])
  testthat::expect_equal(prediction$.prediction, expected)

  if (!prior) {
    newdata_without_weights = fit$data[, colnames(fit$data) != weight_col, drop = FALSE]
    testthat::expect_error(
      predict(fit, newdata = newdata_without_weights, summary = FALSE),
      weight_col,
      fixed = TRUE
    )
    testthat::expect_error(
      log_lik(fit, newdata = newdata_without_weights, summary = FALSE),
      weight_col,
      fixed = TRUE
    )
  }
}


test_pp_eval = function(fit, prior = FALSE) {
  # Test pp_eval
  test_pp_eval_func(fit, fitted, "fitted", prior = prior)
  test_pp_eval_func(fit, predict, "predict", prior = prior)
  test_pp_eval_func(fit, residuals, "residuals", prior = prior)
  test_pp_eval_func(fit, log_lik, "loglik", prior = prior)
  test_pp_eval_weights(fit, prior = prior)

  # Test the other arguments. Inside "try" without further checking because such errors should be caught by the above.
  result_more = try(fitted(
    fit,
    newdata = fit$data[sample(nrow(fit$data), 3), ],
    summary = FALSE,
    probs = c(0.1, 0.5, 0.999),
    prior = prior,
    varying = "cp",
    ndraws = 2,
    arma = FALSE
  ), silent = TRUE)

  if (is.data.frame(result_more)) {
    #testthat::expect_true(nrow(result_more) == nrows * ndraws * 2)  # nrows * ndraws * nchains
    testthat::expect_true(sum(is.na(result_more)) == 0)

    group_effects = get_fit_model_tables(fit)$group_effects
    all_group_cols = unique(stats::na.omit(group_effects$group_col))
    selected_cp = unpack_varying(fit, pars = "cp")
    rhs_cols = intersect(
      setdiff(get_rhs_vars(fit$model), all_group_cols),
      colnames(fit$data)
    )
    expected_colnames_more = c(
      # Tidybayes stuff
      ".chain", ".iteration", ".draw",

      # Model parameters
      fit$pars$population,
      fit$pars$group,

      # Predictors
      fit$pars$trials,
      fit$pars$x,
      rhs_cols,
      selected_cp$cols,
      "data_row",
      ".epred"  # dot-prefixed for summary = FALSE
    )

    testthat::expect_true(dplyr::setequal(colnames(result_more), expected_colnames_more))  # Exactly these columns regardless of order
  }

  # Population-only evaluation should not require any grouping columns.
  if (length(fit$pars$group) > 0) {
    varying_cols = stats::na.omit(unique(get_fit_model_tables(fit)$group_effects$group_col))
    population_newdata = fit$data[, colnames(fit$data) %notin%
      c(fit$pars$y, varying_cols), drop = FALSE]
    population_result = fitted(
      fit,
      newdata = population_newdata,
      summary = FALSE,
      probs = FALSE,
      prior = prior,
      varying = FALSE,
      arma = FALSE,
      ndraws = 2
    )

    testthat::expect_equal(
      sort(unique(population_result$data_row)),
      seq_len(nrow(population_newdata))
    )
    testthat::expect_false(any(varying_cols %in% colnames(population_result)))
  }

  # Test pp_check
  if (length(fit$pars$group) > 0) {
    varying_col = na.omit(get_fit_model_tables(fit)$group_effects$group_col)[1]  # Just use the first column
    pp_default = try(suppressWarnings(pp_check(fit, facet_by = varying_col, ndraws = 2, prior = prior)), silent = TRUE)
  } else {
    pp_default = try(suppressWarnings(pp_check(fit, ndraws = 2, prior = prior)), silent = TRUE)
  }

  if (inherits(pp_default, "try-error")) {
    error_message = as.character(pp_default)

    # Expected count-family errors from the intentionally tiny data, extremely
    # short sampling, and broad prior predictive distributions.
    if (fit$family$family %in% c("poisson", "negbinomial")) {
      testthat::expect_true(stringr::str_detect(error_message, "Modelled extremely large count mean"))
    } else {
      # Fail otherwise
      testthat::expect_true(error_message)
    }
  } else {
    testthat::expect_s3_class(pp_default, "ggplot")
  }
}



# Rutine for testing a list of erroneous models
test_bad = function(models, ...) {
  for (model in models) {
    test_name = paste0(as.character(substitute(models)), ":
    ", paste0(model, collapse=", "))

    testthat::test_that(test_name, {
      suppressWarnings(testthat::expect_error(test_runs(model, sample = FALSE, ...)))  # should err before sampling
    })
  }
}


# Routine for testing a list of good models
test_good = function(models, ...) {
  for (model in models) {
    test_name = paste0(as.character(substitute(models)), ":
    ", paste0(model, collapse=", "))

    testthat::test_that(test_name, {
      test_runs(model, ...)
    })
  }
}
