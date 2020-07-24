# Use all mcp functions on a given model to check that it does not result in
# errors.
test_runs = function(model,
                     data = data_gauss,
                     prior = list(),
                     family = gaussian(),
                     par_x = "x",
                     sample = TRUE) {

  # Without sampling, on a data.frame.
  empty = mcp(
    model = model,
    data = data,
    prior = prior,
    family = family,
    par_x = par_x,
    sample = FALSE
  )

  # With (very brief!) sampling, on a tibble
  # Just to leverage JAGS code checking and the mcpfit data structure
  if (sample == TRUE) {
    # If sample = FALSE, it should pass/fail with the above. If TRUE,
    # check for correct types in data structure
    testthat::expect_true(is.list(empty$model), model)
    testthat::expect_true(is.data.frame(empty$data), model)
    testthat::expect_true(is.list(empty$prior), model)
    testthat::expect_true(class(empty$family) == "family", model)
    testthat::expect_true(is.null(empty$samples), model)
    testthat::expect_true(is.null(empty$loglik), model)
    testthat::expect_true(is.null(empty$loo), model)
    testthat::expect_true(is.null(empty$waic), model)
    testthat::expect_true(is.list(empty$pars), model)
    testthat::expect_true(is.character(empty$pars$population), model)
    testthat::expect_true((is.character(empty$pars$varying) | is.null(empty$pars$varying)), model)
    testthat::expect_true(is.character(empty$pars$x), model)
    testthat::expect_true(is.character(empty$pars$y), model)
    testthat::expect_true(is.character(empty$jags_code), model)
    testthat::expect_true(is.function(empty$simulate), model)
    testthat::expect_true(is.list(empty$.other), model)

    # Should work for tibbles as well. So do this sometimes
    if (rbinom(1, 1, 0.5) == 1)
      data = tibble::as_tibble(data)

    # Initiates and shut downs sampling much faster with no multicore
    options(mc.cores = 1)

    # Capture (expected) messages and warnings
    quiet_out = purrr::quietly(mcp)(  # Do not print to console
      model = model,
      data = data,
      family = family,
      sample = "both",  # prior and posterior to check hypotheses
      par_x = par_x,
      adapt = 6,
      iter = 18,  # loo fails if this is too low. TO DO: require next version of loo when it is out.
      chains = 2,  # 1 or 2
      cores = 1  # run serial for faster init. Parallel can be trused to just work.
    )

    # Allow for known messages and wornings that does not signify errors
    accepted_warnings = c("Adaptation incomplete")  # due to very small test datasets
    accepted_messages = c(
      "Finished sampling in",
      "The current implementation of autoregression can be fragile",
      "Autoregression currently assumes homoskedasticity",
      "You are using ar\\(\\) together"
    )

    for (warn in quiet_out$warnings) {
      if (!any(stringr::str_starts(warn, accepted_warnings))) {
        testthat::fail("Got an unknown warning: ", warn)
      }
    }
    for (msg in quiet_out$messages) {
      if (!any(stringr::str_starts(msg, accepted_messages))) {
        testthat::fail("Got an unknown message: ", msg)
      }
    }

    # Assign globally so errors can be inspected upon hard fail
    fit = quiet_out$result


    # Test criterions. Will warn about very few samples
    if (!is.null(fit$mcmc_post)) {
      fit$loo = suppressWarnings(loo(fit))
      fit$waic = suppressWarnings(waic(fit))
      testthat::expect_true(loo::is.psis_loo(fit$loo))
      testthat::expect_true(loo::is.waic(fit$waic))
    }

    # Test hypothesis
    test_hypothesis(fit)

    for (col in c("mcmc_post", "mcmc_prior")) {
      # To test the prior, try setting mcmc_post = NULL to force use of prior
      # (mcmclist_samples checks for NULL)
      if (col == "mcmc_prior")
        fit$mcmc_post = NULL

      # Check that samples are the correct format
      testthat::expect_true(is.list(fit[[col]]), model)
      testthat::expect_true(coda::is.mcmc(fit[[col]][[1]]), model)
      testthat::expect_true(all(fit$pars$population %in% colnames(fit[[col]][[1]])))

      # Test mcpfit functions
      varying_cols = na.omit(fit$.other$ST$cp_group_col)
      test_summary(fit, varying_cols)
      test_plot(fit, varying_cols)  # default plot
      test_plot_pars(fit)  # bayesplot call
      test_pp_eval(fit)
    }
  }
}


# Tests if summary(fit) and ranef(fit) work as expected
test_summary = function(fit, varying_cols) {
  summary_cols = c('name','mean','lower','upper','Rhat','n.eff')
  result = purrr::quietly(summary)(fit)$result  # Do not print to console
  output = purrr::quietly(summary)(fit)$output  # Do not print to console
  testthat::expect_true(all(colnames(result) %in% summary_cols))  # All columns
  testthat::expect_true(all(result$name %in% fit$pars$population))  # All parameters

  # If there are varying effects
  if (length(varying_cols) > 0) {
    testthat::expect_match(output, "ranef\\(")  # noticed about varying effects
    varying = ranef(fit)
    testthat::expect_true(is.character(varying$name))
    testthat::expect_true(is.numeric(varying$mean))

    group_level_counts = lapply(varying_cols, function(col) length(fit$data[, col]))
    n_unique_data = sum(unlist(group_level_counts))
    testthat::expect_true(nrow(varying) == n_unique_data)  # TO DO: should fail if there are multiple groups
  }
}

# Test the regular plot, including faceting
test_plot = function(fit, varying_cols) {
  q_fit = rbinom(1, 1, 0.5) == 1  # add quantiles sometimes
  q_predict = rbinom(1, 1, 0.5) == 1  # add quantiles sometimes
  # To facet or not to facet
  if (length(varying_cols) > 0) {
    gg = try(plot(fit, facet_by = varying_cols[1], q_fit = q_fit, q_predict = q_predict, lines = 3, nsamples = NULL), silent = TRUE)  # just take the first
  } else {
    gg = try(plot(fit, q_fit = q_fit, q_predict = q_predict, lines = 3, nsamples = NULL), silent = TRUE)
  }
  # Is it a ggplot or a known error?
  if (inherits(gg, "try-error")) {
    # (the error is an artefact of very small test data --> wide posteriors.)
    if (fit$family$family == "poisson") {

      expected_error = "Problem with \`mutate\\(\\)\` input \`predict\`"
    } else if (any(stringr::str_detect(fit$pars$sigma, "^sigma_.*_.*$"))) {  # for slopes on sigma
      expected_error = "Modelled negative sigma"
    } else {
      expected_error = ">>>>do_not_expect_any_errors<<<<<"
    }
    error_message = attr(gg, "condition")$message
    is_expected = any(stringr::str_starts(error_message, expected_error))
    expect_true(is_expected)
  } else {
    testthat::expect_true(ggplot2::is.ggplot(gg))
  }
}

# Test plot() calls to bayesplot
test_plot_pars = function(fit) {
  gg = plot_pars(fit, type = "dens_overlay")
  testthat::expect_true(ggplot2::is.ggplot(gg))
}



test_hypothesis = function(fit) {
  # Function to test both directional and point hypotheses
  run_test_hypothesis = function(fit, base) {
    hypotheses = c(
      paste0(base, " > 1"),  # Directional
      paste0(base, " = -1")  # Savage-Dickey (point)
    )
    result = hypothesis(fit, hypotheses)
    testthat::expect_true(is.data.frame(result) & nrow(result) == 2)
  }

  # Test single pop effect
  run_test_hypothesis(fit, fit$pars$population[1])

  # Test multiple pop effect
  if (length(fit$pars$population) > 1)
    run_test_hypothesis(fit, paste0(fit$pars$population[1] , " + ", fit$pars$population[2]))

  # Varying
  if (!is.null(fit$pars$varying)) {
    mcmc_vars = colnames(mcmclist_samples(fit)[[1]])
    varying_starts = paste0("^", fit$pars$varying[1])
    varying_col_ids = stringr::str_detect(mcmc_vars, varying_starts)
    varying_cols = paste0("`", mcmc_vars[varying_col_ids], "`")  # Add these for varying

    # Test single varying effect
    run_test_hypothesis(fit, varying_cols[1])

    # Test multiple varying effects
    if (length(varying_cols) > 1)
      run_test_hypothesis(fit, paste0(varying_cols[1], " + ", varying_cols[2]))
  }
}


test_pp_eval_func = function(fit, func) {
  # Settings
  expected_colnames = c(
    fit$pars$x,
    fit$pars$trials,
    na.omit(unique(fit$.other$ST$cp_group_col)),  # varying effects
    as.character(substitute(func)), "error", "Q2.5", "Q97.5"  # substitute-stuff just gets the func name as string
  )
  if (length(fit$pars$arma) > 0 || as.character(substitute(func)) == "residuals")
    expected_colnames = c(expected_colnames, fit$pars$y)

  # Run and test
  result = try(func(fit), silent = TRUE)
  if (inherits(result, "try-error")) {
    error_message = attr(result, "condition")$message
    expected_error_poisson = "Problem with \`mutate\\(\\)\` input"  # OK: a side-effect of the small data and short sampling.
    testthat::expect_true(fit$family$family != "poisson" | stringr::str_starts(error_message, expected_error_poisson))  # Only test message for poisson
  } else {
    testthat::expect_true(is.data.frame(result))
    testthat::expect_equal(nrow(result), nrow(fit$data))  # Returns same number of rows as data
    testthat::expect_true(sum(is.na(result)) == 0)  # No missing values)
    testthat::expect_true(dplyr::setequal(colnames(result), expected_colnames))  # Exactly these columns regardless of order
    testthat::expect_true(all(result[, fit$pars$x] == fit$data[, fit$pars$x]))  # Output should have same order as input
  }
}

test_pp_eval = function(fit) {
  # Test pp_eval
  test_pp_eval_func(fit, fitted)
  test_pp_eval_func(fit, predict)
  test_pp_eval_func(fit, residuals)

  # Test the other arguments. Inside "try" without further checking because such errors should be caught by the above.
  result_more = try(fitted(
    fit,
    newdata = fit$data[sample(nrow(fit$data), 3), ],
    summary = FALSE,
    probs = c(0.1, 0.5, 0.999),
    prior = TRUE,
    nsamples = 2,
    arma = FALSE
  ), silent = TRUE)

  if (is.data.frame(result_more)) {
    #testthat::expect_true(nrow(result_more) == nrows * nsamples * 2)  # nrows * nsamples * nchains
    testthat::expect_true(sum(is.na(result_more)) == 0)

    expected_colnames_more = c(
      # Tidybayes stuff
      ".chain", ".iteration", ".draw",

      # Model parameters
      fit$pars$population,
      fit$pars$varying,

      # Predictors
      fit$pars$x,
      fit$pars$trials,
      na.omit(unique(fit$.other$ST$cp_group_col)),  # varying effects
      "data_row",
      "fitted"
    )

    is_equal = dplyr::setequal(colnames(result_more), expected_colnames_more)  # Exactly these columns regardless of order
    testthat::expect_true(is_equal)
  }

  # Test pp_check
  if (length(fit$pars$varying) > 0) {
    varying_col = na.omit(fit$.other$ST$cp_group_col)[1]  # Just use the first column
    pp_default = try(pp_check(fit, facet_by = varying_col, nsamples = 2), silent = TRUE)
  } else {
    pp_default = try(pp_check(fit, nsamples = 2), silent = TRUE)
  }

  if (inherits(pp_default, "try-error")) {
    error_message = attr(pp_default, "condition")$message
    expected_error_poisson = "Problem with \`mutate\\(\\)\` input"  # OK: a side-effect of the small data and short sampling.
    testthat::expect_true(fit$family$family != "poisson" || stringr::str_starts(error_message, expected_error_poisson))  # Only test message for poisson
  } else {
    testthat::expect_true(ggplot2::is.ggplot(pp_default))
  }
}



# Rutine for testing a list of erroneous models
test_bad = function(models, ...) {
  for (model in models) {
    test_name = paste0(as.character(substitute(models)), ":
    ", paste0(model, collapse=", "))

    testthat::test_that(test_name, {
      testthat::expect_error(test_runs(model, sample = FALSE, ...))  # should err before sampling
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

