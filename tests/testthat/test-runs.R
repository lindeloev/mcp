
options(mc.cores = 1)

############
# DATASETS #
############

# Samples and checks data structure.
# Meant to be used with testthat::expect_true()
data_gauss = data.frame(
  # y should be continuous
  y = 1:5,
  ok_y = rnorm(5),  # test underscore and decimals
  bad_y_char = c("a", "b", "c", "d", "e"),
  bad_y_factor = factor(1:5),

  # x should be continuous
  x = -1:3,
  ok_x = rnorm(5),  # test underscore and decimals
  bad_x_char = c("a", "b", "c", "d", "e"),
  bad_x_factor = factor(1:5),

  # varying effects should be categorical-ish
  id = c("a", "b", "c", "d", "e"),
  ok_id_factor = factor(c(-3, 0, 5, 9, 1.233243)),  # It's a factor, so decimals are OK
  ok_id_integer = -2:2,  # interval
  bad_id = rnorm(5)  # decimal numbers
)

# Only needs to test binomial-specific stuff
data_binomial = data.frame(
  # y should be a natural number > 0
  y = c(1, 0, 100, 3, 5),
  y_bad_numeric = c(-1, 5.1, 10, 3, 5),  # negative, decimal,

  y_bern = c(0, 1, 0, 1, 1),

  # trials should be a natural number 0 <= N <= y
  N = c(1, 1, 100, 6, 10),
  N_bad_numeric = c(-1, 1.1, 99, 6, 10),  # smaller than y, decimal, negative
  N_bad_factor = factor(c(1, 0, 50, 6, 10)),
  N_bad_char = c("1", "1", "100", "6", "10"),

  # x
  x = -1:3,

  # Varying effects
  id = c("a", "b", "c", "d", "e")
)


##################
# TEST FUNCTIONS #
##################

test_runs = function(segments,
                    data = data_gauss,
                    prior = list(),
                    family = gaussian(),
                    par_x = "x",
                    sample = TRUE) {

  # Without sampling, on a data.frame.
  empty = mcp(
    segments = segments,
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
    testthat::expect_true(is.list(empty$segments), segments)
    testthat::expect_true(is.data.frame(empty$data), segments)
    testthat::expect_true(is.list(empty$prior), segments)
    testthat::expect_true(class(empty$family) == "family", segments)
    testthat::expect_true(is.null(empty$samples), segments)
    testthat::expect_true(is.null(empty$loglik), segments)
    testthat::expect_true(is.null(empty$loo), segments)
    testthat::expect_true(is.null(empty$waic), segments)
    testthat::expect_true(is.list(empty$pars), segments)
    testthat::expect_true(is.character(empty$pars$population), segments)
    testthat::expect_true((is.character(empty$pars$varying) | is.null(empty$pars$varying)), segments)
    testthat::expect_true(is.character(empty$pars$x), segments)
    testthat::expect_true(is.character(empty$pars$y), segments)
    testthat::expect_true(is.character(empty$jags_code), segments)
    testthat::expect_true(is.function(empty$simulate), segments)
    testthat::expect_true(is.list(empty$.other), segments)

    # Should work for tibbles as well. So do this sometimes
    if (rbinom(1, 1, 0.5) == 1)
      data = tibble::as_tibble(data)

    # Capture (expected) messages and warnings
    quiet_out = purrr::quietly(mcp)(  # Do not print to console
      segments = segments,
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
    fit <<- quiet_out$result

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
      # (get_samples checks for NULL)
      if (col == "mcmc_prior")
        fit$mcmc_post = NULL

      # Check that samples are the correct format
      testthat::expect_true(is.list(fit[[col]]), segments)
      testthat::expect_true(coda::is.mcmc(fit[[col]][[1]]), segments)
      testthat::expect_true(all(fit$pars$population %in% colnames(fit[[col]][[1]])))

      # Test mcpfit functions
      varying_cols = na.omit(fit$.other$ST$cp_group_col)
      test_summary(fit, varying_cols)
      test_plot(fit, varying_cols)  # default plot
      test_plot_pars(fit)  # bayesplot call
    }
  }
}


# Tests if summary(fit) and ranef(fit) work as expected
test_summary = function(fit, varying_cols) {
  summary_cols = c('name','mean','lower','upper','Rhat','n.eff','ts_se')
  result = purrr::quietly(summary)(fit)$result  # Do not print to console
  output = purrr::quietly(summary)(fit)$output  # Do not print to console
  testthat::expect_true(all(colnames(result) %in% summary_cols))  # All columns
  testthat::expect_true(all(result$name %in% fit$pars$population))  # All parameters

  # If there are varying effects
  if (length(varying_cols) > 0) {
    testthat::expect_match(output, "ranef\\(")  # noticed about varying effects
    varying = ranef(fit)  # Do not print to console
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
    gg = try(plot(fit, facet_by = varying_cols[1], q_fit = q_fit, q_predict = q_predict), silent = TRUE)  # just take the first
  } else {
    gg = try(plot(fit, q_fit = q_fit, q_predict = q_predict), silent = TRUE)
  }
  # Is it a ggplot or a known error?
  if (inherits(gg, "try-error")) {
    # (the error is an artefact of very small test data --> wide posteriors.)
    if (fit$family$family == "poisson") {
      expected_error = "Modelled extremely large value"
    } else if (any(stringr::str_detect(fit$pars$sigma, "^sigma_.*_.*$"))) {  # for slopes on sigma
      expected_error = "Modelled negative sigma"
    } else {
      expected_error = ">>>>do_not_expect_any_errors<<<<<"
    }
    expect_true(any(stringr::str_starts(attr(gg, "condition")$message, expected_error)))
  } else {
    testthat::expect_s3_class(gg, c("gg", "ggplot"))
  }
}

# Test plot() calls to bayesplot
test_plot_pars = function(fit) {
  gg = plot_pars(fit, type = "dens_overlay")
  testthat::expect_s3_class(gg, c("gg", "ggplot"))
}



test_hypothesis = function(fit) {
  # Function to test both directional and point hypotheses
  run_test_hypothesis = function(fit, base) {
    hypotheses = c(
      paste0(base, " > 1"),  # Directional
      paste0(base, " = -1")  # Savage-Dickey (point)
    )
    result = hypothesis(fit, hypotheses)  # Do not print to console
    testthat::expect_true(is.data.frame(result) & nrow(result) == 2)
  }

  # Test single pop effect
  run_test_hypothesis(fit, fit$pars$population[1])

  # Test multiple pop effect
  if (length(fit$pars$population) > 1)
    run_test_hypothesis(fit, paste0(fit$pars$population[1] , " + ", fit$pars$population[2]))

  # Varying
  if (!is.null(fit$pars$varying)) {
    mcmc_vars = colnames(get_samples(fit)[[1]])
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



# Ruitine for testing a list of erroneous segments
test_bad = function(segments_list, title, ...) {
  for (segments in segments_list) {
    test_name = paste0(title, ":
    ", paste0(segments, collapse=", "))

    testthat::test_that(test_name, {
      testthat::expect_error(test_runs(segments, sample = FALSE, ...))  # should err before sampling
    })
  }
}


# Routine for testing a list of good segments
test_good = function(segments_list, title, ...) {
  for (segments in segments_list) {
    test_name = paste0(title, ":
    ", paste0(segments, collapse=", "))

    testthat::test_that(test_name, {
      test_runs(segments, ...)
    })
  }
}



###############
# TEST PRIORS #
###############
good_prior_segments = list(
  y ~ 1 + x,
  1 + (1|id) ~ rel(1) + rel(x),
  rel(1) ~ 0
)
good_prior = list(
  list(
    int_2 = "int_1",
    cp_1 = "dnorm(3, 10)",
    x_2 = "-0.5"
  )
)

for (prior in good_prior) {
  test_name = paste0("Good priors: ", paste0(prior, collapse=", "))
  testthat::test_that(test_name, {
    test_runs(good_prior_segments, prior = prior)
  })
}



##########
# TEST Y #
##########
bad_y = list(
  list( ~ 1),  # No y
  list((1|id) ~ 1),  # y cannot be varying
  list(1 ~ 1),  # 1 is not y
  list(y ~ 1,  # Two y
       a ~ 1 ~ 1),
  list(y ~ 1,  # Intercept y
       1 ~ 1 ~ 1),
  list(bad_y_char ~ 1),  # Character y
  list(bad_y_factor ~ 1)  # Factor y
)

test_bad(bad_y, "Bad y")


good_y = list(
  list(y ~ 1),  # Regular
  list(y ~ 1,  # Explicit and implicit y and cp
       y ~ 1 ~ 1,
       rel(1) + (1|id) ~ rel(1) + x,
       ~ 1),
  list(ok_y ~ 1)  # decimal y
)

test_good(good_y, "Good y")



###################
# TEST INTERCEPTS #
###################
bad_intercepts = list(
  list(y ~ rel(0)),  # rel(0) not supported
  list(y ~ rel(1)),  # Nothing to be relative to here
  list(y ~ 2),  # 2 not supported
  list(y ~ 1,
       1 ~ rel(0))  # rel(0) not supported
)

test_bad(bad_intercepts, "Bad intercepts")


good_intercepts = list(
  #list(y ~ 0),  # would be nice if it worked, but mcmc.list does not behave well with just one variable
  list(ok_y ~ 1),  # y can be called whatever
  list(y ~ 0,  # Multiple segments
       1 ~ 1,
       1 ~ 0,
       1 ~ 1),
  list(y ~ 1,  # Chained relative intercepts
       1 ~ rel(1),
       1 ~ rel(1))
)

test_good(good_intercepts, "Good intercepts")




###############
# TEST SLOPES #
###############

bad_slopes = list(
  list(y ~ rel(x)),  # Nothing to be relative to
  list(y ~ x + y),  # Two slopes
  list(y ~ x,  # Two slopes
       1 ~ y),
  list(y ~ 1,  # Relative slope after no slope
       1 ~ rel(x)),
  list(y ~ bad_x_char),  # not numeric x
  list(y ~ bad_x_factor),  # not numeric x
  list(y ~ 1,
       1 ~ log(x)),  # should fail explicitly because negative x
  list(y ~ 1,
       1 ~ sqrt(x))  # should fail explicitly because negative x
)

test_bad(bad_slopes, "Bad slopes")



good_slopes = list(
  list(y ~ 0 + x),  # Regular
  list(y ~ 0 + x,  # Multiple on/off
       1 ~ 0,
       1 ~ 1 + x),
  list(y ~ x,  # Chained relative slopes
       1 ~ 0 + rel(x),
       1 ~ rel(x)),
  list(y ~ 0 + x + I(x^2) + I(x^3),  # Test "non-linear" x
       1 ~ 0 + exp(x) + abs(x),
       1 ~ 0 + sin(x) + cos(x) + tan(x)),
  list(y ~ ok_x)  # alternative x
)

test_good(good_slopes, "Good slopes", par_x = NULL)



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
       rel(1) ~ 1),  # Nothing to be relative to yet
  list(y ~ 1,
       1 + (1|bad_id) ~ 1)  # decimal group
)

test_bad(bad_cps, "Bad change points")


good_cps = list(
  list(y ~ 0 + x,  # Regular cp
       1 ~ 1),
  list(y ~ 1,  # Implicit cp
       ~ 1,
       ~ 0),
  list(y ~ 0,  # Varying
       1 + (1|id) ~ 1),
  list(y ~ 0,  # Chained varying and relative cp
       y ~ 1 ~ 1,
       rel(1) + (1|id) ~ 0,
       rel(1) + (1|id) ~ 0,
       ~ x),
  list(y ~ 1,
       (1|id) ~ 0),  # Intercept is implicit. I don't like it, but OK.
  list(y ~ 1,
       1 + (1|id) ~ 1,
       1 + (1|ok_id_integer) ~ 1,  # multiple groups and alternative data
       1 + (1|ok_id_factor) ~ 1)  # alternative group data
)

test_good(good_cps, "Good change points")




#################
# TEST VARIANCE #
#################
bad_variance = list(
  list(y ~ 1 + sigma(rel(1))),  # no sigma to be relative to
  list(y ~ 1,
       y ~ 1 + sigma(rel(x))),  # no sigma slope to be relative to
  list(y ~ 1 + sigma(q))  # variable does not exist
)

test_bad(bad_variance, "Bad variance")


good_variance = list(
  list(y ~ 1 + sigma(1)),
  list(y ~ 1 + sigma(x + I(x^2))),
  list(y ~ 1 + sigma(1 + sin(x))),
  list(y ~ 1,
       ~ 0 + sigma(rel(1)),  # test relative intercept
       ~ x + sigma(x),
       ~ 0 + sigma(rel(x))),  # test relative slope
  list(y ~ 1,
      1 + (1|id) ~ rel(1) + I(x^2) + sigma(rel(1) + x))  # Test with varying change point and more mcp stuff
)

test_good(good_variance, "Good variance")




#############
# TEST ARMA #
#############
# We can assume that it will fail for the same mis-specifications on the formula
# ar(order, [formula]), since the formula runs through the exact same code as
# sigma and ct.
bad_arma = list(
  list(y ~ ar(0)),  # currently not implemented
  list(y ~ ar(-1)),  # must be positive
  list(y ~ ar(1.5)),  # Cannot be decimal
  list(y ~ ar(1) + ar(2)),  # Only one per segment
  list(y ~ ar("1")),  # Should not take strings
  list(y ~ ar(1 + x)),  # must have order
  list(y ~ ar(x))  # must have order
)

test_bad(bad_arma, "Bad ARMA")


good_arma = list(
  list(y ~ ar(1)),  # simple
  list(y ~ ar(5)),  # higher order
  list(y ~ ar(1, 1 + x + I(x^2) + exp(x))),  # complicated regression
  list(y ~ ar(1),
       ~ ar(2, 0 + x)),  # change in ar
  list(y ~ 1,
       ~ 0 + ar(2)),  # onset of AR
  list(y ~ 1,
       1 + (1|id) ~ rel(1) + I(x^2) + ar(2, rel(1) + x)),  # varying change point
  list(y ~ ar(1) + sigma(1 + x),
       ~ ar(2, 1 + I(x^2)) + sigma(1)),  # With sigma
  list(y ~ ar(1),
       ~ ar(2, rel(1)))  # Relative to no variance. Perhaps alter this behavior so it becomes illegal?
)

test_good(good_arma, "Good ARMA")






#################
# TEST BINOMIAL #
#################

bad_binomial = list(
  # Misspecification of y and trials
  list(y ~ 1),  # no trials
  list(y | N ~ 1),  # wrong format
  list(trials(N) | y ~ 1),  # Wrong order
  list(y | trials() ~ 1),  # trials missing
  list(trials(N) ~ 1),  # no y
  list(y | trials(N) ~ 1 + x,
       y | N ~ 1 ~ 1),  # misspecification in later segment

  # Bad data
  list(y_bad_numeric | trials(N) ~ 1),
  list(y | trials(N_bad_numeric) ~ 1),
  list(y | trials(N_bad_factor) ~ 1),
  list(y | trials(N_bad_char) ~ 1),

  # Does not work with sigma
  list(y | trials(N) ~ 1 + sigma(1))
)

test_bad(bad_binomial, "Bad binomial",
         data = data_binomial,
         family = binomial())


good_binomial = list(
  list(y | trials(N) ~ 1),  # one segment
  list(y | trials(N) ~ 1 + x,  # specified multiple times and with rel()
       y | trials(N) ~ 1 ~ rel(1) + rel(x),
       rel(1) ~ 0),
  list(y | trials(N) ~ 1,  # With varying
       1 + (1|id) ~ 1),
  list(y | trials(N) ~ 1 + ar(1))  # Simple AR(1)
  #list(y | trials(N) ~ 1,
  #     1 ~ N)  # N can be both trials and slope. TO DO: Fails in this test because par_x = "x"
)

test_good(good_binomial, "Good binomial",
          data = data_binomial,
          family = binomial())




##################
# TEST BERNOULLI #
##################
# This is rather short since most is tested via binomial
bad_bernoulli = list(
  # Misspecification of y and trials
  list(y_bern | trials(N) ~ 1),  # trials
  list(y_bern ~ 1 + x,
       y_bern | trials(N) ~ 1 ~ 1),  # misspecification in later segment

  # Bad data
  list(y_bad_numeric ~ 1),
  list(y ~ 1),  # binomial response

  # Does not work with sigma
  list(y_bern ~ 1 + sigma(1))
)

test_bad(bad_bernoulli, "Bad Bernoulli",
         data = data_binomial,
         family = bernoulli())


good_bernoulli = list(
  list(y_bern ~ 1),  # one segment
  list(y_bern ~ 1 + x,  # specified multiple times and with rel()
       y_bern ~ 1 ~ rel(1) + rel(x),
       rel(1) ~ 0),
  list(y_bern ~ 1,  # With varying
       1 + (1|id) ~ 1)
)

test_good(good_bernoulli, "Good Bernoulli",
          data = data_binomial,
          family = bernoulli())



################
# TEST POISSON #
################
# Like binomial, but without the trials()

bad_poisson = list(
  # Misspecification of y and trials
  list(y | trials(N) ~ 1),  # bad response format
  list(y ~ 1 + x,
       y | trials(N) ~ 1 ~ 1),  # misspecification in later segment

  # Bad data
  list(y_bad_numeric ~ 1),

  # Does not work with sigma
  list(y ~ 1 + sigma(1))
)

test_bad(bad_poisson, "Bad Poisson",
         data = data_binomial,
         family = poisson())


good_poisson = list(
  list(y ~ 1),  # one segment
  list(y ~ 1 + x,  # specified multiple times and with rel()
       y  ~ 1 ~ rel(1) + rel(x),
       rel(1) ~ 0),
  list(y ~ 1,  # With varying
       1 + (1|id) ~ 1),
  list(y ~ 1 + ar(1),
       ~ 1 + x + ar(2, 1 + x + I(x^3)))
)

test_good(good_poisson, "Good Poisson",
          data = data_binomial,
          family = poisson())
