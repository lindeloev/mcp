#########
# SETUP #
#########

# Samples and checks data structure.
# Meant to be used with testthat::expect_true()
test_data = data.frame(
  # y should be continuous
  y = 1:5,
  ok_y = rnorm(5),  # test underscore and decimals
  bad_y_char = c("a", "b", "c", "d", "e"),
  bad_y_factor = factor(1:5),

  # x should be continuous
  x = 1:5,
  ok_x = rnorm(5),  # test underscore and decimals
  bad_x_char = c("a", "b", "c", "d", "e"),
  bad_x_factor = factor(1:5),

  # change point should be categorical-ish
  id = c("a", "b", "c", "d", "e"),
  ok_id_factor = factor(c(-3, 0, 5, 9, 1.233243)),
  ok_id_integer = -2:2,  # interval
  bad_id = rnorm(5)  # decimal numbers
)

check_empty = function(segments, par_x = NULL) {
  # Sample
  fit = mcp(
    segments = segments,
    data = test_data,
    par_x = par_x,
    sample = FALSE
  )

  # Do it again, but this time on a tibble
  fit = mcp(
    segments = segments,
    data = tibble::as_tibble(test_data),
    par_x = par_x,
    sample = FALSE
  )

  # Check for correct types in data structure
  is.null(fit$samples) &
    is.null(fit$loglik) &
    is.null(fit$loo) &
    is.null(fit$waic) &
    is.list(fit$data) &
    is.list(fit$prior) &
    is.list(fit$segments) &
    is.list(fit$pars) &
    is.character(fit$pars$population) &
    (is.character(fit$pars$varying) | is.null(fit$pars$varying)) &
    is.character(fit$pars$x) &
    is.character(fit$pars$y) &
    is.character(fit$jags_code) &
    is.function(fit$func_y) &
    is.list(fit$.other) &
    all.equal(fit$segments, segments)
}

# Samples a model. Meant to be used with testthat::expect_error()
check_error = function(segments) {
  mcp(
    segments = segments,
    data = test_data,
    par_x = "x",
    sample = FALSE
  )
}


##########
# TEST Y #
##########

test_that("Bad y", {
  bad_y = list(
    list(  ~ 1),  # No y
    list((1|id) ~ 1),  # y cannot be varying
    list(1 ~ 1),  # 1 is not y
    list(y ~ 1,  # Two y
         a ~ 1 ~ 1),
    list(y ~ 1,  # Intercept y
         1 ~ 1 ~ 1),
    list(bad_y_char ~ 1),  # Character y
    list(bad_y_factor ~ 1)  # Factor y
  )

  for (segments in bad_y) {
    expect_error(
      check_error(segments),
      info = paste0("segments: ", segments)
    )
  }
})

test_that("Good y", {
  bad_y = list(
    list(y ~ 1),  # Regular
    list(y ~ 1,  # Explicit and implicit y
         y ~ 1 ~ 1,
         rel(1) + (1|id) ~ rel(1) + x,
         ~ 1),
    list(ok_y ~ 1)  # decimal y
  )

  for (segments in bad_y) {
    expect_true(
      check_empty(segments, par_x = "x"),
      info = paste0("segments: ", segments)
    )
  }
})


###################
# TEST INTERCEPTS #
###################

test_that("bad intercepts", {
  bad_intercepts = list(
    list(y ~ rel(0)),  # rel(0) not supported
    list(y ~ rel(1)),  # Nothing to be relative to here
    list(y ~ 2),  # 2 not supported
    list(y ~ 1,
         1 ~ rel(0))  # rel(0) not supported
  )

  for (segments in bad_intercepts) {
    expect_error(
      check_error(segments),
      info = paste0("segments: ", segments)
    )
  }
})

test_that("good intercepts", {
  good_intercepts = list(
    list(y ~ 0),
    list(ok_y ~ 1),  # y can be called whatever
    list(y ~ 0,  # Multiple segments
         1 ~ 1,
         1 ~ 0,
         1 ~ 1),
    list(y ~ 1,  # Chained relative intercepts
         1 ~ rel(1),
         1 ~ rel(1))
  )

  for (segments in good_intercepts) {
    expect_true(
      check_empty(segments, par_x = "x"),
      info = paste0("segments: ", segments)
    )
  }
})



###############
# TEST SLOPES #
###############

test_that("bad slopes", {
  bad_slopes = list(
    list(y ~ rel(x)),  # Nothing to be relative to
    list(y ~ x + y),  # Two slopes
    list(y ~ x,  # Two slopes
         1 ~ y),
    list(y ~ 1,  # Relative slope after no slope
         1 ~ rel(x)),
    list(y ~ bad_x_char),  # not numeric x
    list(y ~ bad_x_factor)  # not numeric x
  )

  for (segments in bad_slopes) {
    expect_error(
      check_error(segments),
      info = paste0("segments: ", segments)
    )
  }
})

test_that("good slopes", {
  good_slopes = list(
    list(y ~ 0 + x),  # Regular
    list(y ~ 0 + x,  # Multiple on/off
         1 ~ 0,
         1 ~ 1 + x),
    list(y ~ x,  # Chained relative slopes
         1 ~ 0 + rel(x),
         1 ~ rel(x)),
    list(y ~ ok_x)  # alternative x
  )

  for (segments in good_slopes) {
    expect_true(
      check_empty(segments),
      info = paste0("segments: ", segments)
    )
  }
})



######################
# TEST CHANGE POINTS #
######################

test_that("bad change points", {
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

  for (segments in bad_cps) {
    expect_error(
      check_error(segments),
      info = paste0("segments: ", segments)
    )
  }
})

test_that("good change points", {
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

  for (segments in good_cps) {
    expect_true(
      check_empty(segments, par_x = "x"),
      info = paste0("segments: ", segments)
    )
  }
})
