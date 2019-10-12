###############################
# TEST SEGMENT SPECIFICATIONS #
###############################

# TO DO: CHECK VALUES TOO
check_empty = function(fit, segments) {
  # Correct types
  is.null(fit$samples) &
    is.null(fit$loglik) &
    is.null(fit$loo) &
    is.null(fit$waic) &
    is.null(fit$data) &
    is.list(fit$prior) &
    is.list(fit$segments) &
    is.list(fit$pars) &
    is.character(fit$jags_code) &
    is.function(fit$func_y) &
    all.equal(fit$segments, segments)
}



test_that("bad intercepts", {
  bad_intercepts = list(
    list(y ~ rel(0)),
    list(y ~ rel(1)),
    list( ~ 1),
    list(y ~ 2),
    list(y ~ 1, 1 ~ rel(0))  # Not (yet) supported
  )

  for(segments in bad_intercepts) {
    expect_error(mcp(
      segments = segments,
      param_x = "x",
      sample = FALSE
    ), info = paste0("segments: ", segments))
  }
})


test_that("bad slopes", {
  bad_slopes = list(
    list(y ~ rel(x)),
    list(y ~ x + y),
    list(y ~ x, 1 ~ y),
    list(y ~ x, 0 ~ 1),
    list(y ~ x, q ~ 1)
  )

  for(segments in bad_slopes) {
    expect_error(mcp(
      segments = segments,
      sample = FALSE
    ), info = paste0("segments: ", segments))
  }
})


test_that("good intercepts", {
  good_intercepts = list(
    list(y ~ 0),
    list(x ~ 1),
    list(stuff ~ 0, 1 ~ 1, 1 ~ 0, 1 ~ 1)
  )

  for(segments in good_intercepts) {
    fit = mcp(
      segments = segments,
      param_x = "horse",
      sample = FALSE
    )
    expect_true(
      check_empty(fit, segments),
      info = paste0("segments: ", segments)
    )
  }
})


test_that("good slopes", {
  good_slopes = list(
    list(y ~ x),
    list(stuff ~ goat, 1 ~ 0, 1 ~ 1 + goat)
  )

  for(segments in good_slopes) {
    fit = mcp(
      segments = segments,
      sample = FALSE
    )
    expect_true(
      check_empty(fit, segments),
      info = paste0("segments: ", segments)
    )
  }
})


