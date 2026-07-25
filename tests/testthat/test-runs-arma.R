#############
# TEST ARMA #
#############
# We can assume that it will fail for the same mis-specifications on the formula
# ar(order, [formula]), since the formula runs through the exact same code as
# sigma and mu.
bad_arma = list(
  list(y ~ ar(0)),  # currently not implemented
  list(y ~ ar(-1)),  # must be positive
  list(y ~ ar(1.5)),  # Cannot be decimal
  list(y ~ ar(1) + ar(2)),  # Only one per segment
  list(y ~ ar("1")),  # Should not take strings
  list(y ~ ar(1 + x)),  # must have order
  list(y ~ ar(x)),  # must have order
  list(y ~ ma(0)),
  list(y ~ ma(-1)),
  list(y ~ ma(1.5)),
  list(y ~ ma(1) + ma(2))
)

test_bad(bad_arma)


good_arma = list(
  list(y ~ ar(1)),  # simple
  list(y ~ ma(1)),
  list(y ~ ar(1) + ma(1)),
  list(y ~ ar(3)),  # higher order
  list(y ~ ar(1, 1 + x + I(x^2) + exp(x))),  # complicated regression
  list(y ~ ar(1),
       ~ ar(2, 0 + x)),  # change in ar
  list(y ~ 1,
       ~ 0 + ar(2)),  # onset of AR
  list(y ~ 1,
       1 + (1|id) ~ 1 + I(x^2) + ar(2, 1 + x)),  # varying change point
  list(y ~ ar(1) + sigma(1 + x),
       ~ ar(2, 1 + I(x^2)) + sigma(1)),  # With sigma
  list(y ~ ar(1),
       ~ ar(2, 1)),
  list(y | weights(weights_ok) ~ 1 + ar(1),  # With weights
       ~ 0 + ar(2, 1 + x)),
  list(y ~ 0 + ar(1),
       ~ 0 + ar(2))  # mu is ~0 across all segments; only ar() gives structure
)

test_good(good_arma)


test_that("ar() and ma() require a supported default-link model", {
  data = data.frame(x = 1:4, y = 1:4)

  expect_error(
    mcp(
      list(y ~ 1 + ar(1) + ma(1)),
      data,
      family = gaussian(link = "log"),
      par_x = "x",
      sample = FALSE
    ),
    "ar() and ma() currently support gaussian(), binomial(), poisson(), and negbinomial() with their default links; bernoulli() and non-default links are not supported.",
    fixed = TRUE
  )
})
