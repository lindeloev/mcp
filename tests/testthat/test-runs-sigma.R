#################
# TEST VARIANCE #
#################
bad_variance = list(
  list(y ~ 1 + sigma(q))  # variable does not exist
)

test_bad(bad_variance)


good_variance = list(
  list(y ~ 1 + sigma(1)),
  list(y ~ 1 + sigma(x + I(x^2))),
  list(y ~ 1 + sigma(1 + sin(x))),
  list(y ~ 1,
       1 + (1|id) ~ 1 + I(x^2) + sigma(1 + x)),  # Test with varying change point and more mcp stuff
  list(y | weights(weights_ok) ~ 1 + sigma(1 + x),  # With weights
       ~ 0 + sigma(1 + x))
)

test_good(good_variance)
