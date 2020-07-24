bad_poisson = list(
  # Misspecification of y and trials
  list(y | trials(N) ~ 1),  # bad response format
  list(y ~ 1 + x,
       y | trials(N) ~ 1 ~ 1),  # misspecification in later segment

  # Bad data
  list(y_bad_numeric ~ 1),

  # Does not work with sigma
  list(y ~ 1 + sigma(1)),

  # Does not work with weights
  list(y | weights(weights_ok) ~ 1)
)

test_bad(bad_poisson,
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

test_good(good_poisson,
          data = data_binomial,
          family = poisson())
