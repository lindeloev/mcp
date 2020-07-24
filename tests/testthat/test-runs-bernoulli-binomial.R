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
  list(y | trials(N) ~ 1 + sigma(1)),

  # Weights not implemented yet
  list(y | trials(N) + weights(weights_ok) ~ 1)
)

test_bad(bad_binomial,
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

test_good(good_binomial,
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
  list(y_bern ~ 1 + sigma(1)),

  # Weights not implemented yet
  list(y | trials(N) + weights(weights_ok) ~ 1)
)

test_bad(bad_bernoulli,
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

test_good(good_bernoulli,
          data = data_binomial,
          family = bernoulli())

