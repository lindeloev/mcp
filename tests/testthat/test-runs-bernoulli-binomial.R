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

test_bad(bad_binomial,
         data = data_binomial,
         family = binomial())


good_binomial = list(
  list(y | trials(N) ~ 1),  # one segment
  list(y | trials(N) ~ 1 + x,  # specified multiple times
       y | trials(N) ~ 1 ~ 1 + x,
       ~ 0),
  list(y | trials(N) ~ 1,  # With varying
       1 + (1|id) ~ 1),
  list(y | trials(N) ~ 1 + ar(1) + ma(1)),
  list(y | trials(N) ~ 1,
       1 ~ N),  # N can be both trials and slope
  list(y | trials(N) + weights(weights_ok) ~ 1)  # With weights
)

test_good(good_binomial,
          data = data_binomial,
          family = binomial())

test_that("binomial and Bernoulli links share weakly regularizing defaults", {
  families = list(
    mcpfamily(binomial("logit")),
    mcpfamily(binomial("probit")),
    bernoulli("logit"),
    bernoulli("probit")
  )
  expected = c(
    "dt(0, 1.5, 3)",
    "dt(0, 1.5, 3)",
    "dt(0, 1.5 / predictor_scale(), 3)"
  )

  for (family in families)
    expect_equal(family$default_prior$prior, expected)
})

test_that("binomial responses cannot exceed trials", {
  invalid_data = data_binomial
  invalid_data$y[1] = invalid_data$N[1] + 1

  expect_error(
    mcp(
      list(y | trials(N) ~ 1),
      data = invalid_data,
      family = binomial(),
      par_x = "x",
      sample = FALSE
    ),
    "responses in 'y' cannot exceed trials in 'N'. Found invalid data in row(s): 1.",
    fixed = TRUE
  )
})




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

  # Bernoulli does not take trials
  list(y | trials(N) + weights(weights_ok) ~ 1)
)

test_bad(bad_bernoulli,
         data = data_binomial,
         family = bernoulli())


good_bernoulli = list(
  list(y_bern ~ 1),  # one segment
  list(y_bern ~ 1 + x,  # specified multiple times
       y_bern ~ 1 ~ 1 + x,
       1 ~ 0),
  list(y_bern ~ 1,  # With varying
       1 + (1|id) ~ 1),
  list(y_bern | weights(weights_ok) ~ 1),  # With weights
  list(y_bern ~ 1 + ar(1) + ma(1))  # With AR and MA
)

test_good(good_bernoulli,
          data = data_binomial,
          family = bernoulli())

