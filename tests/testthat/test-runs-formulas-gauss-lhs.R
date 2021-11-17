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
       1 + (1|bad_id) ~ 1)  # decimal group
)

test_bad(bad_cps)




good_cps = list(
  list(y ~ 0 + x,  # Regular cp
       1 ~ 1),
  list(y ~ 1,  # Implicit cp
       ~ 1,
       ~ 0),
  list(y ~ 0,  # Varying
       1 + (1|id) ~ 1),
  list(y ~ 0,  # Chained varying cps
       y ~ 1 ~ 1,
       1 + (1|id) ~ 0,
       1 + (1|id) ~ 0,
       ~ x),
  list(y ~ 1,
       (1|id) ~ 0),  # Intercept is implicit. I don't like it, but OK.
  list(y ~ 1,
       1 + (1|id) ~ 1,
       1 + (1|ok_id_integer) ~ 1,  # multiple groups and alternative data
       1 + (1|ok_id_factor) ~ 1)  # alternative group data
)

test_good(good_cps)



# Test detection of par_x
bad_par_x = list(
  list(y ~ 0),  # no par_x
  list(y ~ 1),  # no par_x
  list(y ~ 1 + bad_x_char,
       ~ 0 + bad_x_factor),  # only invalid par_x
  list(y ~ bad_x_char),  # Has to be continuous
  list(y ~ bad_x_factor),  # Has to be continuous
  list(y ~ x + ok_x),  # multiple in one segment
  list(y ~ x,
       ~ ok_x)  # Multiple in two segments
)

test_bad(bad_par_x, par_x = NULL)


good_par_x = list(
  list(y ~ 1 + x),
  list(y ~ 0,
       ~ 1,
       ~ 0 + ok_x)
)

test_good(good_par_x, par_x = NULL)



################
# TEST WEIGHTS #
################
bad_weights = list(
  list(y + weights(weights_ok) ~ 1),  # weights added
  list(weights(y) ~ 1),  # just wrong :-)
  list(y | weights_ok ~ 1),  # Has to be weights(weights_ok)
  list(y | weights(weights_bad) ~ 1),  # Bad weights
  list(y | weights(weights_ok) ~ 1,
       y | weights(weights_bad) ~ 1 ~ 1)  # Different weights
)

test_bad(bad_weights)

good_weights = list(
  list(y | weights(weights_ok) ~ 1),  # Regular
  list(y | weights(weights_ok) ~ 1,
       ~ 1 + x + I(x^2),
       1 + (1|id) ~ 1)  # With multiple segments and functions and varying
)

test_good(good_weights)
