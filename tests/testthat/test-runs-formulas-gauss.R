#################
# TEST RESPONSE #
#################
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

test_bad(bad_y)


good_y = list(
  list(y ~ 1),  # Regular
  list(y ~ 1,  # Explicit and implicit y and cp
       y ~ 1 ~ 1,
       rel(1) + (1|id) ~ rel(1) + x,
       ~ 1),
  list(ok_y ~ 1)  # decimal y
)

test_good(good_y)



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
       1 + (1|id) ~ rel(1))  # With multiple segments and functions and varying
)

test_good(good_weights)


###################
# TEST INTERCEPTS #
###################
bad_intercepts = list(
  list(y ~ rel(0)),  # rel(0) not supported
  list(y ~ rel(1)),  # Nothing to be relative to here
  list(y ~ 2),  # 2 not supported
  list(y ~ 1,
       ~ rel(0))  # rel(0) not supported
)

test_bad(bad_intercepts)


good_intercepts = list(
  #list(y ~ 0),  # would be nice if it worked, but mcmc.list does not behave well with just one variable
  list(ok_y ~ 1),  # y can be called whatever
  list(y ~ 0,  # Multiple segments
       ~ 1,
       ~ 0,
       ~ 1),
  list(y ~ 1,  # Chained relative intercepts
       ~ rel(1),
       ~ rel(1))
)

test_good(good_intercepts)


###############
# TEST SLOPES #
###############
bad_slopes = list(
  list(y ~ rel(x)),  # Nothing to be relative to
  list(y ~ x + y),  # Two slopes
  list(y ~ x,  # Two slopes
       ~ y),
  list(y ~ 1,  # Relative slope after no slope
       ~ rel(x)),
  list(y ~ bad_x_char),  # not numeric x
  list(y ~ bad_x_factor),  # not numeric x
  list(y ~ 1,
       ~ log(x)),  # should fail explicitly because negative x
  list(y ~ 1,
       ~ sqrt(x))  # should fail explicitly because negative x
)

test_bad(bad_slopes)



good_slopes = list(
  list(y ~ 0 + x),  # Regular
  list(y ~ 0 + x,  # Multiple on/off
       ~ 0,
       ~ 1 + x),
  list(y ~ x,  # Chained relative slopes
       ~ 0 + rel(x),
       ~ rel(x)),
  list(y ~ 0 + x + I(x^2) + I(x^3),  # Test "non-linear" x
       ~ 0 + exp(x) + abs(x),
       ~ 0 + sin(x) + cos(x) + tan(x)),
  list(y ~ ok_x)  # alternative x
)

test_good(good_slopes, par_x = NULL)



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

test_bad(bad_cps)


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

test_good(good_cps)
