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
  bad_id = rnorm(5),  # decimal numbers

  weights_ok = c(0.1, 1, 2, 1, 1),
  weights_bad = c(-0.1, 1, 2, 1, 1)  # With negative
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
  id = c("a", "b", "c", "d", "e"),

  weights_ok = c(0.1, 1, 2, 1, 1)  # Actually not OK since it's not implemented yet
)
