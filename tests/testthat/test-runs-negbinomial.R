bad_negbinomial = list(
  list(y | trials(N) ~ 1),
  list(y ~ 1 + sigma(1)),
  list(y | weights(weights_ok) ~ 1),
  list(y_bad_numeric ~ 1)
)

test_bad(
  bad_negbinomial,
  data = data_binomial,
  family = negbinomial()
)


good_negbinomial = list(
  list(y ~ 1 + x),
  list(
    y ~ 1 + x,
    ~ 1 + x
  ),
  list(y ~ 1 + x + shape(1 + x)),
  list(
    y ~ 1 + x + shape(1),
    ~ 0 + x + shape(1)
  ),
  list(y ~ 1 + ar(1) + ma(1))
)

test_good(
  good_negbinomial,
  data = data_binomial,
  family = negbinomial()
)


test_that("negative-binomial links are explicit and currently log-only", {
  expect_error(negbinomial(link = "identity"), '`link` must be one of "log"')
  expect_error(negbinomial(link_shape = "identity"), '`link_shape` must be one of "log"')
})
