test_that("S3 methods work on archetypal Gaussian fit", {
  test_runs(list(y ~ 1 + x, ~ 0 + x), test_s3 = TRUE)
})

test_that("S3 methods work on archetypal group-level-effects fit", {
  test_runs(list(y ~ 1 + (1 | id), ~ 0 + x), test_s3 = TRUE)
})

test_that("S3 methods work on archetypal AR fit", {
  test_runs(list(y ~ 1 + ar(1), ~ 0 + ar(2, 1 + x)), test_s3 = TRUE)
})

test_that("S3 methods work on archetypal Binomial fit", {
  test_runs(
    list(y | trials(N) ~ 1 + x, ~ 1),
    data = data_binomial,
    family = binomial(),
    test_s3 = TRUE
  )
})

test_that("S3 methods work on archetypal Poisson fit", {
  test_runs(
    list(y ~ 1 + x, ~ 1),
    family = poisson(),
    test_s3 = TRUE
  )
})
