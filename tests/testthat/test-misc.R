# Link functions
test_that("Link functions", {
  expect_equal(ilogit(1), 0.73105858)
  expect_equal(logit(0.7), 0.84729786)
  expect_equal(phi(1), 0.84134475)
  expect_equal(probit(0.7), 0.52440051)
})

test_that("families", {
  expect_true(is.mcpfamily(bernoulli()))
  expect_false(is.mcpfamily(gaussian()))
})

# Code
test_that("formula tools", {
  expect_equal(sd_to_prec("dnorm(5, 1)"), "dnorm(5, 1/(1)^2) ")
  expect_equal(sd_to_prec("dt(3,2,1)"), "dt(3, 1/(2)^2, 1) ")
})


########################
# MCPFIT CLASS-METHODS #
########################
# Test on new fit
demo_settings = mcp_example("demo", sample = FALSE)
demo_fit2 = mcp(demo_settings$model, demo_settings$data, adapt = 2500, iter = 4000, cores = 3)

test_that("Simple mcpfit methods", {
  expect_equal(niterations(demo_fit2), 12000)
  expect_equal(nchains(demo_fit), 3)

  expect_true(is.mcpfit(demo_fit))
  expect_false(is.mcpfit(mtcars))

})

# hypothesis()
test_that("hypothesis()", {
  actual_hypothesis1 = hypothesis(demo_fit2, "cp_1 > 27")
  expected_hypothesis1 = data.frame(hypothesis = "cp_1 - 27 > 0", mean = 4.37, lower = -3.12, upper = 12.34, p = 0.87, BF = 7.3)
  expect_equal(actual_hypothesis1, expected_hypothesis1, tolerance = 0.2)

  actual_hypothesis2 = hypothesis(demo_fit2, "(cp_1 > 27 | cp_1 < 25) & time_3 > -0.2")
  expected_hypothesis2 = data.frame(hypothesis = "(cp_1 > 27 | cp_1 < 25) & time_3 > -0.2", mean = NA, lower = NA, upper = NA, p = 0.166, BF = 0.199)
  expect_equal(actual_hypothesis2, expected_hypothesis2, tolerance = 0.03)
})
