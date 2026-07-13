test_that("AR and MA terms are parsed as separate model components", {
  data = data.frame(x = 1:6, y = 1:6)
  rhs_table = get_rhs_table(
    list(y ~ 1 + x + ar(2, 1 + x) + ma(1, 0 + x)),
    data,
    mcpfamily(gaussian()),
    par_x = "x"
  )

  expect_setequal(unique(rhs_table$dpar), c("mu", "sigma", "ar", "ma"))
  expect_equal(unique(rhs_table$order[rhs_table$dpar == "ar"]), 1:2)
  expect_equal(unique(rhs_table$order[rhs_table$dpar == "ma"]), 1)
  expect_equal(sum(rhs_table$dpar == "ar"), 4)
  expect_equal(sum(rhs_table$dpar == "ma"), 1)
  expect_true(all(grepl("^ar[12]_", rhs_table$code_name[rhs_table$dpar == "ar"])))
  expect_true(all(grepl("^ma1_", rhs_table$code_name[rhs_table$dpar == "ma"])))
})


test_that("AR and MA each allow one term per segment", {
  data = data.frame(x = 1:6, y = 1:6)
  family = mcpfamily(gaussian())

  expect_error(
    get_rhs_table(list(y ~ ar(1) + ar(2)), data, family, par_x = "x"),
    "Only one of these allowed per segment"
  )
  expect_error(
    get_rhs_table(list(y ~ ma(1) + ma(2)), data, family, par_x = "x"),
    "Only one of these allowed per segment"
  )
  expect_no_error(
    get_rhs_table(list(y ~ ar(1) + ma(1)), data, family, par_x = "x")
  )
})


test_that("MA order errors describe MA syntax", {
  data = data.frame(x = 1:6, y = 1:6)

  expect_error(
    get_rhs_table(
      list(y ~ ma(x)),
      data,
      mcpfamily(gaussian()),
      par_x = "x"
    ),
    "Must be ma(order) or ma(order, formula)",
    fixed = TRUE
  )
})


test_that("mcp rejects parsed MA terms until they are implemented", {
  data = data.frame(x = 1:6, y = 1:6)

  expect_error(
    mcp(list(y ~ 1 + ma(1)), data, par_x = "x", sample = FALSE),
    "ma() syntax is recognized, but moving-average terms are not implemented yet.",
    fixed = TRUE
  )
})
