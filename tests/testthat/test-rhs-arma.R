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


test_that("GARMA boundaries default to 0.1 and can vary by segment", {
  data = data.frame(x = 1:6, y = 1:6)
  rhs_table = get_rhs_table(
    list(
      y ~ ar(1),
      ~ ar(1, 1 + x, boundary = 0.2)
    ),
    data,
    mcpfamily(gaussian()),
    par_x = "x"
  )

  boundaries = rhs_table %>%
    dplyr::filter(.data$dpar == "ar") %>%
    dplyr::distinct(.data$segment, .data$boundary)
  expect_equal(boundaries$boundary, c(0.1, 0.2))

  expect_error(
    get_rhs_table(list(y ~ ar(1, boundary = 0)), data, mcpfamily(gaussian()), "x"),
    "`boundary` in ar() must be one number between 0 and 1.",
    fixed = TRUE
  )
})


test_that("AR and MA share one boundary within a segment", {
  data = data.frame(x = 1:6, y = 1:6)
  family = mcpfamily(gaussian())

  rhs_table = get_rhs_table(
    list(y ~ ar(1, boundary = 0.2) + ma(1)),
    data,
    family,
    par_x = "x"
  )
  expect_equal(unique(rhs_table$boundary[rhs_table$dpar == "ar"]), 0.2)
  expect_equal(unique(rhs_table$boundary[rhs_table$dpar == "ma"]), 0.2)

  expect_error(
    get_rhs_table(
      list(y ~ ar(1, boundary = 0.2) + ma(1, boundary = 0.3)),
      data,
      family,
      par_x = "x"
    ),
    "ar() and ma() must use the same `boundary` within a segment.",
    fixed = TRUE
  )
})


test_that("mcp builds MA parameters and JAGS code", {
  data = data.frame(x = 1:6, y = 1:6)
  fit = mcp(list(y ~ 1 + ma(1)), data, par_x = "x", sample = FALSE)

  expect_true("ma1_1" %in% fit$pars$arma)
  expect_match(fit$jags_code, "ma1_\\[i_\\] \\* resid_ma_\\[i_ - 1\\]")
})
