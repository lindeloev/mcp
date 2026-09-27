test_that("mcp examples simulate values near their default priors", {
  examples = c("ar", "binomial", "demo", "group_cp", "group_mu", "intercepts", "missing", "multiple", "quadratic", "sigma")
  # Quadratics confined to a short segment are far from priors scaled to sd((x - min(x))^2).
  # Pending decision in "dev/plan - priors and sampling.md".
  ignore = list(multiple = "xE2_3", quadratic = "xE2_2")
  for (name in examples) {
    fit = suppressMessages(mcp_example(name, sample = FALSE, plot = FALSE))
    simulated = attr(fit$data[[mcp_columns(fit)$response]], "simulated")
    expect_false(is.null(simulated), label = name)
    expect_prior_distance(fit, simulated, label = paste0("mcp_example(\"", name, "\")"), ignore = if (is.null(ignore[[name]])) character() else ignore[[name]])
  }
})


test_that("prior_distance() parses dt() and dnorm() priors and skips others", {
  fit = list(prior = list(
    a = "dt(10, 2, 3)",
    b = "dnorm(0, 0.5) T(0, )",
    c = "dirichlet(1)",
    d = "x_2",
    e = 5
  ))
  distances = prior_distance(fit, list(a = 16, b = 1, c = 0.3, d = 1, e = 5))
  expect_equal(distances$parameter, c("a", "b"))
  expect_equal(distances$distance, c(3, 2))
})
