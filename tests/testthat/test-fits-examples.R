testthat::skip("This time-consuming test is only run locally before release.")

# Runs through all examples in mcp::mcp_example() and verifies that the parameters are recovered
# Do not run "trigonometric"; it is know to have unstable recovery
examples = c("ar", "binomial", "demo", "intercepts", "multiple", "quadratic", "variance", "varying")
for (example in examples) {
  message("Now running example mcp_example('", example, "')")
  fit = mcp_example(example)
  test_matches_simulated(fit)
}
