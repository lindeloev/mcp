# COMMENT THIS LINE TO ENABLE EXTENSIVE TESTS. THEN REVERT. DO NOT COMMIT.
# It's a slow test so we do it rarely.
testthat::skip("This time-consuming test is only run locally before release.")

# Runs through all examples in mcp::mcp_example() and verifies that the parameters are recovered
examples = c("ar", "binomial", "demo", "group", "intercepts", "multiple", "quadratic", "variance", "varying")
for (example in examples) {
  message("Now running example mcp_example('", example, "')")
  fit = mcp_example(example, plot = FALSE)
  test_matches_simulated(fit)
}
