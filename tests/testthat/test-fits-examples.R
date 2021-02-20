testthat::skip_if(is.null(options("test_mcp_fits")[[1]]),
                  "This time-consuming test is only run locally before release.")

# Runs through all examples in mcp::mcp_example() and verifies that the parameters are recovered
# Do not run "trigonometric"; it is know to have unstable recovery
examples = c("ar", "binomial", "demo", "intercepts", "multiple", "quadratic", "variance", "varying")
for (example in examples) {
  message("Now running example mcp_example('", example, "')")
  fit = mcp_example(example, sample = TRUE)$fit
  summaries = rbind(fixef(fit, width = 0.97), ranef(fit, width = 0.97)) %>%
    dplyr::filter(is.na(sim) == FALSE)

  new_lower = summaries$lower - 0.1*(summaries$mean - summaries$lower)
  new_upper = summaries$upper - 0.1*(summaries$mean - summaries$upper)
  correctly_estimated = summaries$match == "OK" | (summaries$sim > new_lower & summaries$sim < new_upper)
  testthat::expect_true(all(correctly_estimated))
}
