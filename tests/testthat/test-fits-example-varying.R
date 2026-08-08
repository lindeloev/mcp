if (Sys.getenv("MCP_TEST_LEVEL") != "release") testthat::skip("Time-consuming example fit recovery tests are only run when MCP_TEST_LEVEL='release'.")
testthat::skip_if_not_installed("vdiffr")
testthat::test_that("mcp example varying", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  test_mcp_example("varying")
})
