library(testthat)

# Avoid testthat deleting the vdiffr reference snapshots in
# tests/testthat/_snaps/fits-examples/ whenever the example-fit test that
# owns them isn't actually run (i.e. at any level other than "release"). See
# CONTRIBUTING.md for details.
if (Sys.getenv("MCP_TEST_LEVEL") != "release" && !isTRUE(as.logical(Sys.getenv("CI", "false")))) {
  Sys.setenv(CI = "true")
}

test_check("mcp")
