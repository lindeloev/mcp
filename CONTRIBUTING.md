# Contributing to mcp

Thank you for your interest in contributing to `mcp`! This document outlines guidelines and developer workflows for running tests and contributing code.

## Submitting code
Create pull request to the `dev` branch. Before that:
* Run `Sys.setenv(MCP_TEST_LEVEL = "release"); devtools::test()`
* Go through relevant `mcp::release_questions()` - skip those only relevant for CRAN upload.
* Once pull request is created, check the workflows initiated on github that all tests pass.

AI-based bug finding is welcome and is actively used during development of `mcp`. See `dev/promot*` for prompts.

See `dev/DECISIONS.md` for decisions made during development to balance functionality/bug-free against code simplicity. It is better to have simple/readable code than to add 200 lines of code to fix a rare edge case.


## Branches
* `main` is the current CRAN version. It renders site https://lindeloev.github.io/mcp/ which is the public documentation of the CRAN version, tracked by crawlers etc.
* `dev` is current dev release. It should generally be bug-free but the API may change until merged into `main` (on new CRAN release). `dev` renders site to https://lindeloev.github.io/mcp/dev/ which is not tracked by crawlers etc.
* Other branches are feature branches which should be merged into `dev` if intended for eventual release.


## About test levels
`Sys.setenv(MCP_TEST_LEVEL = "release"); devtools::test()` is intended for development testing. 
 * Runs a few heavier tests involving full model fits (see `tests/testthat/fits-*.R` scripts)
 * Renders and compares `mcp_example()` plots to reference (using `vdiffr`), etc. For speedup, also set e.g. `Sys.setenv(NCPUs = 11)`. If plots mismatch, use `testthat::snapshot_review()` to review.

`Sys.setenv(CI = "true"); devtools::test()` (without MCP_TEST_LEVEL set) is the defaul CRAN/github test setting where short test times are required. It includes essential and fast tests.
