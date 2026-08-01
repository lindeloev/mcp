# Contributing to mcp

Thank you for your interest in contributing to `mcp`! This document outlines guidelines and developer workflows for running tests and contributing code.

## Running Tests

`mcp` uses `testthat` (edition 3) for unit and integration testing.

To balance fast feedback during local development and CRAN checks with exhaustive statistical validation, the test suite supports two execution levels set via the `MCP_TEST_LEVEL` environment variable. JAGS sampling is fast even with the full model matrix, so both levels always run the full suite of formulas/models with brief MCMC sampling - only long-running statistical recovery/visual checks are reserved for `"release"`:

| Command | Typical Runtime | What it Runs |
|  :--- | :--- | :--- |
| *(default, i.e. `MCP_TEST_LEVEL` unset)* | ~20-40s | Formula parsing, JAGS code generation, and brief 18-iteration MCMC sampling across all ~150 models, plus S3 method validation (`loo`/`waic`/`summary`/`plot_pars`/`hypothesis`) on 5 representative archetypes. This is what CRAN and AI agents should run. |
| `Sys.setenv(MCP_TEST_LEVEL = "release")` | Extended | Everything above, plus full parameter recovery tests (`test-fits-*.R`), example fit/plot validation (`test-fits-examples.R`), and S3 method tests across all ~150 formula variants. Only needed before a release/checkpoint. |

By default, `testthat` uses 2 parallel worker processes. For maximum test execution speed on multi-core machines (e.g., 11 workers), set `TESTTHAT_CPUS`:

```R
Sys.setenv(TESTTHAT_CPUS = 11)
devtools::test()
```

Or pass `cpus` directly to `devtools::test()`:
```R
devtools::test(cpus = 11)
```

### A Note on `test-fits-examples.R` and `vdiffr` Snapshots

`test-fits-examples.R` (which fits every `mcp_example()` and compares plots
against `vdiffr` reference snapshots in `tests/testthat/_snaps/fits-examples/`)
only actually runs at `MCP_TEST_LEVEL = "release"`; it's `skip()`-ped otherwise.
testthat deletes any snapshot file it doesn't see referenced during a run
(`testthat:::SnapshotReporter$end_reporter()`), *unless* it detects it's
running on CI (`Sys.getenv("CI") == "true"`) - so without protection, every
default-level run would silently delete the reference SVGs from disk.
`tests/testthat.R` and `tests/testthat/helper-runs.R` both set `CI = "true"`
automatically whenever `MCP_TEST_LEVEL != "release"` to prevent this. If you
call `devtools::test()`/`testthat::test_local()` directly from an R console
(bypassing `tests/testthat.R`), set it yourself for the same protection:

```R
Sys.setenv(CI = "true")
devtools::test()
```

### Rendering and Reviewing Reference Plots

`mcp` includes visual testing tools for reviewing generated plots and managing snapshot regression tests:

1. **Exporting Example Plots to PDF**:
   During `"release"` mode testing, you can export all rendered plots from `mcp_example()` into a multi-page PDF for review:
   ```R
   Sys.setenv(MCP_TEST_LEVEL = "release", MCP_MAKE_PLOT_TEST_FILE = "plots_review.pdf")
   devtools::test()
   ```

2. **Updating `vdiffr` Visual Snapshots**:
   Visual snapshots for example plots (`demo`, `group`, `variance`) are checked using [`vdiffr`](https://vdiffr.r-lib.org/). When plot aesthetics or layout changes are introduced intentionally, you can review and update the reference SVG snapshots stored in `tests/testthat/_snaps/`:
   
   - **Interactive review (Shiny app)**:
     ```R
     vdiffr::manage_cases()
     ```
     This launches a Shiny UI where you can visually inspect side-by-side plot diffs and individually accept or reject new reference images.
   
   - **Accept all updated snapshots**:
     ```R
     testthat::snapshot_accept()
     ```

### Environment Variables Summary

* `MCP_TEST_LEVEL`: Unset (default) or `"release"`.
* `TESTTHAT_CPUS`: Number of parallel test worker processes to spawn (e.g., `11`).
* `MCP_MAKE_PLOT_TEST_FILE`: Path to export rendered example plots into a multi-page PDF during release testing.
