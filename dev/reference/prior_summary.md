# Summarise priors used by an mcp model

Shows the effective, resolved priors on the familiar SD/scale
parameterization rather than JAGS precision. Use `verbose = TRUE` to
also see the symbolic rule, its description, source, and kind.

## Usage

``` r
prior_summary(fit, verbose = FALSE)
```

## Arguments

- fit:

  An `mcpfit` object.

- verbose:

  Logical. Include rule, description, source, and kind.

## Value

A tibble with one row per model parameter, ordered and labeled the same
way as
[`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md):
change points first, then `mu`, then the other distributional
parameters, then `ar`/`ma` components - each with `segment` and `dpar`
columns.
