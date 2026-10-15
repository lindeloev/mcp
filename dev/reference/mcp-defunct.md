# Defunct functions in mcp

These functions were exported in mcp v0.3.4 and removed in v0.4.0. They
now throw an error pointing to their replacements.

## Usage

``` r
criterion(...)

get_segment_table(...)
```

## Arguments

- ...:

  Ignored.

## Details

- `criterion()`: Use
  [`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
  or
  [`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md).

- `get_segment_table()`: Use `fit = mcp(..., sample = FALSE)` followed
  by
  [`mcp_pars()`](https://lindeloev.github.io/mcp/dev/reference/mcp_pars.md)
  or
  [`mcp_columns()`](https://lindeloev.github.io/mcp/dev/reference/mcp_columns.md).
