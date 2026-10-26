# Deprecated functions in mcp

These functions were exported in mcp v0.3.4 and are deprecated in
v0.4.0. They still work but warn and point to their replacements.

## Usage

``` r
phi(q)
```

## Arguments

- q:

  Vector of quantiles.

## Value

A numeric vector of probabilities.

## Details

- `phi()`: The inverse probit link. Use
  [`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html).
