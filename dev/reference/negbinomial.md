# Negative Binomial for mcp

Parameterized as `mu` (the conditional mean) and `shape` (the same
quantity as `size` in
[`rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html)). Thus
`Var(y) = mu + mu^2 / shape`, which approaches the Poisson variance as
`shape` approaches infinity.

## Usage

``` r
negbinomial(link = "log", link_shape = "log")
```

## Arguments

- link:

  Link function for `mu`.

- link_shape:

  Link function for `shape`.

## Details

`shape(1)` is added implicitly and is constant across segments unless a
`shape()` formula is supplied. For example, `y ~ 1 + x + shape(1 + x)`
models both the mean and shape. Regression coefficients for both dpars
are on their link scales.
