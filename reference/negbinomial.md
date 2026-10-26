# Negative Binomial for mcp

Parameterized as `mu` (the conditional mean) and `shape` (the same
quantity as `size` in
[`rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html)). Thus
`Var(y) = mu + mu^2 / shape`, which approaches the Poisson variance as
`shape` approaches infinity. [Read
more](https://lindeloev.github.io/mcp/articles/poisson.html#negative-binomial-extension).

## Usage

``` r
negbinomial(link = "log", link_shape = "log")
```

## Arguments

- link:

  Link function for `mu`.

- link_shape:

  Link function for `shape` when it is modeled using `shape()`.

## Details

`shape` is required by the likelihood, so it includes an intercept in
segment 1. Later segments that do not include `shape()` are joined. For
example, `y ~ 1 + x + shape(1 + x)` models both the mean and shape.
Regression coefficients for both dpars are on their link scales, e.g.,
`shape_Intercept_1` and `shape_x_1`.

Without `shape()`, `shape_1` is the shape itself with a
`dinvgamma(0.4, 0.3)` prior, as in `brms`, much like the implicit
`sigma_1` in [`gaussian()`](https://rdrr.io/r/stats/family.html).

## Examples

``` r
# Fit an overdispersed count model with the default log links
data = data.frame(time = 1:6, count = c(1, 2, 8, 3, 12, 5))
fit = mcp(list(count ~ 1), data, family = negbinomial(), par_x = "time", sample = FALSE)
mcp_pars(fit)$parameter  # Show the mean and shape parameters
#> [1] "Intercept_1" "shape_1"    
```
