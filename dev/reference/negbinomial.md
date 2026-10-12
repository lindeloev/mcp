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

  Link function for `shape` when it is modeled using `shape()`.

## Details

`shape` is required by the likelihood, so it includes an intercept in
segment 1. Later segments that do not include `shape()` are joined. For
example, `y ~ 1 + x + shape(1 + x)` models both the mean and shape.
Regression coefficients for both dpars are on their link scales.

Without `shape()`, `shape_1` is the shape itself with an
`inv_gamma(0.4, 0.3)` prior, as in `brms`, much like the implicit
`sigma_1` in [`gaussian()`](https://rdrr.io/r/stats/family.html).

## Examples

``` r
# Fit an overdispersed count model with the default log links
data = data.frame(time = 1:6, count = c(1, 2, 8, 3, 12, 5))
fit = mcp(list(count ~ 1), data, family = negbinomial(), par_x = "time", sample = FALSE)
mcp_pars(fit)  # Show the mean and shape parameters
#> # A tibble: 2 × 9
#>   parameter part  scope role  segment dpar  order group_col population_parameter
#>   <chr>     <chr> <chr> <chr>   <int> <chr> <int> <chr>     <chr>               
#> 1 Intercep… pred… popu… fixe…       1 mu       NA NA        NA                  
#> 2 shape_1   pred… popu… dpar…       1 shape    NA NA        NA                  
```
