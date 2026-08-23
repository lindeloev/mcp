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

## Examples

``` r
# Fit an overdispersed count model with the default log links
data = data.frame(time = 1:6, count = c(1, 2, 8, 3, 12, 5))
fit = mcp(list(count ~ 1), data, family = negbinomial(), par_x = "time", sample = FALSE)
mcp_pars(fit)  # Show the mean and shape parameters
#> # A tibble: 2 × 9
#>   name        part     scope role  segment dpar  order group_col population_name
#>   <chr>       <chr>    <chr> <chr>   <int> <chr> <int> <chr>     <chr>          
#> 1 Intercept_1 predict… popu… fixe…       1 mu       NA NA        NA             
#> 2 shape_1     predict… popu… dpar…       1 shape    NA NA        NA             
```
