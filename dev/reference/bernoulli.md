# Bernoulli Family for mcp

Bernoulli Family for mcp

## Usage

``` r
bernoulli(link = "logit")
```

## Arguments

- link:

  Link function.

## Examples

``` r
# Fit a binary-response model with a probit link
data = data.frame(time = 1:6, y = c(0, 0, 0, 1, 1, 1))
fit = mcp(list(y ~ 1), data, family = bernoulli(link = "probit"), par_x = "time", sample = FALSE)
mcp_pars(fit)  # Show the parameters of the fitted Bernoulli model
#> # A tibble: 1 × 8
#>   parameter   part      scope role  segment dpar  group_col population_parameter
#>   <chr>       <chr>     <chr> <chr>   <int> <chr> <chr>     <chr>               
#> 1 Intercept_1 predictor popu… coef…       1 mu    NA        NA                  
```
