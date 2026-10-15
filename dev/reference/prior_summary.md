# Summarise priors used by an mcp model

**\[experimental\]**

## Usage

``` r
prior_summary(object, ...)

# S3 method for class 'mcpfit'
prior_summary(object, verbose = FALSE, ...)
```

## Arguments

- object:

  An `mcpfit` object.

- ...:

  Must be empty. Reserved for future use.

- verbose:

  Logical. Include rule, description, source, and kind.

## Value

A tibble with one row per model parameter, ordered and labeled the same
way as
[`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md):
change points first, then `mu`, then the other distributional
parameters, then AR/MA lags - each with `segment` and `dpar` columns.

## Details

Shows the specified prior distributions, resolved to numbers and on the
familiar SD/scale parameterization rather than JAGS precision. `bounds`
are the lower and upper bounds of each prior, i.e., the arguments of
[`dunif()`](https://rdrr.io/r/stats/Uniform.html) or
[`T()`](https://rdrr.io/r/base/logical.html), or the support of the
distribution. Change points are also "ordered": their specified priors
are restricted to ordered change points. Use
`plot_pars(fit, prior = TRUE)` to see the marginal priors after this
restriction. `verbose = TRUE` also displays the underlying rule,
description, source, and kind. Numbers are rounded to 5 significant
digits for display. `fit$prior` holds the exact values used in the
model. Parameter names within a prior describe its dependence on other
parameters, e.g., the default prior on `Intercept_1` given `x_1`.

`prior_summary()` shares its generic with rstantools (used by `brms` and
`rstanarm`), so it works for all of these models regardless of package
load order.

## See also

[mcp-priors](https://lindeloev.github.io/mcp/dev/reference/mcp-priors.md)
for how to specify priors, and [Default
priors](https://lindeloev.github.io/mcp/articles/default_priors.html)
for the rules behind them.

## Examples

``` r
prior_summary(demo_fit)  # Show the resolved priors and their bounds
#> # A tibble: 7 × 5
#>   parameter   segment dpar  prior                                        bounds 
#>   <chr>         <int> <chr> <chr>                                        <chr>  
#> 1 cp_1              2 NA    uniform(min = 0.61863, max = 99.741)         [0.618…
#> 2 cp_2              3 NA    uniform(min = 0.61863, max = 99.741)         [0.618…
#> 3 Intercept_1       1 mu    normal(mean = 13.825, sd = 17.383)           none   
#> 4 time_2            2 mu    normal(mean = 0, sd = 0.62124)               none   
#> 5 Intercept_3       3 mu    normal(mean = 13.825, sd = 17.383)           none   
#> 6 time_3            3 mu    normal(mean = 0, sd = 0.62124)               none   
#> 7 sigma_1           1 sigma student_t(df = 3, location = 0, scale = 5.7) [0, In…
prior_summary(demo_fit, verbose = TRUE)  # Include their rules and sources
#> # A tibble: 7 × 9
#>   parameter   segment dpar  prior          bounds rule  description source kind 
#>   <chr>         <int> <chr> <chr>          <chr>  <chr> <chr>       <chr>  <chr>
#> 1 cp_1              2 NA    uniform(min =… [0.61… unif… Flat Diric… defau… dist…
#> 2 cp_2              3 NA    uniform(min =… [0.61… unif… Flat Diric… defau… dist…
#> 3 Intercept_1       1 mu    normal(mean =… none   norm… Mean inter… defau… dist…
#> 4 time_2            2 mu    normal(mean =… none   norm… Autoscaled… defau… dist…
#> 5 Intercept_3       3 mu    normal(mean =… none   norm… Mean inter… defau… dist…
#> 6 time_3            3 mu    normal(mean =… none   norm… Autoscaled… defau… dist…
#> 7 sigma_1           1 sigma student_t(df … [0, I… stud… Positive r… defau… dist…
```
