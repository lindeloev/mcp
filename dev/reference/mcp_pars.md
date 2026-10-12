# Model parameters

**\[experimental\]**

## Usage

``` r
mcp_pars(fit, scope = NULL, role = NULL)
```

## Arguments

- fit:

  An `mcpfit` object.

- scope:

  Optional parameter scope(s): `"population"` or `"group"`.

- role:

  Optional parameter role(s): `"change_point"`, `"coefficient"`
  (regression coefficients of all distributional parameters, e.g., `mu`,
  `sigma`, and `ar1`), or `"group_sd"`. Combine with `scope` to select
  population-level parameters or their group-level deviations, e.g.,
  `scope = "group", role = "change_point"` for group-level change
  points.

## Value

A tibble with one row per parameter definition. `parameter` is the model
parameter name; `dpar` is its distributional parameter (e.g., `"mu"`,
`"sigma"`, or `"ar1"`; see
[mcp-formula](https://lindeloev.github.io/mcp/dev/reference/mcp-formula.md)),
as accepted by the `dpar` argument of
[`fitted.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
[`fixef.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
etc., and `NA` for change points; `group_col` and `population_parameter`
describe group-level effects.

## Details

Return the canonical parameter definitions for an `mcpfit` object. This
works before sampling and is the stable way to discover model
parameters.

## Examples

``` r
# Show every parameter in the model
mcp_pars(demo_fit)
#> # A tibble: 7 × 8
#>   parameter   part      scope role  segment dpar  group_col population_parameter
#>   <chr>       <chr>     <chr> <chr>   <int> <chr> <chr>     <chr>               
#> 1 cp_1        cp        popu… chan…       2 NA    NA        NA                  
#> 2 cp_2        cp        popu… chan…       3 NA    NA        NA                  
#> 3 Intercept_1 predictor popu… coef…       1 mu    NA        NA                  
#> 4 time_2      predictor popu… coef…       2 mu    NA        NA                  
#> 5 Intercept_3 predictor popu… coef…       3 mu    NA        NA                  
#> 6 time_3      predictor popu… coef…       3 mu    NA        NA                  
#> 7 sigma_1     predictor popu… coef…       1 sigma NA        NA                  

# Select population-level coefficients
mcp_pars(demo_fit, scope = "population", role = "coefficient")
#> # A tibble: 5 × 8
#>   parameter   part      scope role  segment dpar  group_col population_parameter
#>   <chr>       <chr>     <chr> <chr>   <int> <chr> <chr>     <chr>               
#> 1 Intercept_1 predictor popu… coef…       1 mu    NA        NA                  
#> 2 time_2      predictor popu… coef…       2 mu    NA        NA                  
#> 3 Intercept_3 predictor popu… coef…       3 mu    NA        NA                  
#> 4 time_3      predictor popu… coef…       3 mu    NA        NA                  
#> 5 sigma_1     predictor popu… coef…       1 sigma NA        NA                  
```
