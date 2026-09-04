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

  Optional parameter role(s), such as `"fixed_effect"`, `"dpar_effect"`,
  `"arma"`, `"group_sd"`, or `"group_deviation"`.

## Value

A data frame with one row per parameter definition. `name` is the model
parameter name; `dpar` identifies its distributional parameter or
component (`"cp"`, `"ar"`, or `"ma"`); `order` gives an AR/MA lag;
`group_col` and `population_name` describe group-level effects.

## Details

Return the canonical parameter definitions for an `mcpfit` object. This
works before sampling and is the stable way to discover model
parameters.

## Examples

``` r
# Show every parameter in the model
mcp_pars(demo_fit)
#> # A tibble: 7 × 9
#>   name        part     scope role  segment dpar  order group_col population_name
#>   <chr>       <chr>    <chr> <chr>   <int> <chr> <int> <chr>     <chr>          
#> 1 cp_1        cp       popu… chan…       2 cp       NA NA        NA             
#> 2 cp_2        cp       popu… chan…       3 cp       NA NA        NA             
#> 3 Intercept_1 predict… popu… fixe…       1 mu       NA NA        NA             
#> 4 time_2      predict… popu… fixe…       2 mu       NA NA        NA             
#> 5 Intercept_3 predict… popu… fixe…       3 mu       NA NA        NA             
#> 6 time_3      predict… popu… fixe…       3 mu       NA NA        NA             
#> 7 sigma_1     predict… popu… dpar…       1 sigma    NA NA        NA             

# Select population-level coefficients
mcp_pars(demo_fit, scope = "population", role = "fixed_effect")
#> # A tibble: 4 × 9
#>   name        part     scope role  segment dpar  order group_col population_name
#>   <chr>       <chr>    <chr> <chr>   <int> <chr> <int> <chr>     <chr>          
#> 1 Intercept_1 predict… popu… fixe…       1 mu       NA NA        NA             
#> 2 time_2      predict… popu… fixe…       2 mu       NA NA        NA             
#> 3 Intercept_3 predict… popu… fixe…       3 mu       NA NA        NA             
#> 4 time_3      predict… popu… fixe…       3 mu       NA NA        NA             
```
