# Model data columns

**\[experimental\]**

## Usage

``` r
mcp_columns(fit)
```

## Arguments

- fit:

  An `mcpfit` object.

## Value

A named list with `par_x`, `response`, and `series`, plus any
family-defined response-auxiliary column roles.

## Details

Return the resolved data-column roles for an `mcpfit` object. In
particular, `par_x` is the change-point predictor chosen automatically
or supplied to
[`mcp()`](https://lindeloev.github.io/mcp/reference/mcp.md).
Family-defined response auxiliaries, such as `trials` and `weights`, are
included when relevant.

## Examples

``` r
# Show the predictor, response, and auxiliary-column roles
mcp_columns(demo_fit)
#> $par_x
#> [1] "time"
#> 
#> $response
#> [1] "response"
#> 
#> $series
#> NULL
#> 
#> $weights
#> NULL
#> 
```
