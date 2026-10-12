# Print mcplist

Shows a list in a more condensed format using `str(list)`.

## Usage

``` r
# S3 method for class 'mcplist'
print(x, ...)
```

## Arguments

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object, or an `mcpsummary` returned by
  [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md).

- ...:

  Must be empty. Reserved for future use.

## Value

`x`, invisibly.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
print(demo_fit$model)
#> List of 3
#>  $ response ~ 1
#>  $ ~0 + time
#>  $ ~1 + time
print(demo_fit$prior)
#> List of 7
#>  $ cp_1       :"dirichlet(1)"
#>  $ cp_2       :"dirichlet(1)"
#>  $ Intercept_1:"dnorm(13.825120048974153, 17.382500430759464)"
#>  $ time_2     :"dnorm(0, 0.62123709509255609)"
#>  $ Intercept_3:"dnorm(13.825120048974153, 17.382500430759464)"
#>  $ time_3     :"dnorm(0, 0.62123709509255609)"
#>  $ sigma_1    :"dt(0, 5.7, 3) T(0, )"
```
