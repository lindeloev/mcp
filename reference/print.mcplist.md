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
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object.

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
#>  $ response ~ 1 ~ 0 + time
#>  $ response ~ 1 ~ 1 + time
print(demo_fit$prior)
#> List of 7
#>  $ cp_1       :"dirichlet(1)"
#>  $ cp_2       :"dirichlet(1)"
#>  $ Intercept_1:"dt(14.2, 6, 3)"
#>  $ time_2     :"dt(0, 0.06053105, 3)"
#>  $ Intercept_3:"dt(14.2, 6, 3)"
#>  $ time_3     :"dt(0, 0.06053105, 3)"
#>  $ sigma_1    :"dt(0, 6, 3) T(0.001, )"
```
