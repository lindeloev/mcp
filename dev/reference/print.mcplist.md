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
#> <environment: 0x55fd5c0d7c58>
#>  $ response ~ 1 ~ 0 + time
#> <environment: 0x55fd5c0d7c58>
#>  $ response ~ 1 ~ 1 + time
#> <environment: 0x55fd5c0d7c58>
print(demo_fit$prior)
#> List of 7
#>  $ cp_1       :"dt(3.189323, 47.38506, 1) T(3.189323, 97.95945)"
#>  $ cp_2       :"dt(3.189323, 47.38506, 1) T(cp_1, 97.95945)"
#>  $ Intercept_1:"dt(10.4, 12.8, 3)"
#>  $ time_2     :"dt(0, 0.1350637, 3)"
#>  $ Intercept_3:"dt(10.4, 12.8, 3)"
#>  $ time_3     :"dt(0, 0.1350637, 3)"
#>  $ sigma_1    :"dt(0, 12.8, 3) T(0, )"
```
