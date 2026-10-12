# Extract Model Information from an `mcpfit`

Standard R accessors for the model formulas, family, fitting data, and
number of observations stored in an `mcpfit`.

## Usage

``` r
# S3 method for class 'mcpfit'
family(object, ...)

# S3 method for class 'mcpfit'
nobs(object, ...)

# S3 method for class 'mcpfit'
model.frame(formula, ...)

# S3 method for class 'mcpfit'
formula(x, segment = NULL, ...)
```

## Arguments

- object, x, formula:

  An `mcpfit` object.

- ...:

  Must be empty. Reserved for future use.

- segment:

  `NULL` to return all segment formulas, or a positive integer selecting
  one segment.

## Value

[`formula()`](https://rdrr.io/r/stats/formula.html) returns the complete
list of segment formulas, or one formula when `segment` is supplied.
[`family()`](https://rdrr.io/r/stats/family.html) returns an
`mcpfamily`. [`model.frame()`](https://rdrr.io/r/stats/model.frame.html)
returns the data retained in the fit.
[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
observed response values (excluding missing responses).

## Examples

``` r
formula(demo_fit)  # Show all segment formulas
#> List of 3
#>  $ response ~ 1
#>  $ ~0 + time
#>  $ ~1 + time
formula(demo_fit, segment = 2)  # Show the formula for segment 2
#> ~0 + time
family(demo_fit)  # Show the response family and link
#> Family: gaussian
#> Links: mu = identity; sigma = identity
head(model.frame(demo_fit))  # Show the top rows of fitting data
#>    response     time
#> 1 12.235516 76.33986
#> 2  6.351706 83.51711
#> 3 28.049953 60.18529
#> 4 15.681983 74.72964
#> 5 16.213643 85.88256
#> 6 22.322824 40.05069
nobs(demo_fit)  # Count observed response rows
#> [1] 100
```
