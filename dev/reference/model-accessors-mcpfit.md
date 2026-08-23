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

  Currently unused.

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
fitting-data rows.

## Examples

``` r
formula(demo_fit)  # Show all segment formulas
#> List of 3
#>  $ response ~ 1
#> <environment: 0x55ee2832faa0>
#>  $ response ~ 1 ~ 0 + time
#> <environment: 0x55ee2832faa0>
#>  $ response ~ 1 ~ 1 + time
#> <environment: 0x55ee2832faa0>
formula(demo_fit, segment = 2)  # Show the formula for segment 2
#> response ~ 1 ~ 0 + time
#> <environment: 0x55ee2832faa0>
family(demo_fit)  # Show the response family and link
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
head(model.frame(demo_fit))  # Show the top rows of fitting data
#>     response     time
#> 1 -3.0084198 91.48060
#> 2 -7.8768640 93.70754
#> 3 16.3029101 28.61395
#> 4 -0.0373553 83.04476
#> 5 27.4463185 64.17455
#> 6 22.0610004 51.90959
nobs(demo_fit)  # Count observed response rows
#> [1] 100
```
