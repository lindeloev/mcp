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
fitting-data rows.

## Examples

``` r
formula(demo_fit)  # Show all segment formulas
#> List of 3
#>  $ response ~ 1
#> <environment: 0x5629ddcea9c8>
#>  $ response ~ 1 ~ 0 + time
#> <environment: 0x5629ddcea9c8>
#>  $ response ~ 1 ~ 1 + time
#> <environment: 0x5629ddcea9c8>
formula(demo_fit, segment = 2)  # Show the formula for segment 2
#> response ~ 1 ~ 0 + time
#> <environment: 0x5629ddcea9c8>
family(demo_fit)  # Show the response family and link
#> Family: gaussian
#> Links: mu = identity; sigma = identity
head(model.frame(demo_fit))  # Show the top rows of fitting data
#>    response     time
#> 1 32.842651 68.35820
#> 2 -1.160003 87.29038
#> 3 27.564248 69.01173
#> 4 10.062971 11.59361
#> 5 14.056859 19.50091
#> 6 18.292640 46.12009
nobs(demo_fit)  # Count observed response rows
#> [1] 100
```
