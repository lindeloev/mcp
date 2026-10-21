# Index `mcpfit` objects

Index variables, iterations, chains, and draws.

## Usage

``` r
# S3 method for class 'mcpfit'
variables(x, ...)

# S3 method for class 'mcpfit'
niterations(x, ...)

# S3 method for class 'mcpfit'
nchains(x, ...)

# S3 method for class 'mcpfit'
ndraws(x, ...)

variables(x, ...)

ndraws(x)

nchains(x)

niterations(x)
```

## Arguments

- x:

  An `mcpfit` object or a posterior draws object.

- ...:

  Must be empty. Reserved for future use.

## Value

`variables()` returns the names of the draws' variables, e.g., `"cp_1"`
and `"cp_1_id[John]"`. The other methods return an integer count of
iterations, chains, or draws.

## Details

These methods require posterior draws. To index prior draws, use e.g.
`ndraws(as_draws(fit, prior = TRUE))`.

## Functions

- `variables(mcpfit)`: Variable names of an `mcpfit` object, as in
  `as_draws(fit)`.

- `niterations(mcpfit)`: Number of iterations per chain of an `mcpfit`
  object.

- `nchains(mcpfit)`: Number of chains of an `mcpfit` object.

## Examples

``` r
variables(demo_fit)
#> [1] "cp_1"        "cp_2"        "Intercept_1" "time_2"      "Intercept_3" "time_3"      "sigma_1"    
niterations(demo_fit)
#> [1] 500
nchains(as_draws(demo_fit, prior = TRUE))
#> [1] 2
```
