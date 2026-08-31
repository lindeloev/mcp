# Index `mcpfit` objects

Index variables, iterations, chains, and draws.

## Usage

``` r
# S3 method for class 'mcpfit'
niterations(x, ...)

# S3 method for class 'mcpfit'
nchains(x, ...)

# S3 method for class 'mcpfit'
ndraws(x, ...)

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

An integer count of iterations, chains, or draws.

## Functions

- `niterations(mcpfit)`: Number of iterations per chain of an `mcpfit`
  object.

- `nchains(mcpfit)`: Number of chains of an `mcpfit` object.

## Examples

``` r
niterations(demo_fit)
#> [1] 1000
nchains(as_draws(demo_fit, prior = TRUE))
#> [1] 2
```
