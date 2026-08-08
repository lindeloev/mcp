# Index `mcpfit` objects

Index variables, iterations, chains, and draws.

## Usage

``` r
# S3 method for class 'mcpfit'
niterations(object, ...)

# S3 method for class 'mcpfit'
nchains(object, ...)

# S3 method for class 'mcpfit'
ndraws(object, ...)

ndraws(x)

nchains(x)

niterations(x)
```

## Arguments

- object:

  An `mcpfit` object.

- ...:

  Currently ignored.

## Functions

- `niterations(mcpfit)`: Number of iterations per chain of an `mcpfit`
  object.

- `nchains(mcpfit)`: Number of chains of an `mcpfit` object.

## Examples

``` r
niterations(demo_fit)
#> [1] 1000
nchains(as_draws(demo_fit, prior = TRUE))
#> [1] 3
```
