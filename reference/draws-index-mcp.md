# Index `mcpfit` objects

Index variables, iterations, chains, and draws.

## Usage

``` r
# S3 method for class 'mcpfit'
niterations(object, prior = FALSE, ...)

# S3 method for class 'mcpfit'
nchains(object, prior = FALSE, ...)
```

## Arguments

- object:

  An `mcpfit` object.

- prior:

  TRUE/FALSE. Plot using prior samples? Useful for
  `mcp(..., sample = "both")`

- ...:

  Currently ignored.

## Functions

- `niterations(mcpfit)`: Total number of iterations of an `mcpfit`
  object.

- `nchains(mcpfit)`: Number of chains of an `mcpfit` object.

## Examples

``` r
niterations(demo_fit)
#> [1] 3000
nchains(demo_fit, prior = TRUE)
#> [1] 3
```
