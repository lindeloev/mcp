# Extract MCMC Draws from `mcpfit` Objects

Extract posterior or prior draws using posterior, tidybayes, or coda S3
generics.

## Usage

``` r
# S3 method for class 'mcpfit'
as_draws(x, prior = FALSE, ...)
```

## Arguments

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object.

- prior:

  Logical. Extract prior draws (`TRUE`) instead of posterior draws
  (`FALSE`)?

- ...:

  Passed to posterior or tidybayes format conversion functions.

## Value

A posterior `draws` object or a coda `mcmc.list` object.
