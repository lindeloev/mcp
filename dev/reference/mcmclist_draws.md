# Internal function to get draws.

Returns posterior draws, if available. If not, then prior draws. If not,
then throw an informative error. This is useful for summary and
plotting, that works on both.

## Usage

``` r
mcmclist_draws(fit, prior = FALSE, message = TRUE, error = TRUE)
```

## Arguments

- fit:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object

- prior:

  Logical. Summarise prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)?

- message:

  TRUE: gives a message if returning prior draws. FALSE = no message

- error:

  TRUE: err if there are no draws. FALSE: return NULL
