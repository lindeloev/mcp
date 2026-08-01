# Internal function to get samples.

Returns posterior samples, if available. If not, then prior samples. If
not, then throw an informative error. This is useful for summary and
plotting, that works on both.

## Usage

``` r
mcmclist_samples(fit, prior = FALSE, message = TRUE, error = TRUE)
```

## Arguments

- fit:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object

- prior:

  TRUE/FALSE. Summarise prior instead of posterior?

- message:

  TRUE: gives a message if returning prior samples. FALSE = no message

- error:

  TRUE: err if there are no samples. FALSE: return NULL
