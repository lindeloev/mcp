# Add Log-Likelihood to an mcpfit Object.

Add Log-Likelihood to an mcpfit Object.

## Usage

``` r
add_loglik(
  x,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated()
)
```

## Arguments

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object.

- varying:

  One of:

  - `TRUE` All varying effects (`fit$pars$varying`).

  - `FALSE` No varying effects ([`c()`](https://rdrr.io/r/base/c.html)).

  - `"cp"` or `"predictor"`: All varying effects belonging to that part
    of the model.

  - Character vector: Only include specified varying parameters - see
    `fit$pars$varying`.

- arma:

  Whether to include AR and MA effects.

  - `TRUE` Compute the GARMA residual recurrence. Requires the response
    variable in `newdata`.

  - `FALSE` Disregard AR and MA effects. For `family = gaussian()`,
    [`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
    uses only `sigma` for residuals.

- ndraws:

  Integer or `NULL`. Number of posterior draws used for the
  log-likelihood or information criterion. `NULL` uses all draws.

- nsamples:

  Deprecated. Use `ndraws` instead.

## Value

An `mcpfit` object with `fit$loglik` filled as an (Nchains \* Ndraws) by
N-observed-responses data matrix with chain number as rownames. Rows
with missing responses are excluded from the likelihood columns.

## See also

loo.mcpfit waic.mcpfit

## Examples

``` r
# \donttest{
demo_fit = add_loglik(demo_fit)
# }
```
