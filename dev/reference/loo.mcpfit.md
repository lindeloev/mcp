# Information Criteria for Model Comparison

Compare models using
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html) and
[`loo_model_weights`](https://mc-stan.org/loo/reference/loo_model_weights.html).
more in [`loo`](https://mc-stan.org/loo/reference/loo.html).

## Usage

``` r
# S3 method for class 'mcpfit'
loo(
  x,
  ...,
  by_row = FALSE,
  pointwise = lifecycle::deprecated(),
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated()
)

# S3 method for class 'mcpfit'
waic(
  x,
  ...,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated()
)
```

## Arguments

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object.

- ...:

  Currently ignored

- by_row:

  `TRUE` calls
  [`loo.function`](https://mc-stan.org/loo/reference/loo.html) to
  compute log-likelihood contributions row-by-row
  (observation-by-observation), which is slower but more memory
  efficient. `FALSE` (default) computes the full log-likelihood matrix
  at once. Note that both modes calculate pointwise (observation-level)
  PSIS-LOO cross-validation.

- pointwise:

  Deprecated alias for `by_row`.

- varying:

  Group-level effects. One of:

  - `TRUE` All group-level deviations.

  - `FALSE` No group-level deviations
    ([`c()`](https://rdrr.io/r/base/c.html)).

  - `"cp"` or `"predictor"`: All group-level deviations belonging to
    that part of the model.

  - Character vector: Only include specified group-level parameters.

- arma:

  Whether to include AR and MA effects.

  - `TRUE` Compute the GARMA residual recurrence. Requires the response
    variable in `newdata`.

  - `FALSE` Disregard AR and MA effects. For `family = gaussian()`,
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    uses only `sigma` for residuals. For posterior evaluation of the
    original data, retained JAGS imputations supply missing GARMA
    histories. In models with group-level effects, this currently
    requires including all such effects (`varying = TRUE`).

- ndraws:

  Integer or `NULL`. Number of posterior draws used for the
  log-likelihood or information criterion. `NULL` uses all draws.

- nsamples:

  Deprecated. Use `ndraws` instead.

## Value

a `loo` or `psis_loo` object.

## Details

Observationwise PSIS-LOO and WAIC are problematic for AR/MA models
because both treat individual conditional likelihood terms as validation
units. In PSIS-LOO, a held-out response also remains in the conditioning
history of later terms. Prefer leave-future-out or blocked
cross-validation, which are not currently implemented in mcp. When a
missing response enters a later observed GARMA history,
[`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
`loo()`, and `waic()` are unavailable with `arma = TRUE`: the
observed-data likelihood requires integrating over that missing history,
which mcp does not currently implement. Setting `arma = FALSE` evaluates
the model without its AR/MA contribution rather than repairing the
criterion for the fitted serial model.

## Functions

- `waic(mcpfit)`: Computes WAIC on mcpfit objects

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# \donttest{
# Define two models and sample them
# future::plan(future::multisession, workers = 3)  # Uncomment for parallel sampling
data = mcp_example_data("intercepts")  # Get some simulated data.
model1 = list(y ~ 1 + x, ~ 1)
model2 = list(y ~ 1 + x)  # Without a change point
fit1 = mcp(model1, data)
#> Warning: Some parameters may not have converged well:
#>   * rhat > 1.01 or ess_bulk < 400 or ess_tail < 400: cp_1
#> Inspect `summary(fit)` and `plot_pars(fit)`, and consider increasing `iter`/`warmup` or simplifying the model before trusting these results.
fit2 = mcp(model2, data)

# Compute LOO for each and compare (works for waic(fit) too)
loo1 = loo(fit1)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
loo2 = loo(fit2)
loo::loo_compare(loo1, loo2)
#>   model elpd_diff se_diff p_worse diag_diff      diag_elpd
#>  model1       0.0     0.0      NA           3 k_psis > 0.7
#>  model2      -2.0     3.2    0.73   N < 100               
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
