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
  pointwise = FALSE,
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

- pointwise:

  `TRUE` calls calls
  [`loo.function`](https://mc-stan.org/loo/reference/loo.html) which is
  slower but more memory efficient. `FALSE` calls the default
  [`loo`](https://mc-stan.org/loo/reference/loo.html).

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
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    uses only `sigma` for residuals.

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
cross-validation, which are not currently implemented in mcp.

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
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 60
#>    Unobserved stochastic nodes: 5
#>    Total graph size: 991
#> 
#> Initializing model
#> 
#> Finished sampling in 1.4 seconds
#> Warning: Some parameters may not have converged well:
#>   * Rhat > 1.01: cp_1 and x_1
#>   * ess_bulk or ess_tail < 400: Intercept_2 and cp_1 and x_1
#> Inspect `summary(fit)` and `plot_pars(fit)`, and consider increasing `iter`/`adapt` or simplifying the model before trusting these results.
fit2 = mcp(model2, data)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 60
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 568
#> 
#> Initializing model
#> 
#> Finished sampling in 0.7 seconds

# Compute LOO for each and compare (works for waic(fit) too)
fit1$loo = loo(fit1)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
fit2$loo = loo(fit2)
loo::loo_compare(fit1$loo, fit2$loo)
#>   model elpd_diff se_diff p_worse diag_diff      diag_elpd
#>  model1       0.0     0.0      NA           2 k_psis > 0.7
#>  model2      -1.5     3.2    0.68   N < 100               
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
