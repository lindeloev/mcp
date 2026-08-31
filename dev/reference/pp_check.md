# Posterior Predictive Checks For Mcpfit Objects

Plot posterior (default) or prior (`prior = TRUE`) predictive checks.
This is convenience wrapper around the `bayesplot::ppc_*()` methods.

## Usage

``` r
pp_check(
  object,
  type = "dens_overlay",
  facet_by = NULL,
  newdata = NULL,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = 100,
  nsamples = lifecycle::deprecated(),
  ...
)
```

## Arguments

- object:

  An `mcpfit` object.

- type:

  String specifying the type of plot. See
  [`bayesplot::available_ppc()`](https://mc-stan.org/bayesplot/reference/PPC-overview.html)
  for available plot types, omitting the `"ppc_"` prefix (e.g.,
  `"dens_overlay"`, `"ribbon"`, `"intervals"`).

- facet_by:

  Name of a grouping column used by group-level effects.

- newdata:

  A `tibble` or a `data.frame` containing predictors in the model.

  - If `NULL` (default), the original data is used.

  - For models with [`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()`:
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    or
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    conditions on the response history, so `newdata` must include the
    response. For
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    missing response histories are supported only in the original fitted
    data, using retained posterior imputations.
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    is unavailable when a missing response enters a later observed
    history. Use
    [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
    to generate fresh response series recursively from predictor-only
    `newdata`.

  - For models with `y | weights()`: Require the weights column except
    for
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    and
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).

- prior:

  Logical. Evaluate prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)? Useful for `mcp(..., sample = "both")`.

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

  Integer or `NULL`. Number of posterior draws to return/summarise. If
  there are group-level effects, this is the number of draws from each
  group. `NULL` means "all". More draws trade speed for accuracy.

- nsamples:

  Deprecated. Use `ndraws` instead.

- ...:

  Further arguments passed to `bayesplot::ppc_type(y, yrep, ...)`

## Value

A `ggplot2` object for single plots. Enriched by `patchwork` for faceted
plots.

## Details

Missing responses are omitted from the observed-data check. For GARMA
models, each replicated response series is generated recursively rather
than conditioning on the observed response history. LOO predictive
checks use posterior draws and the original fitted data, so they do not
support `prior = TRUE` or `newdata`.

## See also

[`plot.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
[`pp_eval`](https://lindeloev.github.io/mcp/dev/reference/pp_eval.md)

plot_pars plot_dpar pp_check

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# \donttest{
pp_check(demo_fit)

pp_check(demo_fit, type = "ecdf_overlay")

#pp_check(some_group_fit, type = "loo_intervals", facet_by = "id")
# }
```
