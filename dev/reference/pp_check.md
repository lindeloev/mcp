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

  One of
  `bayesplot::available_ppc("grouped", invert = TRUE) %>% stringr::str_remove("ppc_")`

- facet_by:

  Name of a column in data modeled as varying effect(s).

- newdata:

  A `tibble` or a `data.frame` containing predictors in the model.
  Weighted Gaussian predictions and log-likelihoods also require the
  weights column. If `NULL` (default), the original data is used.

- prior:

  TRUE/FALSE. Plot using prior samples? Useful for
  `mcp(..., sample = "both")`

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

  Number of posterior draws. Note that you may want to use all draws for
  summary geoms, e.g., `pp_check(fit, type = "ribbon", ndraws = NULL)`.
  LOO checks always evaluate all posterior draws to preserve their PSIS
  weights; where supported, `ndraws` is passed to bayesplot to control
  the number of plotted samples.

- nsamples:

  Deprecated. Use `ndraws` instead.

- ...:

  Further arguments passed to `bayesplot::ppc_type(y, yrep, ...)`

## Value

A `ggplot2` object for single plots. Enriched by `patchwork` for faceted
plots.

## Details

Missing responses are omitted from the observed-data check. LOO
predictive checks use posterior draws and the original fitted data, so
they do not support `prior = TRUE` or `newdata`.

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

#pp_check(some_varying_fit, type = "loo_intervals", facet_by = "id")
# }
```
