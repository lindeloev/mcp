# Fits and Predictions given Draws and data

Fits and Predictions given Draws and data

## Usage

``` r
pp_eval(
  object,
  newdata = NULL,
  summary = TRUE,
  type = "fitted",
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  dpar = "epred",
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  samples_format = "tidy",
  scale = "response",
  .include_fitted = FALSE,
  nsamples = lifecycle::deprecated()
)
```

## Arguments

- object:

  An `mcpfit` object.

- newdata:

  A `tibble` or a `data.frame` containing predictors in the model.
  Weighted Gaussian predictions and log-likelihoods also require the
  weights column. If `NULL` (default), the original data is used.

- summary:

  Summarise at each x-value

- type:

  One of:

  - `"fitted"`: return expected values. When `dpar` is the name of a
    dpar (e.g., `"mu"` or `"sigma"`), the expected value for just this
    dpar is returned. See also
    [`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).

  - `"predict"`: return predicted values (e.g.,
    `y_predict = rnorm(N, y_fitted, sigma_fitted)` for
    `family = gaussian()`). See also
    [`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).

  - `"residuals"`: observed y-values minus the fitted values. See also
    [`residuals()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).

  - `"loglik"`: return the log-likelihood for each sample for each data
    point. See also
    [`log_lik()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).
    Requires `scale = "response"`.

- probs:

  Vector of quantiles. Only in effect when `summary == TRUE`.

- rate:

  Boolean. For binomial models, plot on raw data (`rate = FALSE`) or
  response divided by number of trials (`rate = TRUE`). If FALSE, linear
  interpolation on trial number is used to infer trials at a particular
  x.

- prior:

  TRUE/FALSE. Plot using prior samples? Useful for
  `mcp(..., sample = "both")`

- dpar:

  What distributional parameter to evaluate. This is only relevant when
  `type == "fitted"`. E.g.,

  - `"epred"` (default): Expected value of the full model (or `NULL` for
    compatibility with brms etc.).

  - `"mu"`: The central tendency which is often the mean after applying
    the link function.

  - `"sigma"`: The standard deviation of the residuals.

  - `"ar1"`, `"ar2"`, `"ma1"`, `"ma2"`, etc. depending on which AR or MA
    coefficient you want to evaluate.

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

  Integer or `NULL`. Number of posterior draws to return/summarise. If
  there are varying effects, this is the number of draws from each
  varying group. `NULL` means "all". Ignored if both are `FALSE`. More
  samples trade speed for accuracy.

- samples_format:

  One of "tidy" or "matrix". Controls the output format when
  `summary == FALSE`. See more under "value"

- scale:

  One of

  - `"response"`: return on the observed scale, i.e., after applying the
    inverse link function.

  - `"linear"`: return on the parameter scale (where the linear trends
    are modelled). A linear scale is only applicable when
    `type == "fitted"` and `dpar` is not `NULL`.

- .include_fitted:

  Internal. Include fitted values with unsummarised predictions.

- nsamples:

  Deprecated. Use `ndraws` instead.

## Value

- If `summary = TRUE`: A `tibble` with the posterior mean for each row
  in `newdata`, If `newdata` is `NULL`, the data in `fit$data` is used.

- If `summary = FALSE` and `samples_format = "tidy"`: A `tidybayes`
  `tibble` with all the posterior draws (`Nd`) evaluated at each row in
  `newdata` (`Nn`), i.e., with `Nd x Nn` rows. If there are varying
  effects, the returned data is expanded with the relevant levels for
  each row.

  The return columns are:

  - Predictors from `newdata`.

  - Draw descriptors: ".chain", ".iteration", ".draw" (see the
    `posterior` and `tidybayes` packages), and `data_row`, the row
    number in the evaluated `newdata`.

  - Draw values: one column for each parameter in the model.

  - The estimate. Either ".epred", ".prediction", ".residual", or
    ".loglik" (matching tidybayes/ggdist conventions).

- If `summary = FALSE` and `samples_format = "matrix"`: An `N_draws` X
  `nrows(newdata)` matrix with fitted/predicted values (depending on
  `type`). This format is used by `brms` and it's useful as `yrep` in
  `bayesplot::ppc_*` functions.

## See also

[`fitted.mcpfit`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
[`predict.mcpfit`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
[`residuals.mcpfit`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
