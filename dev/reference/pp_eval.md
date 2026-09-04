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
  draws_format = "tidy",
  scale = "response",
  .include_fitted = FALSE,
  .include_dpars = FALSE,
  .garma_replicate = FALSE,
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated()
)
```

## Arguments

- object:

  An `mcpfit` object.

- newdata:

  A `tibble` or a `data.frame` containing predictors in the model.

  - If `NULL` (default), the original data is used.

  - For models with [`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()`:
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and posterior
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    condition on the response history, so `newdata` must include the
    response. For
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    missing response histories are supported only in the original fitted
    data, using retained posterior imputations. Prior
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    and
    [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
    generate fresh response series recursively, so their `newdata` need
    only contain predictors.
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    is unavailable when a missing response enters a later observed
    history.

  - For models with `y | weights()`: Require the weights column except
    for
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    and
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).

- summary:

  Summarise at each x-value

- type:

  One of:

  - `"fitted"`: return the expected response. When `dpar` names a
    distributional parameter (e.g., `"mu"` or `"sigma"`), that parameter
    is returned instead. See also
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).

  - `"predict"`: return predicted values (e.g.,
    `y_predict = rnorm(N, y_fitted, sigma_fitted)` for
    `family = gaussian()`). See also
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).

  - `"residuals"`: observed y-values minus the fitted values. See also
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).

  - `"loglik"`: return the log-likelihood for each draw for each data
    point. See also
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).
    Requires `scale = "response"`.

- probs:

  Vector of quantiles. Only in effect when `summary == TRUE`.

- rate:

  Logical scalar. For binomial models, return counts (`rate = FALSE`) or
  the observed or expected success proportion (`rate = TRUE`).
  Predictions and count-scale fitted values require a trials column in
  `newdata`. Distributional parameters such as `dpar = "mu"` evaluate
  the parameter itself (e.g., success probability) and are unaffected by
  `rate`.

- prior:

  Logical. Evaluate prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)? Useful for `mcp(..., sample = "both")`.

- dpar:

  What distributional parameter to evaluate. This is only relevant when
  `type == "fitted"`. E.g.,

  - `"epred"` (default): Expected response from the full model (or
    `NULL` for compatibility with brms etc.).

  - `"mu"`: The conditional mean (or success probability per trial for
    binomial/bernoulli models), on the link or response scale.

  - `"sigma"`: The standard deviation of the residuals.

  - `"ar1"`, `"ar2"`, `"ma1"`, `"ma2"`, etc. depending on which AR or MA
    coefficient you want to evaluate.

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

- draws_format:

  One of "tidy" or "matrix". Controls the output format when
  `summary == FALSE` (for
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  and
  [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)).
  [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  always returns tidy output.

- scale:

  One of

  - `"response"`: return on the observed scale, i.e., after applying the
    inverse link function.

  - `"linear"`: return on the linear-predictor (link) scale, where the
    linear trends are modeled. A linear scale is only applicable when
    `type == "fitted"` and `dpar` is not `NULL`.

- .include_fitted:

  Internal. Include fitted values with unsummarised predictions.

- .include_dpars:

  Internal. Include distributional parameters and response data as
  attributes with unsummarised predictions.

- .garma_replicate:

  Internal. For GARMA predictions, generate each response history
  recursively instead of conditioning on observed responses.

- nsamples:

  Deprecated. Use `ndraws` instead.

- samples_format:

  Deprecated. Use `draws_format` instead. See more under "value"

## Value

- If `summary = TRUE`: A data frame with the draw mean and SD (`sd`) for
  each row in `newdata`. With posterior draws (the default), `sd` is the
  posterior predictive SD for `type = "predict"` and the posterior SD of
  the evaluated quantity otherwise. With `prior = TRUE`, these are the
  analogous prior summaries. If `newdata` is `NULL`, the data in
  `fit$data` is used.

- If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes`
  `tibble` with all the posterior draws (`Nd`) evaluated at each row in
  `newdata` (`Nn`), i.e., with `Nd x Nn` rows. If there are group-level
  effects, the returned data is expanded with the relevant levels for
  each row.

  The return columns are:

  - Predictors from `newdata`, plus its response column when supplied.

  - Draw descriptors: ".chain", ".iteration", ".draw" (see the
    `posterior` and `tidybayes` packages), and `data_row`, the row
    number in the evaluated `newdata`.

  - Draw values: one column for each parameter in the model.

  - The estimate. Either ".epred", ".prediction", ".residual", or
    ".loglik" (matching tidybayes/ggdist conventions).

- If `summary = FALSE` and `draws_format = "matrix"`: An `N_draws` X
  `nrows(newdata)` matrix with fitted/predicted values (depending on
  `type`). This format is used by `brms` and it's useful as `yrep` in
  `bayesplot::ppc_*` functions.

## See also

[`fitted.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
[`predict.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
[`residuals.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
