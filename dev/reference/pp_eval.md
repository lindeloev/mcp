# Fits and Predictions given Draws and data

Fits and Predictions given Draws and data

## Usage

``` r
pp_eval(
  object,
  newdata = NULL,
  summary = TRUE,
  type = "fitted",
  width = 0.95,
  rate = TRUE,
  prior = FALSE,
  dpar = "epred",
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = "response",
  .return_fitted = FALSE,
  .return_dpars = FALSE,
  .garma_replicate = FALSE
)
```

## Arguments

- object:

  An `mcpfit` object.

- newdata:

  A `tibble` or a `data.frame` containing predictors in the model.

  - If `NULL` (default), the original data is used.

  - Grouping factors must only contain levels from the fitted data. Use
    `group_effects = FALSE` to predict without group-level effects,
    e.g., for new groups.

  - For models with [`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()`:
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    condition on the response history by default, so `newdata` must
    include the response. For
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    missing response histories are supported only in the original fitted
    data, using retained posterior imputations as histories. Predictions
    are fresh response draws, including at missing rows. With
    `conditional = FALSE`,
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    and
    [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
    generate fresh response series recursively, so their `newdata` need
    only contain predictors and required response auxiliaries.
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

- width:

  Width of the central interval (strictly between 0 and 1), returned in
  the `.lower` and `.upper` columns. `FALSE` omits the interval. Only in
  effect when `summary == TRUE`. To get custom quantiles, summarise the
  draws from `summary = FALSE`, e.g., using
  [`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html)
  (see examples).

- rate:

  Logical scalar. For binomial models, return counts (`rate = FALSE`,
  the default for
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  and
  [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md))
  or the observed or expected success proportion (`rate = TRUE`).
  Predictions and count-scale fitted values require a trials column in
  `newdata`. Distributional parameters such as `dpar = "mu"` evaluate
  the parameter itself (e.g., success probability) and are unaffected by
  `rate`.

- prior:

  Logical. Evaluate prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default). The selected draws must be available; prior-only
  fits require `prior = TRUE`.

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

- group_effects:

  Group-level effects to include. One of:

  - `TRUE` All group-level deviations.

  - `FALSE` or `NA` No group-level deviations, i.e., population-level
    predictions. Like `re_formula = NA` in `{brms}`.

  - `"cp"` or `"predictor"`: All group-level deviations belonging to
    that part of the model.

  - Character vector of group-level parameter names, e.g.,
    `c("cp_1_id", "x_2_id")`. `*` is a wildcard, e.g., `"state*_1_id"`
    selects all levels of a factor and `"*_id"` selects all effects by
    `id`. See names using `ranef(fit)`. Can be combined, e.g.,
    `c("cp", "sigma_*")`.

- arma:

  Whether to include AR and MA effects.

  - `TRUE` Compute the GARMA residual recurrence. Requires the response
    variable in `newdata`.

  - `FALSE` Disregard AR and MA effects. For `family = gaussian()`,
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    uses only `sigma` for residuals. For posterior evaluation of the
    original data, retained JAGS imputations supply missing GARMA
    histories. In models with group-level effects, this currently
    requires all such effects (`group_effects = TRUE`).

- ndraws:

  Integer or `NULL`. Number of posterior draws to return/summarise.
  `NULL` means "all". More draws trade speed for accuracy.

- draws_format:

  One of "tidy" or "matrix". Controls the output format when
  `summary == FALSE`.

- scale:

  One of

  - `"response"`: return on the response scale, i.e., after applying the
    inverse link function.

  - `"linear"`: return on the linear-predictor (link) scale, where the
    linear trends are modeled. A linear scale is only applicable when
    `type == "fitted"` and `dpar` is not `NULL`.

- .return_fitted:

  Internal. Include fitted values with unsummarised predictions.

- .return_dpars:

  Internal. Include distributional parameters and response data as
  attributes with unsummarised predictions.

- .garma_replicate:

  Internal. For GARMA predictions, generate each response history
  recursively instead of conditioning on observed responses.

## Value

- If `summary = TRUE`: A tibble with the columns of `newdata` (or
  `fit$data` if `newdata` is `NULL`), followed by summaries for each
  row:

  - `.mean`: The mean of the draws.

  - `.sd`: The SD of the draws. With posterior draws (the default), this
    is the posterior predictive SD for `type = "predict"` and the
    posterior SD of the evaluated quantity otherwise.

  - `.lower` and `.upper`: The bounds of the central interval given in
    `width`.

  With `prior = TRUE`, these are the analogous prior summaries. mcp's
  own columns are dot-prefixed to avoid collisions with data columns.

- If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes`
  `tibble` with all the posterior draws (`Nd`) evaluated at each row in
  `newdata` (`Nn`), i.e., with `Nd x Nn` rows. If there are group-level
  effects, the returned data is expanded with the relevant levels for
  each row.

  The return columns are:

  - Predictors from `newdata`, plus its response column when supplied.

  - Draw descriptors: ".chain", ".iteration", ".draw" (see the
    `posterior` and `tidybayes` packages), and `.data_row`, the row
    number in the evaluated `newdata`.

  - Draw values: one column for each parameter in the model.

  - The estimate (matching tidybayes conventions): `.epred` for expected
    responses, `.linpred` with `scale = "linear"`, the `.dpar` for other
    dpars (e.g., `.mu` or `.sigma`), `.prediction`, `.residual`, or
    `.loglik`.

- If `summary = FALSE` and `draws_format = "matrix"`: An `N_draws` X
  `nrows(newdata)` matrix with fitted/predicted values (depending on
  `type`). This format is used by `brms` and it's useful as `yrep` in
  `bayesplot::ppc_*` functions.

## See also

[`fitted.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
[`predict.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
[`residuals.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
