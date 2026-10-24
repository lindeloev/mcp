# Plot full fits

Plot prior or posterior model draws on top of data.

## Usage

``` r
# S3 method for class 'mcpfit'
plot(
  x,
  interval = "none",
  width = 0.95,
  facet_by = NULL,
  color_by = "auto",
  lines = 25,
  geom_data = "point",
  cp_dens = TRUE,
  rate = TRUE,
  prior = FALSE,
  arma = TRUE,
  ndraws = 500,
  scale = "response",
  at = NULL,
  q_fit = lifecycle::deprecated(),
  q_predict = lifecycle::deprecated(),
  nsamples = lifecycle::deprecated(),
  ...
)

plot_dpar(
  x,
  dpar = "epred",
  interval = "none",
  width = 0.95,
  facet_by = NULL,
  color_by = "auto",
  lines = 25,
  cp_dens = TRUE,
  prior = FALSE,
  arma = TRUE,
  ndraws = 500,
  scale = "response",
  at = NULL,
  ...
)
```

## Arguments

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object

- interval:

  Which central intervals to plot as dashed lines:

  - `"none"` (default): No intervals.

  - `"fitted"`: Interval of the fitted values (expected response), as in
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).
    Red lines.

  - `"predict"`: Posterior predictive interval, as in
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).
    Green lines.

  - `c("fitted", "predict")`: Both. `plot_dpar()` supports only `"none"`
    and `"fitted"`.

  For other quantiles, see [the plotting
  article](https://lindeloev.github.io/mcp/articles/plotting.html#making-your-own-plots)
  for computing them from `fitted(..., summary = FALSE)` or
  `predict(..., summary = FALSE)` and adding them as a ggplot2 layer.

- width:

  Width of the central interval(s) given in `interval` (strictly between
  0 and 1). For example, `width = 0.8` plots the 10% and 90% quantiles.

- facet_by:

  Character vector. Names of categorical data columns to split to
  facets. Can be grouping-factor or categorical predictor columns.

- color_by:

  A character vector naming categorical predictor or grouping columns to
  color by. Multiple columns are combined as an interaction. Curves and
  quantiles remain separate for grouping columns not mapped to color.
  `"auto"` colors by the categorical predictor if the model has exactly
  one and it is not in `facet_by`. `NULL` disables color.

- lines:

  Positive integer or `FALSE`. The number of fitted lines (draws). It is
  the number of joint posterior draws shown for every curve. FALSE or
  `lines = 0` plots no lines. Note that lines always plot fitted
  values - not predicted. For posterior predictive intervals, use
  `interval = "predict"`.

- geom_data:

  String. One of "point", "line" (good for time-series), or FALSE (do
  not plot).

- cp_dens:

  TRUE/FALSE. Plot posterior densities of the change point(s)?

- rate:

  Logical scalar. For binomial models, plot the observed or expected
  success proportion (`rate = TRUE`, the default here) or counts
  (`rate = FALSE`). This default is the opposite of
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  and
  [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  which return counts by default; plotting proportions keeps curves
  comparable across varying trial counts. Predictions and count-scale
  fitted values require a trials column in `newdata`/`at`. Other
  families do not support `rate = TRUE`, so the default only applies to
  binomial and Bernoulli models.

- prior:

  Logical. Evaluate prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default). The selected draws must be available; prior-only
  fits require `prior = TRUE`.

- arma:

  How to include AR and MA effects.

  - `TRUE`: Condition on the observed response history, i.e.,
    one-step-ahead evaluation of the GARMA residual recurrence. Requires
    the response variable in `newdata`.

  - `FALSE`: Disregard AR and MA effects. For `family = gaussian()`,
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    uses only `sigma` for residuals.

  - `"replicate"`: Generate fresh response series recursively instead of
    conditioning on the observed responses. Only for predictions:
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`posterior_predict()`](https://lindeloev.github.io/mcp/dev/reference/posterior_epred.mcpfit.md),
    and
    [`pp_check()`](https://lindeloev.github.io/mcp/dev/reference/pp_check.md)
    (where it is the default). Applies equally to prior and posterior
    draws.

  For posterior evaluation of the original data, retained JAGS
  imputations supply missing GARMA histories. In models with group-level
  effects, this currently requires all such effects
  (`group_effects = TRUE`).

- ndraws:

  Integer or `NULL`. Number of posterior draws to return/summarise.
  `NULL` means "all". More draws trade speed for accuracy.

- scale:

  One of

  - `"response"`: return on the response scale, i.e., after applying the
    inverse link function.

  - `"linear"`: return on the linear-predictor (link) scale, where the
    linear trends are modeled.

- at:

  Named list setting additional continuous predictors to fixed values.
  They default to their observed means. Family-specific response
  auxiliaries can be supplied as explicit scalar design values. Passed
  to
  [`interpolate_newdata()`](https://lindeloev.github.io/mcp/dev/reference/interpolate_newdata.md).

- q_fit, q_predict:

  Deprecated. Use `interval` and `width` instead.

- nsamples:

  Deprecated. Use `ndraws` instead.

- ...:

  Must be empty. Reserved for future use.

- dpar:

  What distributional parameter to evaluate. E.g.,

  - `"epred"` (default): Expected response from the full model (or
    `NULL` for compatibility with brms etc.).

  - `"mu"`: The conditional mean (or success probability per trial for
    binomial/bernoulli models), on the link or response scale.

  - `"sigma"`: The standard deviation of the residuals.

  - `"ar1"`, `"ar2"`, `"ma1"`, `"ma2"`, etc. depending on which AR or MA
    coefficient you want to evaluate.

## Value

A ggplot2 object.

## Details

`plot()` evaluates the model on posterior draws like
[`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
and
[`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).
These represent the (joint) posterior distribution. Interval summaries
are based on `ndraws` draws, and change-point densities are simple
densities on all available draws.

## Functions

- `plot_dpar()`: Plot distributional parameters

## See also

plot_pars plot_dpar pp_check

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# Typical usage. demo_fit is an mcpfit object.
plot(demo_fit)

# \donttest{
plot(demo_fit, prior = TRUE)  # The prior


plot(demo_fit, lines = 0, interval = "fitted")  # 95% central interval without lines

plot(demo_fit, interval = "predict", width = 0.8)  # 80% posterior predictive interval

plot_dpar(demo_fit, dpar = "sigma", lines = 100)  # Residual standard deviation on y


# Show a panel for each group-level effect
# plot(fit, facet_by = "my_column")

# Customize plots using regular ggplot2
library(ggplot2)
plot(demo_fit) + theme_bw(15) + ggtitle("Great plot!")

# }
```
