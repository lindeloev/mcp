# Plot full fits

Plot prior or posterior model draws on top of data.

## Usage

``` r
# S3 method for class 'mcpfit'
plot(
  x,
  q_fit = FALSE,
  q_predict = FALSE,
  facet_by = NULL,
  color_by = NULL,
  lines = 25,
  geom_data = "point",
  cp_dens = TRUE,
  rate = TRUE,
  prior = FALSE,
  arma = TRUE,
  ndraws = 500,
  scale = "response",
  at = NULL,
  samples = lifecycle::deprecated(),
  nsamples = lifecycle::deprecated(),
  ...
)

plot_dpar(
  x,
  dpar = "epred",
  q_fit = FALSE,
  facet_by = NULL,
  color_by = NULL,
  lines = 25,
  cp_dens = TRUE,
  prior = FALSE,
  arma = TRUE,
  ndraws = 500,
  scale = "response",
  at = NULL,
  samples = lifecycle::deprecated(),
  nsamples = lifecycle::deprecated(),
  ...
)
```

## Arguments

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object

- q_fit:

  Whether to plot quantiles of the posterior (fitted value).

  - `TRUE` Add 2.5% and 97.5% quantiles. Corresponds to
    `q_fit = c(0.025, 0.975)`.

  - `FALSE` No quantiles

  - A vector of quantiles. For example, `q_fit = 0.5` plots the median
    and `q_fit = c(0.2, 0.8)` plots the 20% and 80% quantiles.

- q_predict:

  Same as `q_fit`, but for the posterior predictive interval.

- facet_by:

  Character vector. Names of categorical data columns to split to
  facets. Can be grouping-factor or categorical predictor columns.

- color_by:

  A character vector naming categorical predictor or grouping columns to
  color by. If both `color_by` and `facet_by` are omitted, a sole
  categorical predictor is colored automatically. Set `color_by = NULL`
  explicitly to disable this. Multiple columns are combined as an
  interaction. Curves and quantiles remain separate for grouping columns
  not mapped to color.

- lines:

  Positive integer or `FALSE`. The number of fitted lines (draws). It is
  the number of joint posterior draws shown for every curve. FALSE or
  `lines = 0` plots no lines. Note that lines always plot fitted
  values - not predicted. For posterior predictive intervals, see the
  `q_predict` argument.

- geom_data:

  String. One of "point", "line" (good for time-series), or FALSE (do
  not plot).

- cp_dens:

  TRUE/FALSE. Plot posterior densities of the change point(s)?

- rate:

  Logical scalar. For binomial models, return counts (`rate = FALSE`) or
  the observed or expected success proportion (`rate = TRUE`).
  Predictions and count-scale fitted values require a trials column in
  `newdata`.

- prior:

  Logical. Evaluate prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)? Useful for `mcp(..., sample = "both")`.

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

- scale:

  One of

  - `"response"`: return on the observed scale, i.e., after applying the
    inverse link function.

  - `"linear"`: return on the linear-predictor (link) scale, where the
    linear trends are modeled. A linear scale is only applicable when
    `type == "fitted"` and `dpar` is not `NULL`.

- at:

  Named list setting additional continuous predictors to fixed values.
  They default to their observed means. Family-specific response
  auxiliaries can be supplied as explicit scalar design values. Passed
  to
  [`interpolate_newdata()`](https://lindeloev.github.io/mcp/dev/reference/interpolate_newdata.md).

- samples, nsamples:

  Deprecated. Use `lines` instead.

- ...:

  Must be empty. Reserved for future use.

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

## Value

A ggplot2 object.

## Details

`plot()` uses `fit$simulate()` on posterior draws. These represent the
(joint) posterior distribution. Interval summaries are based on CDFs
computed on `ndraws` draws, and change-point densities are simple
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


plot(demo_fit, lines = 0, q_fit = TRUE)  # 95% central interval without lines

plot(demo_fit, q_predict = c(0.1, 0.9))  # 80% posterior predictive interval

plot_dpar(demo_fit, dpar = "sigma", lines = 100)  # Residual standard deviation on y


# Show a panel for each group-level effect
# plot(fit, facet_by = "my_column")

# Customize plots using regular ggplot2
library(ggplot2)
plot(demo_fit) + theme_bw(15) + ggtitle("Great plot!")

# }
```
