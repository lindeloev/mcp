# Underlies `plot()` and `plot_dpar()`

Underlies
[`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
and
[`plot_dpar()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)

## Usage

``` r
get_plot(
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
  dpar = "epred",
  arma = TRUE,
  ndraws = 1000,
  scale = "response",
  at = NULL,
  .grouping = "auto",
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

  - A vector of quantiles. For example, `quantiles = 0.5` plots the
    median and `quantiles = c(0.2, 0.8)` plots the 20% and 80%
    quantiles.

- q_predict:

  Same as `q_fit`, but for the prediction interval.

- facet_by:

  Character vector. Names of categorical data columns to split to
  facets. Can be grouping-factor or categorical predictor columns.

- color_by:

  A character vector naming categorical or varying-effect data columns
  to color by. If both `color_by` and `facet_by` are omitted, a sole
  categorical predictor is colored automatically. Set `color_by = NULL`
  explicitly to disable this. Multiple columns are combined as an
  interaction. Curves and quantiles remain separate for grouping columns
  not mapped to color.

- lines:

  Positive integer or `FALSE`. The number of fitted lines (draws). It is
  the number of joint posterior draws shown for every curve. FALSE or
  `lines = 0` plots no lines. Note that lines always plot fitted
  values - not predicted. For prediction intervals, see the `q_predict`
  argument.

- geom_data:

  String. One of "point", "line" (good for time-series), or FALSE (do
  not plot).

- cp_dens:

  TRUE/FALSE. Plot posterior densities of the change point(s)? Currently
  does not respect `facet_by`. This will be added in the future.

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

- arma:

  Whether to include AR and MA effects.

  - `TRUE` Compute the GARMA residual recurrence. Requires the response
    variable in `newdata`.

  - `FALSE` Disregard AR and MA effects. For `family = gaussian()`,
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    uses only `sigma` for residuals.

- ndraws:

  Integer or `NULL`. Number of posterior draws to return/summarise. If
  there are varying effects, this is the number of draws from each
  varying group. `NULL` means "all". Ignored if both are `FALSE`. More
  samples trade speed for accuracy.

- scale:

  One of

  - `"response"`: return on the observed scale, i.e., after applying the
    inverse link function.

  - `"linear"`: return on the parameter scale (where the linear trends
    are modelled). A linear scale is only applicable when
    `type == "fitted"` and `dpar` is not `NULL`.

- at:

  Named list setting additional continuous predictors to fixed values.
  They default to their observed means. Passed to
  [`interpolate_newdata()`](https://lindeloev.github.io/mcp/dev/reference/interpolate_newdata.md).

- .grouping:

  Internal. Whether grouping arguments were omitted, mapped, or
  explicitly disabled.

- nsamples:

  Deprecated. Use `ndraws` instead.

- ...:

  Currently ignored.

## Value

A ggplot2 object.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
