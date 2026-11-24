# Tibble with All Combinations of Predictors

**\[experimental\]**

## Usage

``` r
interpolate_newdata(
  fit,
  by = NULL,
  x_values = NULL,
  at = NULL,
  arma = NULL,
  prior = FALSE
)
```

## Arguments

- fit:

  An `mcpfit` object.

- by:

  Character vector of categorical or group-level columns to evaluate
  separately. Categorical model predictors are always evaluated
  separately.

- x_values:

  Numeric vector of x-values to evaluate at.

- at:

  Named list setting additional continuous predictors to fixed values.
  They default to their observed means. Family response auxiliaries can
  also be supplied as explicit scalar design values; e.g.,
  `at = list(N = 20)`.

- arma:

  How AR/MA terms will be evaluated (see
  [`predict.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)).
  Defaults to `TRUE` when the fit includes
  [`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()` terms.

  - `TRUE`: Preserve the observed rows and response history for
    conditional AR/MA evaluation.

  - `FALSE`: Interpolate unconditional trends.

  - `"replicate"`: Preserve the observed rows, because AR/MA terms
    depend on the previous rows.

- prior:

  Logical. Increase the resolution around the change points of the prior
  draws (`TRUE`) instead of the posterior draws (`FALSE`, default).

## Value

A tibble with

- A column for the change-point variable (`par_x`)

- All combinations of categorical predictor levels

- Fixed values for additional continuous predictors

## Details

This function synthesizes predictors for all combinations of predictor
values. It is used internally in
[`plot.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
and may be useful if you want to build your own custom plot.

The `par_x` variable will be interpolated with higher resolution around
the change points where the values can change abruptly, but lower
resolution in between to speed up the computation.

Categorical variables and requested grouping factors are combined
factorially (all level combinations). Additional continuous predictors
are held at their observed means, or at values supplied through `at`.
Family-specific response auxiliaries are not interpolated. Supply
binomial trial counts as a scalar in `at` or use `newdata` for a varying
design.

Likelihood weights are needed only when evaluating
[`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
and therefore need not be supplied here.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# \donttest{
# Get predictors for a fit
newdata = interpolate_newdata(demo_fit)

# Fit summary
head(fitted(demo_fit, newdata))
#> # A tibble: 6 × 6
#>    time  .row .mean   .sd .lower .upper
#>   <dbl> <int> <dbl> <dbl>  <dbl>  <dbl>
#> 1 0.619     1  9.92 0.658   8.60   11.2
#> 2 1.62      2  9.92 0.658   8.60   11.2
#> 3 2.62      3  9.92 0.658   8.60   11.2
#> 4 3.62      4  9.92 0.658   8.60   11.2
#> 5 4.62      5  9.92 0.658   8.60   11.2
#> 6 5.62      6  9.92 0.658   8.60   11.2

# Predictions for each draw
prediction = predict(demo_fit, newdata, summary = FALSE)
head(prediction)
#> # A tibble: 6 × 13
#>    time  .row .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3  time_3 sigma_1 .prediction
#>   <dbl> <int>  <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>   <dbl>   <dbl>       <dbl>
#> 1 0.619     1      1          1     1  32.5  69.6       10.2   0.582        10.8 -0.0577    4.49        9.55
#> 2 0.619     1      1          2     2  29.4  72.3        8.48  0.582        13.9 -0.197     4.14        6.36
#> 3 0.619     1      1          3     3  29.0  69.6        9.49  0.521        14.0 -0.199     3.99        9.89
#> 4 0.619     1      1          4     4  32.8  69.7       11.0   0.521        14.8 -0.277     4.05        9.11
#> 5 0.619     1      1          5     5  33.9  69.5       10.2   0.609        15.2 -0.276     4.16       17.9 
#> 6 0.619     1      1          6     6  30.4  70.2        9.94  0.539        15.0 -0.226     4.04        8.00

# Custom plot
library(ggplot2)
plotdata = fitted(demo_fit, newdata)  # `width = .95` by default 
ggplot(plotdata, aes(x = time, y = .mean)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(lwd = 2) +
  geom_point(aes(y = response), data = demo_fit$data)

# }
```
