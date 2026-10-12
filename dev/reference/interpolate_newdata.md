# Returns a tibble with all combos of predictors

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

  Logical. If `TRUE`, preserve the observed response history for
  conditional AR/MA evaluation. Defaults to `TRUE` when the fit includes
  [`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()` terms. Set to
  `FALSE` to interpolate unconditional trends.

- prior:

  Logical. Increase the resolution around the change points of the prior
  draws (`TRUE`) instead of the posterior draws (`FALSE`, default).

## Value

A tibble with

- Cols for par_x

- unique levels combos of factorial vars

- fixed values for additional continuous predictors

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
#> # A tibble: 6 × 5
#>    time .mean   .sd .lower .upper
#>   <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1 0.619  9.95 0.665   8.49   11.2
#> 2 1.62   9.95 0.665   8.49   11.2
#> 3 2.62   9.95 0.665   8.49   11.2
#> 4 3.62   9.95 0.665   8.49   11.2
#> 5 4.62   9.95 0.665   8.49   11.2
#> 6 5.62   9.95 0.665   8.49   11.2

# Predictions for each draw
prediction = predict(demo_fit, newdata, summary = FALSE)
head(prediction)
#> # A tibble: 6 × 13
#>   .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3 time_3
#>    <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>  <dbl>
#> 1      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 2      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 3      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 4      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 5      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 6      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> # ℹ 4 more variables: sigma_1 <dbl>, time <dbl>, .data_row <int>,
#> #   .prediction <dbl>

# Custom plot
library(ggplot2)
plotdata = fitted(demo_fit, newdata)  # `width = .95` by default 
ggplot(plotdata, aes(x = time, y = .mean)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(lwd = 2) +
  geom_point(aes(y = response), data = demo_fit$data)

# }
```
