# Returns a data.frame with all combos of predictors

**\[experimental\]**

## Usage

``` r
interpolate_newdata(fit, by = NULL, x_values = NULL, at = NULL, arma = NULL)
```

## Arguments

- fit:

  An `mcpfit` object.

- by:

  Character vector of categorical or group-level columns to evaluate
  separately. Categorical model predictors are always included.

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

## Value

`data.frame` with

- Cols for par_x

- unique levels combos of factorial vars

- fixed values for additional continuous predictors

## Details

This function synthesizes predictors for all combinations of predictor
values. It is used internally in
[`plot.mcpfit()`](https://lindeloev.github.io/mcp/reference/plot.mcpfit.md)
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
[`log_lik()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
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
#>        time   fitted        sd     Q2.5    Q97.5
#> 1 0.6186323 10.04702 0.6479304 8.782501 11.32416
#> 2 1.6198715 10.04702 0.6479304 8.782501 11.32416
#> 3 2.6211108 10.04702 0.6479304 8.782501 11.32416
#> 4 3.6223501 10.04702 0.6479304 8.782501 11.32416
#> 5 4.6235893 10.04702 0.6479304 8.782501 11.32416
#> 6 5.6248286 10.04702 0.6479304 8.782501 11.32416

# Predictions for each draw
prediction = predict(demo_fit, newdata, summary = FALSE)
head(prediction)
#> # A tibble: 6 × 13
#>   .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3 time_3
#>    <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>  <dbl>
#> 1      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 2      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 3      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 4      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 5      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 6      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> # ℹ 4 more variables: sigma_1 <dbl>, time <dbl>, .prediction <dbl>,
#> #   data_row <int>

# Custom plot
library(ggplot2)
plotdata = fitted(demo_fit, newdata)
ggplot(plotdata, aes(x = time, y = fitted)) +
  geom_ribbon(aes(ymin = `Q2.5`, ymax = `Q97.5`), alpha = 0.3) +
  geom_line(lwd = 2) +
  geom_point(aes(y = response), data = demo_fit$data)

# }
```
