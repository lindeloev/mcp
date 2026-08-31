# Returns a data.frame with all combos of predictors

**\[experimental\]**

## Usage

``` r
interpolate_newdata(
  fit,
  by = NULL,
  x_values = get_x_values(fit, by),
  at = NULL
)
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

## Value

`data.frame` with

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
#>       time   fitted     error     Q2.5    Q97.5
#> 1 3.189323 10.29993 0.7167564 8.874782 11.66981
#> 2 4.146597 10.29993 0.7167564 8.874782 11.66981
#> 3 5.103871 10.29993 0.7167564 8.874782 11.66981
#> 4 6.061145 10.29993 0.7167564 8.874782 11.66981
#> 5 7.018419 10.29993 0.7167564 8.874782 11.66981
#> 6 7.975693 10.29993 0.7167564 8.874782 11.66981

# Predictions for each draw
prediction = predict(demo_fit, newdata, summary = FALSE)
head(prediction)
#> # A tibble: 6 × 13
#>   .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3 time_3
#>    <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>  <dbl>
#> 1      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 2      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 3      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 4      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 5      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 6      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> # ℹ 4 more variables: sigma_1 <dbl>, time <dbl>, data_row <int>,
#> #   .prediction <dbl>

# Custom plot
library(ggplot2)
plotdata = fitted(demo_fit, newdata)
ggplot(plotdata, aes(x = time, y = fitted)) +
  geom_ribbon(aes(ymin = `Q2.5`, ymax = `Q97.5`), alpha = 0.3) +
  geom_line(lwd = 2) +
  geom_point(aes(y = response), data = demo_fit$data)

# }
```
