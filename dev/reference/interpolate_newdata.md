# Returns a data.frame with all combos of predictors

This function synthesizes predictors for all combinations of predictor
values. It is used internally in
[`plot.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
and may be useful if you want to build your own custom plot.

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

  Character vector of categorical or varying-effect columns to evaluate
  separately. Categorical model predictors are always included.

- x_values:

  Numeric vector of x-values to evaluate at.

- at:

  Named list setting additional continuous predictors to fixed values.
  They default to their observed means. For example,
  `at = list(age = 40)`.

## Value

`tibble` with

- Cols for par_x

- unique levels combos of factorial vars

- fixed values for additional continuous predictors

## Details

The `par_x` variable will be interpolated with higher resolution around
the change points where the values can change abruptly, but lower
resolution in between to speed up the computation.

Categorical variables and requested varying-effect groups are combined
factorially (all level combinations). Additional continuous predictors
are held at their observed means, or at values supplied through `at`.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
if (FALSE) { # \dontrun{
# Get predictors for a fit
fit = mcp_example("multiple")
newdata = interpolate_newdata(fit)

# Fit summary
fitted(fit, newdata)

# Predictions for each sample
prediction = predict(fit, newdata, summary = FALSE)
prediction[, c(".chain", ".iteration", ".draw", "x", "group", "z", "predict")]

# Custom plot
library(ggplot2)
newdata = interpolate_newdata(fit)
plotdata = fitted(fit, newdata)
ggplot(plotdata, aes(x = x, y = fitted, color = group)) +
  geom_ribbon(aes(ymin = `Q2.5`, ymax = `Q97.5`, fill = group), alpha = 0.3) +
  geom_line(lwd = 2) +
  geom_point(aes(y = y), data = fit$data)
} # }
```
