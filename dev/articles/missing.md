# Missing responses and imputation

`mcp` allows missing values in the response. JAGS treats the unknown
responses as latent variables and samples them together with the model
parameters. `mcp` retains those posterior imputation draws for later
use. This article shows how to inspect expected responses, posterior
imputations, and their uncertainty.

## An example with missing responses

The built-in example has a change in slope, a strong difference between
two conditions, scattered missing responses, and a short run of missing
responses. First, let’s fit the example and visualize it:

``` r

library(mcp)
library(dplyr)
library(ggplot2)
set.seed(42)
```

``` r

fit = mcp_example("missing")
```

![](missing_files/figure-html/unnamed-chunk-2-1.png)

The ordinary plot shows the observed responses and model estimates.
Missing responses are omitted. See below how they can be visualized.

The underlying data contain some missing responses in `y`:

``` r

fit$data |> filter(is.na(y))
```

    ##     y  x condition
    ## 1  NA  8         B
    ## 2  NA 19         A
    ## 3  NA 27         A
    ## 4  NA 28         B
    ## 5  NA 29         A
    ## 6  NA 30         B
    ## 7  NA 31         A
    ## 8  NA 68         B
    ## 9  NA 84         B
    ## 10 NA 96         B

## Expected responses and imputations

You can see the imputations fairly directly.
[`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
describes plausible values of the missing response itself. These are the
latent response draws sampled jointly with the model parameters. The
imputation comes from trends in non-missing data and what they imply,
given the observed covariates (`condition` and `x`), for the missing
rows.

``` r

# Quantiles of posterior predictive for missing data: 10%, median, 90%
predict(fit, probs = c(0.1, 0.5, 0.9)) |>
  filter(is.na(y))
```

    ##     y  x condition  predict       sd      Q10      Q50      Q90
    ## 1  NA  8         B 35.51722 4.408885 29.97320 35.59747 41.23138
    ## 2  NA 19         A 15.70005 4.409613 10.12722 15.68574 21.24990
    ## 3  NA 27         A 17.26378 4.343011 11.65769 17.20047 22.74389
    ## 4  NA 28         B 39.37180 4.354504 33.84050 39.38414 44.92549
    ## 5  NA 29         A 17.55903 4.343115 12.03646 17.57914 23.12128
    ## 6  NA 30         B 39.73816 4.364214 34.21878 39.76277 45.30338
    ## 7  NA 31         A 17.95321 4.351576 12.41368 17.95779 23.50021
    ## 8  NA 68         B 40.64902 4.399515 35.00915 40.62840 46.26328
    ## 9  NA 84         B 34.26827 4.374458 28.70075 34.26938 39.84238
    ## 10 NA 96         B 29.43116 4.495332 23.81245 29.49694 35.18704

[`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
has the same syntax as
[`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
but returns the posterior *expected* response. Its intervals are
narrower because they do not include response variation:

``` r

# Quantiles of posterior evaluated at missing data: 10%, median, 90%
fitted(fit, probs = c(0.1, 0.5, 0.9)) |>
  filter(is.na(y))
```

    ##     y  x condition   fitted        sd      Q10      Q50      Q90
    ## 1  NA  8         B 35.60043 1.0803038 34.22732 35.59874 36.97791
    ## 2  NA 19         A 15.68747 0.8392979 14.62587 15.68943 16.76656
    ## 3  NA 27         A 17.20067 0.7620321 16.22350 17.20167 18.16512
    ## 4  NA 28         B 39.38345 0.7593611 38.41313 39.38536 40.34947
    ## 5  NA 29         A 17.57898 0.7589502 16.60533 17.58348 18.53966
    ## 6  NA 30         B 39.76175 0.7585016 38.78779 39.76350 40.72794
    ## 7  NA 31         A 17.95728 0.7627685 16.98320 17.96328 18.92934
    ## 8  NA 68         B 40.63312 1.0764083 39.34544 40.55922 42.04241
    ## 9  NA 84         B 34.27065 0.8778946 33.12572 34.27083 35.39471
    ## 10 NA 96         B 29.49858 1.2543575 27.88809 29.51210 31.11228

## Visualize imputations on the model plot

Imputations are not added to
[`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
automatically. Since
[`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
returns a ggplot, it is easy to extend. Here an `×` marks the posterior
median and a vertical line shows the central 80% imputation interval:

``` r

imputed = predict(fit, probs = c(0.1, 0.5, 0.9)) |>
  filter(is.na(y))

# Start with mcp plot with prediction interval. Use all draws for less Monte Carlo error.
set.seed(42)
plot(fit, color_by = "condition", q_predict = c(0.1, 0.9)) +

  # Add posterior predictive interval
  geom_linerange(
    data = imputed,
    aes(x = x, ymin = Q10, ymax = Q90),
    inherit.aes = FALSE,
    color = "black"
  ) +

  # Add "x" at posterior median
  geom_point(
    data = imputed,
    aes(x = x, y = Q50),
    inherit.aes = FALSE,
    shape = 4,  # an "x"
    size = 2
  )
```

![](missing_files/figure-html/unnamed-chunk-6-1.png)

This deliberately distinguishes imputed values from the solid observed
points. You can change the quantiles, marker, color, or add a subset of
the missing rows using ordinary `dplyr` and `ggplot2` code.

The plot above reduces each posterior predictive distribution to a
simple interval. Use `summary = FALSE` when you need the individual
draws instead. Here we plot the posterior predictive distribution for
each missing row:

``` r

imputed_draws = predict(fit, summary = FALSE) |>
  filter(is.na(y))

ggplot(imputed_draws, aes(x = .prediction)) +
  geom_density() + 
  facet_wrap(~x)
```

![](missing_files/figure-html/unnamed-chunk-7-1.png)

### Making probabilistic statements about missing responses

Sometimes, the missing response is of direct interest. For example, you
may want to know the probability that a missing response is above a
threshold. You can use
[`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
with `summary = FALSE` to get all posterior draws and then calculate
probabilities of various statements. For example:

``` r

imputed_draws |>
  filter(data_row == 19) |>
  summarise(
    p_greater = mean(.prediction > 20),
    p_lower = mean(.prediction < 10),
    p_between = mean(.prediction > 10 & .prediction < 20)
  )
```

    ## # A tibble: 1 × 3
    ##   p_greater p_lower p_between
    ##       <dbl>   <dbl>     <dbl>
    ## 1     0.163  0.0973     0.739

### Other continuous predictors

Plotting differs from
[`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
and
[`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
in one respect: The smooth curves in
[`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
are evaluated on data made by
[`interpolate_newdata()`](https://lindeloev.github.io/mcp/dev/reference/interpolate_newdata.md),
which keeps additional continuous predictors (i.e., predictors not
assigned an aesthetic) at their observed means by default, or at values
supplied through the `at`-argument. The plot caption reports these fixed
values.

The imputation markers above are different: `predict(fit)` evaluates the
original data, so every missing row uses its own predictor values. If
another continuous predictor is especially important, an imputation may
therefore lie away from a curve drawn at that predictor’s mean. This is
expected, just as an observed response with unusual predictor values may
lie away from that curve. Use `plot(fit, at = ...)` or
`interpolate_newdata(fit, at = ...)` or make separate plots when those
differences deserve emphasis.

## What else can be modeled?

Missing responses work with the usual `mcp` model features. For example,
you can:

- Include [continuous and categorical
  predictors](https://lindeloev.github.io/mcp/dev/articles/formulas.md)
  so imputations follow observed covariate information.
- Use [group-level
  effects](https://lindeloev.github.io/mcp/dev/articles/group_effects.md)
  so sparsely observed groups borrow information from the population and
  their available observations.
- [Model changing residual standard deviation with
  `sigma()`](https://lindeloev.github.io/mcp/dev/articles/dpar.md) so
  imputation uncertainty can vary over the predictor range.
- Use [binomial, Bernoulli, Poisson, or negative-binomial
  responses](https://lindeloev.github.io/mcp/dev/articles/families.md).
  Their imputations remain on the appropriate discrete response scale;
  binomial predictions can also be returned as rates.

The predictors required by the model must themselves be observed. `mcp`
currently imputes missing responses, not missing predictor values, and
it does not model the process that caused responses to be missing. As
with ordinary regression using incomplete outcomes, interpretation
relies on the missing responses being reasonably explained by the
observed predictors and model structure.

## Model evaluation and time-series histories

Missing rows are not treated as observed contributions to
[`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
[`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
or
[`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md).
For ordinary non-AR/MA models, the remaining observed rows can still be
used for these calculations.

AR/MA models need extra care because a missing response can become part
of the history for later observations. `mcp` retains the JAGS draws so
[`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
and
[`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
can reconstruct that history, but it does not currently integrate over
missing histories for pointwise likelihood calculations. Consequently,
[`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
[`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
and
[`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
are unavailable when a missing response preceeds some observed data.

The [time-series
article](https://lindeloev.github.io/mcp/dev/articles/arma.md) discusses
this. Briefly, use `posterior_predict()` to generate fresh replicated
series from the posteriors (as all other mcp functions do), which do not
depend on the observed response history.

## Original data versus new data

The retained JAGS draws belong to the missing rows in the original
fitted data. Calling `predict(fit)` uses them. For the non-AR/MA example
in this article, predictions for genuinely new data are ordinary
posterior predictions rather than imputations of the original missing
responses:

``` r

newdata = data.frame(
  x = c(20, 80),
  condition = factor(c("A", "B"), levels = levels(fit$data$condition))
)

predict(fit, newdata = newdata, probs = c(0.1, 0.5, 0.9))
```

    ##    x condition  predict       sd      Q10      Q50      Q90
    ## 1 20         A 15.90036 4.343471 10.31987 15.87509 21.43529
    ## 2 80         B 35.84599 4.345775 30.30002 35.86017 41.42434
