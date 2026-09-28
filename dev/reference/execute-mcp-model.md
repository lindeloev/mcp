# Fitted and predicted values of `mcp` models fits

Evaluate the model on data, either summarised (per data-row) or per
draw. You can use draws from the prior (`prior = TRUE`), select a
distributional parameter with `dpar`, and choose the response or
linear-predictor scale with `scale` where applicable.

## Usage

``` r
# S3 method for class 'mcpfit'
predict(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = FALSE,
  prior = FALSE,
  group = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  conditional = TRUE,
  ...
)

# S3 method for class 'mcpfit'
fitted(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = FALSE,
  prior = FALSE,
  dpar = "epred",
  group = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = "response",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  ...
)

log_lik(object, ...)

# S3 method for class 'mcpfit'
log_lik(
  object,
  newdata = NULL,
  summary = FALSE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  group = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "matrix",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  ...
)

# S3 method for class 'mcpfit'
residuals(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  prior = FALSE,
  group = TRUE,
  arma = TRUE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  ...
)
```

## Arguments

- object:

  An `mcpfit` object.

- newdata:

  A `tibble` or a `data.frame` containing predictors in the model.

  - If `NULL` (default), the original data is used.

  - For models with [`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()`:
    `fitted()`, `residuals()`, `log_lik()`, and `predict()` condition on
    the response history by default, so `newdata` must include the
    response. For `fitted()`, `predict()`, and `residuals()`, missing
    response histories are supported only in the original fitted data,
    using retained posterior imputations as histories. Predictions are
    fresh response draws, including at missing rows. With
    `conditional = FALSE`, `predict()` and
    [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
    generate fresh response series recursively, so their `newdata` need
    only contain predictors and required response auxiliaries.
    `log_lik()` is unavailable when a missing response enters a later
    observed history.

  - For models with `y | weights()`: Require the weights column except
    for `fitted()` and `predict()`.

- summary:

  Summarise at each x-value

- probs:

  Vector of quantiles (strictly between 0 and 1). Only in effect when
  `summary == TRUE`.

- rate:

  Logical scalar. For binomial models, return counts (`rate = FALSE`,
  the default for `fitted()` and `predict()`) or the observed or
  expected success proportion (`rate = TRUE`). Predictions and
  count-scale fitted values require a trials column in `newdata`.
  Distributional parameters such as `dpar = "mu"` evaluate the parameter
  itself (e.g., success probability) and are unaffected by `rate`.

- prior:

  Logical. Evaluate prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default). The selected draws must be available; prior-only
  fits require `prior = TRUE`.

- group:

  Group-level effects. One of:

  - `TRUE` All group-level deviations.

  - `FALSE` No group-level deviations
    ([`c()`](https://rdrr.io/r/base/c.html)).

  - `"cp"` or `"predictor"`: All group-level deviations belonging to
    that part of the model.

  - Character vector: Only include specified group-level parameters.

- arma:

  Whether to include AR and MA effects.

  - `TRUE` Compute the GARMA residual recurrence. Requires the response
    variable in `newdata`.

  - `FALSE` Disregard AR and MA effects. For `family = gaussian()`,
    `predict()` uses only `sigma` for residuals. For posterior
    evaluation of the original data, retained JAGS imputations supply
    missing GARMA histories. In models with group-level effects, this
    currently requires all such effects (`group = TRUE`).

- ndraws:

  Integer or `NULL`. Number of posterior draws to return/summarise. If
  there are group-level effects, this is the number of draws from each
  group. `NULL` means "all". More draws trade speed for accuracy.

- draws_format:

  One of "tidy" or "matrix". Controls the output format when
  `summary == FALSE` (for `fitted()`, `predict()`, and `log_lik()`).
  `residuals()` always returns tidy output.

- nsamples:

  Deprecated. Use `ndraws` instead.

- samples_format:

  Deprecated. Use `draws_format` instead. See more under "value"

- varying:

  Deprecated. Use `group` instead.

- conditional:

  Logical. For AR/MA models, condition on observed response histories
  (`TRUE`, default) or generate fresh histories recursively (`FALSE`).
  Applies equally to prior and posterior draws. Predictive checks with
  [`pp_check()`](https://lindeloev.github.io/mcp/dev/reference/pp_check.md)
  generate fresh histories.

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

- scale:

  One of

  - `"response"`: return on the response scale, i.e., after applying the
    inverse link function.

  - `"linear"`: return on the linear-predictor (link) scale, where the
    linear trends are modeled. A linear scale is only applicable when
    `type == "fitted"` and `dpar` is not `NULL`.

## Value

- If `summary = TRUE`: A data frame with the draw mean and SD (`sd`) for
  each row in `newdata`. With posterior draws (the default), `sd` is the
  posterior predictive SD for `type = "predict"` and the posterior SD of
  the evaluated quantity otherwise. With `prior = TRUE`, these are the
  analogous prior summaries. If `newdata` is `NULL`, the data in
  `fit$data` is used.

- If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes`
  `tibble` with all the posterior draws (`Nd`) evaluated at each row in
  `newdata` (`Nn`), i.e., with `Nd x Nn` rows. If there are group-level
  effects, the returned data is expanded with the relevant levels for
  each row.

  The return columns are:

  - Predictors from `newdata`, plus its response column when supplied.

  - Draw descriptors: ".chain", ".iteration", ".draw" (see the
    `posterior` and `tidybayes` packages), and `data_row`, the row
    number in the evaluated `newdata`.

  - Draw values: one column for each parameter in the model.

  - The estimate. Either ".epred", ".prediction", ".residual", or
    ".loglik" (matching tidybayes/ggdist conventions).

- If `summary = FALSE` and `draws_format = "matrix"`: An `N_draws` X
  `nrows(newdata)` matrix with fitted/predicted values (depending on
  `type`). This format is used by `brms` and it's useful as `yrep` in
  `bayesplot::ppc_*` functions.

## Details

`fitted()` and `posterior_epred()` evaluate the same expected responses;
`predict()` and `posterior_predict()` evaluate the same response
distributions. The `posterior_*()` methods return draws-by-observation
matrices, while `fitted()` and `predict()` summarise by default and also
offer tidy draws. For binomial models, the default response scale is
counts. Use `rate = TRUE` for proportions or `fitted(..., dpar = "mu")`
for success probabilities. During migration from v0.3.4, an omitted
`rate` warns once per session per function when counts differ from
proportions. Explicit `rate` settings do not warn.

`residuals(fit)` is equivalent to
`fit$data[[mcp_columns(fit)$response]] - fitted(fit, ...)` (or
`newdata[[mcp_columns(fit)$response]] - fitted(fit, ...)`), but with
fixed arguments for `fitted`:
`rate = FALSE, dpar = 'epred', draws_format = 'tidy'`.

`log_lik()` defaults to an unsummarised draws-by-observation matrix, as
used by `loo` and other posterior workflows. Non-default `group` and
`arma` settings evaluate conditional or counterfactual log-likelihoods
(e.g., omitting random effects or serial correlation); they cannot be
used in
[`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
or
[`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
because estimating information criteria for reduced models requires
refitting.

Missing responses in the original data remain missing in the response
column. `fitted()` returns their expected responses, while `predict()`
generates fresh posterior predictive response draws (including at
missing rows). In GARMA models, retained JAGS imputations supply the
history used for later fitted and predicted rows.

## Functions

- `predict(mcpfit)`: Predictive Distribution

- `fitted(mcpfit)`: Expected response

- `log_lik(mcpfit)`: Pointwise log-likelihood

- `residuals(mcpfit)`: Residual distribution

## See also

`fitted.mcpfit` `predict.mcpfit` `residuals.mcpfit` `log_lik.mcpfit`

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
head(fitted(demo_fit))  # Expected response for each row of demo_fit$data
#>   response     time   fitted        sd     Q2.5    Q97.5
#> 1 17.23552 76.33986 16.92982 0.8918820 15.24031 18.66061
#> 2 11.35171 83.51711 16.19732 0.7069182 14.86974 17.60133
#> 3 28.04995 60.18529 25.24068 0.9314782 23.35380 27.15740
#> 4 20.68198 74.72964 17.09416 0.9683574 15.30233 19.00162
#> 5 21.21364 85.88256 15.95591 0.7234506 14.61884 17.38066
#> 6 22.32282 40.05069 14.56256 0.6387450 13.29038 15.75679
head(residuals(demo_fit))  # Residuals for each row of demo_fit$data
#>   response     time  residuals        sd       Q2.5     Q97.5
#> 1 17.23552 76.33986  0.3056925 0.8918820 -1.4250920  1.995210
#> 2 11.35171 83.51711 -4.8456163 0.7069182 -6.2496242 -3.518034
#> 3 28.04995 60.18529  2.8092689 0.9314782  0.8925555  4.696150
#> 4 20.68198 74.72964  3.5878224 0.9683574  1.6803637  5.379652
#> 5 21.21364 85.88256  5.2577362 0.7234506  3.8329871  6.594805
#> 6 22.32282 40.05069  7.7602636 0.6387450  6.5660303  9.032439
log_lik(demo_fit)[1:3, 1:3]  # Log-likelihood at each demo_fit$data
#>              1         2         3
#> [1,] -2.451191 -3.099468 -2.613516
#> [2,] -2.246574 -2.937759 -2.546826
#> [3,] -2.251453 -2.882468 -2.666775

# All of the above take a range of arguments. E.g.,:
# \donttest{
head(predict(demo_fit))  # Pointwise posterior predictive
#>   response     time  predict       sd      Q2.5    Q97.5
#> 1 17.23552 76.33986 16.78289 3.945613  9.057670 24.78217
#> 2 11.35171 83.51711 16.26018 3.904999  8.405115 23.98530
#> 3 28.04995 60.18529 25.12394 4.087744 17.360018 33.12150
#> 4 20.68198 74.72964 17.22145 4.011694  9.187288 24.98071
#> 5 21.21364 85.88256 15.88880 4.203532  8.161476 23.75329
#> 6 22.32282 40.05069 14.64542 4.031344  6.798112 22.33360
head(predict(demo_fit, probs = c(0.1, 0.5, 0.9)))  # Median and 80% posterior predictive interval.
#>   response     time  predict       sd       Q10      Q50      Q90
#> 1 17.23552 76.33986 16.89906 4.072704 11.812794 16.93335 22.04235
#> 2 11.35171 83.51711 16.26356 4.193334 11.130297 16.19803 21.26364
#> 3 28.04995 60.18529 25.10203 4.039357 20.115521 25.24066 30.36594
#> 4 20.68198 74.72964 17.19236 4.045033 11.954065 17.09782 22.22954
#> 5 21.21364 85.88256 16.11136 4.208352 10.885880 15.95534 21.02687
#> 6 22.32282 40.05069 14.50406 3.786710  9.511571 14.56147 19.61455
head(predict(demo_fit, prior = TRUE))  # Prior predictive
#>   response     time  predict       sd      Q2.5    Q97.5
#> 1 17.23552 76.33986 13.62495 14.24374 -15.56523 40.81678
#> 2 11.35171 83.51711 13.82178 13.11392 -16.42007 40.76014
#> 3 28.04995 60.18529 14.49667 15.95062 -13.63652 42.57939
#> 4 20.68198 74.72964 14.03760 15.00494 -14.76847 41.25856
#> 5 21.21364 85.88256 13.87618 16.04855 -16.96870 40.83192
#> 6 22.32282 40.05069 14.14016 15.36999 -14.18878 41.58789
head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#> # A tibble: 6 × 14
#>   .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3 time_3
#>    <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>  <dbl>
#> 1      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 2      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 3      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 4      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 5      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> 6      1          1     1  30.8  72.1        10.2  0.513        15.7 0.0773
#> # ℹ 5 more variables: sigma_1 <dbl>, response <dbl>, time <dbl>,
#> #   data_row <int>, .epred <dbl>
head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#>   response     time   fitted        sd     Q2.5    Q97.5
#> 1 17.23552 76.33986 3.893444 0.2757184 3.381747 4.465763
#> 2 11.35171 83.51711 3.893444 0.2757184 3.381747 4.465763
#> 3 28.04995 60.18529 3.893444 0.2757184 3.381747 4.465763
#> 4 20.68198 74.72964 3.893444 0.2757184 3.381747 4.465763
#> 5 21.21364 85.88256 3.893444 0.2757184 3.381747 4.465763
#> 6 22.32282 40.05069 3.893444 0.2757184 3.381747 4.465763

# Evaluate at novel data
novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
head(predict(demo_fit, newdata = novel_data, probs = c(0.025, 0.5, 0.975)))
#>   time   predict        sd       Q2.5      Q50    Q97.5
#> 1   -5  9.873604  3.948793   2.278611 10.04644 17.81965
#> 2   20 10.209701  3.834257   2.278611 10.04644 17.81965
#> 3  300 -6.074700 16.263353 -40.535047 -3.98945 21.96670

# Work with missing responses
missing_fit = mcp_example("missing", plot = FALSE)
#> NA values detected in 'y'. JAGS will treat them as latent responses and impute them during sampling.
fitted(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Expected responses for missing y
#>    y  x state   fitted        sd     Q2.5    Q97.5
#> 1 NA  8     B 35.39478 1.0927946 33.30087 37.58317
#> 2 NA 19     A 15.50674 0.8459778 13.85304 17.14750
#> 3 NA 27     A 17.09444 0.7751878 15.57812 18.59652
#> 4 NA 28     B 39.36404 0.7655497 37.86660 40.89542
#> 5 NA 29     A 17.49136 0.7747784 15.97729 18.98954
#> 6 NA 30     B 39.76096 0.7661437 38.26464 41.30071
fitted(missing_fit, summary = FALSE) |> dplyr::filter(is.na(y)) |> head()  # Same, but draws
#> # A tibble: 6 × 14
#>   .chain .iteration .draw  cp_1 Intercept_1   x_1 stateB_1    x_2 sigma_1     y
#>    <int>      <int> <int> <dbl>       <dbl> <dbl>    <dbl>  <dbl>   <dbl> <dbl>
#> 1      1          1     1  59.8        12.1 0.189     21.5 -0.416    4.36    NA
#> 2      1          1     1  59.8        12.1 0.189     21.5 -0.416    4.36    NA
#> 3      1          1     1  59.8        12.1 0.189     21.5 -0.416    4.36    NA
#> 4      1          1     1  59.8        12.1 0.189     21.5 -0.416    4.36    NA
#> 5      1          1     1  59.8        12.1 0.189     21.5 -0.416    4.36    NA
#> 6      1          1     1  59.8        12.1 0.189     21.5 -0.416    4.36    NA
#> # ℹ 4 more variables: x <int>, state <fct>, data_row <int>, .epred <dbl>
predict(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Posterior predictive for missing y
#>    y  x state  predict       sd      Q2.5    Q97.5
#> 1 NA  8     B 35.41014 4.412756 26.745900 44.06132
#> 2 NA 19     A 15.53057 4.372454  6.959857 24.05974
#> 3 NA 27     A 17.02735 4.372192  8.570913 25.61888
#> 4 NA 28     B 39.39124 4.297813 30.844286 47.88610
#> 5 NA 29     A 17.46702 4.316559  8.967339 26.01502
#> 6 NA 30     B 39.72046 4.360684 31.240371 48.28261
# }
```
