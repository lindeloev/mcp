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
  rate = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  ...
)

# S3 method for class 'mcpfit'
fitted(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  dpar = "epred",
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = "response",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  ...
)

# S3 method for class 'mcpfit'
log_lik(
  object,
  newdata = NULL,
  summary = FALSE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "matrix",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  ...
)

# S3 method for class 'mcpfit'
residuals(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated(),
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
    `fitted()`, `residuals()`, `log_lik()`, and posterior `predict()`
    condition on the response history, so `newdata` must include the
    response. For `fitted()`, `predict()`, and `residuals()`, missing
    response histories are supported only in the original fitted data,
    using retained posterior imputations. Prior `predict()` and
    [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
    generate fresh response series recursively, so their `newdata` need
    only contain predictors. `log_lik()` is unavailable when a missing
    response enters a later observed history.

  - For models with `y | weights()`: Require the weights column except
    for `fitted()` and `predict()`.

- summary:

  Summarise at each x-value

- probs:

  Vector of quantiles. Only in effect when `summary == TRUE`.

- rate:

  Logical scalar. For binomial models, return counts (`rate = FALSE`) or
  the observed or expected success proportion (`rate = TRUE`).
  Predictions and count-scale fitted values require a trials column in
  `newdata`. Distributional parameters such as `dpar = "mu"` evaluate
  the parameter itself (e.g., success probability) and are unaffected by
  `rate`.

- prior:

  Logical. Evaluate prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)? Useful for `mcp(..., sample = "both")`.

- varying:

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
    currently requires including all such effects (`varying = TRUE`).

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

  - `"response"`: return on the observed scale, i.e., after applying the
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

`residuals(fit)` is equivalent to
`fit$data[[mcp_columns(fit)$response]] - fitted(fit, ...)` (or
`newdata[[mcp_columns(fit)$response]] - fitted(fit, ...)`), but with
fixed arguments for `fitted`:
`rate = FALSE, dpar = 'epred', draws_format = 'tidy'`.

`log_lik()` defaults to an unsummarised draws-by-observation matrix, as
used by `loo` and other posterior workflows. Non-default `varying` and
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
uses retained JAGS imputations for their posterior response draws. In
GARMA models these imputations also supply the history used for later
fitted and predicted rows.

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
#> Error in family_log_lik(y, dpars, data): could not find function "family_log_lik"

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
#> # ℹ 5 more variables: sigma_1 <dbl>, response <dbl>, time <dbl>, .epred <dbl>,
#> #   data_row <int>
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
#>    y  x condition   fitted        sd     Q2.5    Q97.5
#> 1 NA  8         B 35.60043 1.0803038 33.51261 37.74217
#> 2 NA 19         A 15.68747 0.8392979 14.05176 17.34367
#> 3 NA 27         A 17.20067 0.7620321 15.70353 18.68519
#> 4 NA 28         B 39.38345 0.7593611 37.89912 40.87640
#> 5 NA 29         A 17.57898 0.7589502 16.07714 19.05693
#> 6 NA 30         B 39.76175 0.7585016 38.28224 41.25734
fitted(missing_fit, summary = FALSE) |> dplyr::filter(is.na(y)) |> head()  # Same, but draws
#> # A tibble: 6 × 14
#>   .chain .iteration .draw  cp_1 Intercept_1   x_1 conditionB_1    x_2 sigma_1
#>    <int>      <int> <int> <dbl>       <dbl> <dbl>        <dbl>  <dbl>   <dbl>
#> 1      1          1     1  62.8        14.8 0.135         20.3 -0.497    4.65
#> 2      1          1     1  62.8        14.8 0.135         20.3 -0.497    4.65
#> 3      1          1     1  62.8        14.8 0.135         20.3 -0.497    4.65
#> 4      1          1     1  62.8        14.8 0.135         20.3 -0.497    4.65
#> 5      1          1     1  62.8        14.8 0.135         20.3 -0.497    4.65
#> 6      1          1     1  62.8        14.8 0.135         20.3 -0.497    4.65
#> # ℹ 5 more variables: y <dbl>, x <int>, condition <fct>, .epred <dbl>,
#> #   data_row <int>
predict(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Posterior predictive for missing y
#>    y  x condition  predict       sd      Q2.5    Q97.5
#> 1 NA  8         B 35.51722 4.408885 26.944898 44.27282
#> 2 NA 19         A 15.70005 4.409613  7.131763 24.25305
#> 3 NA 27         A 17.26378 4.343011  8.668505 25.73409
#> 4 NA 28         B 39.37180 4.354504 30.849440 47.91352
#> 5 NA 29         A 17.55903 4.343115  9.046847 26.11032
#> 6 NA 30         B 39.73816 4.364214 31.227062 48.29058
# }
```
