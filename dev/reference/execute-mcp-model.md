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
  Weighted Gaussian predictions and log-likelihoods also require the
  weights column. If `NULL` (default), the original data is used. GARMA
  evaluation with `fitted()`, `predict()`, `residuals()`, or `log_lik()`
  conditions on the response history, so `newdata` must include the
  response. For `fitted()`, `predict()`, and `residuals()`, missing
  response histories are supported only in the original fitted data,
  using retained posterior imputations. `log_lik()` is unavailable when
  a missing response enters a later observed history. Use
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  to generate fresh response series recursively from predictor-only
  `newdata`.

- summary:

  Summarise at each x-value

- probs:

  Vector of quantiles. Only in effect when `summary == TRUE`.

- rate:

  Boolean. For binomial models, return counts (`rate = FALSE`) or the
  observed or expected success proportion (`rate = TRUE`). Predictions
  and count-scale fitted values require a trials column in `newdata`.

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

  Currently ignored.

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

- If `summary = TRUE`: A data frame with the draw mean and SD (`error`)
  for each row in `newdata`. With posterior draws (the default), `error`
  is the posterior predictive SD for `type = "predict"` and the
  posterior SD of the evaluated quantity otherwise. With `prior = TRUE`,
  these are the analogous prior summaries. If `newdata` is `NULL`, the
  data in `fit$data` is used.

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
used by `loo` and other posterior workflows.

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
#>     response     time    fitted     error      Q2.5       Q97.5
#> 1 -3.0084198 91.48060 -3.642469 0.7403902 -5.030038 -2.19864725
#> 2 -7.8768640 93.70754 -4.275666 0.8303914 -5.841824 -2.66551420
#> 3 16.3029101 28.61395 10.957474 1.0249532  9.055578 12.91323407
#> 4 -0.0373553 83.04476 -1.243859 0.6537184 -2.563817  0.01680535
#> 5 27.4463185 64.17455 25.072882 0.9142508 23.309179 26.89778573
#> 6 22.0610004 51.90959 20.116260 0.6241193 18.873557 21.32615735
head(residuals(demo_fit))  # Residuals for each row of demo_fit$data
#>     response     time  residuals     error        Q2.5     Q97.5
#> 1 -3.0084198 91.48060  0.6340487 0.7403902 -0.80977256  2.021618
#> 2 -7.8768640 93.70754 -3.6011982 0.8303914 -5.21134983 -2.035040
#> 3 16.3029101 28.61395  5.3454359 1.0249532  3.38967601  7.247332
#> 4 -0.0373553 83.04476  1.2065039 0.6537184 -0.05416065  2.526462
#> 5 27.4463185 64.17455  2.3734361 0.9142508  0.54853280  4.137139
#> 6 22.0610004 51.90959  1.9447407 0.6241193  0.73484310  3.187443
log_lik(demo_fit)[1:3, 1:3]  # Log-likelihood at each demo_fit$data
#>              1         2         3
#> [1,] -2.338686 -3.019339 -2.902410
#> [2,] -2.264241 -2.855583 -3.105783
#> [3,] -2.268495 -2.812040 -3.040566

# All of the above take a range of arguments. E.g.,:
# \donttest{
head(predict(demo_fit))  # Pointwise posterior predictive
#>     response     time   predict    error       Q2.5     Q97.5
#> 1 -3.0084198 91.48060 -3.724386 3.740772 -11.016678  3.892476
#> 2 -7.8768640 93.70754 -4.216435 3.847017 -11.875743  3.323773
#> 3 16.3029101 28.61395 10.834875 3.847943   3.431532 18.314491
#> 4 -0.0373553 83.04476 -1.150261 3.688370  -8.375954  6.330058
#> 5 27.4463185 64.17455 25.111757 3.962978  17.459322 32.764153
#> 6 22.0610004 51.90959 20.133693 3.699329  12.942790 27.580408
head(predict(demo_fit, probs = c(0.1, 0.5, 0.9)))  # Median and 80% posterior predictive interval.
#>     response     time   predict    error       Q10       Q50        Q90
#> 1 -3.0084198 91.48060 -3.652897 3.822827 -8.555410 -3.630860  1.1589727
#> 2 -7.8768640 93.70754 -4.283365 3.735920 -9.080068 -4.272670  0.3787073
#> 3 16.3029101 28.61395 10.896005 3.837307  6.067081 10.712469 15.8965915
#> 4 -0.0373553 83.04476 -1.263612 3.805213 -6.028461 -1.173006  3.5101741
#> 5 27.4463185 64.17455 25.138853 3.757026 20.383688 25.218862 29.9478814
#> 6 22.0610004 51.90959 20.182447 3.695304 15.449103 20.186691 24.9891866
head(predict(demo_fit, prior = TRUE))  # Prior predictive
#>     response     time   predict    error      Q2.5    Q97.5
#> 1 -3.0084198 91.48060  9.323891 38.63634 -59.73841 82.42399
#> 2 -7.8768640 93.70754  8.560608 37.72807 -67.28144 76.15495
#> 3 16.3029101 28.61395  9.343285 28.78061 -38.81967 64.20029
#> 4 -0.0373553 83.04476  9.159680 35.06269 -58.40516 73.94821
#> 5 27.4463185 64.17455 10.088712 31.33772 -51.15107 72.88376
#> 6 22.0610004 51.90959  8.922906 30.10071 -48.79053 65.31738
head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#> # A tibble: 6 × 14
#>   .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3 time_3
#>    <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>  <dbl>
#> 1      1          1     1  23.6  70.1        10.1  0.360        4.73 -0.331
#> 2      1          1     1  23.6  70.1        10.1  0.360        4.73 -0.331
#> 3      1          1     1  23.6  70.1        10.1  0.360        4.73 -0.331
#> 4      1          1     1  23.6  70.1        10.1  0.360        4.73 -0.331
#> 5      1          1     1  23.6  70.1        10.1  0.360        4.73 -0.331
#> 6      1          1     1  23.6  70.1        10.1  0.360        4.73 -0.331
#> # ℹ 5 more variables: sigma_1 <dbl>, response <dbl>, time <dbl>,
#> #   data_row <int>, .epred <dbl>
head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#>     response     time   fitted     error     Q2.5   Q97.5
#> 1 -3.0084198 91.48060 3.667508 0.2718435 3.174695 4.22363
#> 2 -7.8768640 93.70754 3.667508 0.2718435 3.174695 4.22363
#> 3 16.3029101 28.61395 3.667508 0.2718435 3.174695 4.22363
#> 4 -0.0373553 83.04476 3.667508 0.2718435 3.174695 4.22363
#> 5 27.4463185 64.17455 3.667508 0.2718435 3.174695 4.22363
#> 6 22.0610004 51.90959 3.667508 0.2718435 3.174695 4.22363

# Evaluate at novel data
novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
head(predict(demo_fit, newdata = novel_data, probs = c(0.025, 0.5, 0.975)))
#>   time    predict     error       Q2.5        Q50     Q97.5
#> 1   -5   8.979614  3.834472   1.552215   8.949600  16.45513
#> 2   20   9.312321  3.714050   2.115774   9.266639  16.77730
#> 3  300 -63.049547 15.164624 -90.001744 -63.814536 -32.14780

# Work with missing responses
missing_fit = mcp_example("missing", plot = FALSE)
#> NA values detected in 'y'. JAGS will treat them as latent responses and impute them during sampling.
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 90
#>    Unobserved stochastic nodes: 16
#>    Total graph size: 1635
#> 
#> Initializing model
#> 
#> Finished sampling in 3 seconds
fitted(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Expected responses for missing y
#>    y  x condition   fitted     error     Q2.5    Q97.5
#> 1 NA  8         B 35.59660 1.0822734 33.48860 37.77410
#> 2 NA 19         A 15.67273 0.8439373 14.03083 17.34643
#> 3 NA 27         A 17.18509 0.7712880 15.68295 18.69161
#> 4 NA 28         B 39.37751 0.7529116 37.89127 40.87649
#> 5 NA 29         A 17.56318 0.7702221 16.06121 19.06213
#> 6 NA 30         B 39.75560 0.7532309 38.27336 41.24870
fitted(missing_fit, summary = FALSE) |> dplyr::filter(is.na(y)) |> head()  # Same, but draws
#> # A tibble: 6 × 14
#>   .chain .iteration .draw  cp_1 Intercept_1   x_1 conditionB_1    x_2 sigma_1
#>    <int>      <int> <int> <dbl>       <dbl> <dbl>        <dbl>  <dbl>   <dbl>
#> 1      1          1     1  60.0        12.6 0.180         22.2 -0.410    4.54
#> 2      1          1     1  60.0        12.6 0.180         22.2 -0.410    4.54
#> 3      1          1     1  60.0        12.6 0.180         22.2 -0.410    4.54
#> 4      1          1     1  60.0        12.6 0.180         22.2 -0.410    4.54
#> 5      1          1     1  60.0        12.6 0.180         22.2 -0.410    4.54
#> 6      1          1     1  60.0        12.6 0.180         22.2 -0.410    4.54
#> # ℹ 5 more variables: y <dbl>, x <int>, condition <fct>, data_row <int>,
#> #   .epred <dbl>
predict(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Posterior predictive for missing y
#>    y  x condition  predict    error      Q2.5    Q97.5
#> 1 NA  8         B 35.56253 4.419118 26.817180 44.15292
#> 2 NA 19         A 15.73385 4.291075  7.308634 24.25204
#> 3 NA 27         A 17.18995 4.319638  8.640670 25.73600
#> 4 NA 28         B 39.35814 4.360269 30.829841 47.93401
#> 5 NA 29         A 17.53493 4.342543  8.969956 25.91751
#> 6 NA 30         B 39.72457 4.395872 31.119327 48.44444
# }
```
