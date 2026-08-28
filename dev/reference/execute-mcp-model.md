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
#>    response     time    fitted     error      Q2.5    Q97.5
#> 1 32.842651 68.35820 30.298834 1.0124885 28.321430 32.25606
#> 2 -1.160003 87.29038 -3.289723 0.7599234 -4.782638 -1.80833
#> 3 27.564248 69.01173 30.646721 1.0409779 28.633151 32.68228
#> 4 10.062971 11.59361 10.299931 0.7167564  8.874782 11.66981
#> 5 14.056859 19.50091 10.302634 0.7135642  8.942350 11.66981
#> 6 18.292640 46.12009 18.460976 0.9281008 16.431898 20.12248
head(residuals(demo_fit))  # Residuals for each row of demo_fit$data
#>    response     time  residuals     error       Q2.5     Q97.5
#> 1 32.842651 68.35820  2.5438170 1.0124885  0.5865913  4.521221
#> 2 -1.160003 87.29038  2.1297203 0.7599234  0.6483267  3.622635
#> 3 27.564248 69.01173 -3.0824732 1.0409779 -5.1180353 -1.068903
#> 4 10.062971 11.59361 -0.2369604 0.7167564 -1.6068402  1.188188
#> 5 14.056859 19.50091  3.7542243 0.7135642  2.3870478  5.114508
#> 6 18.292640 46.12009 -0.1683359 0.9281008 -1.8298407  1.860741
log_lik(demo_fit)[1:3, 1:3]  # Log-likelihood at each demo_fit$data
#>              1         2         3
#> [1,] -3.167203 -2.408024 -2.306866
#> [2,] -2.956388 -2.324250 -2.270318
#> [3,] -2.747960 -2.362946 -2.384176

# All of the above take a range of arguments. E.g.,:
# \donttest{
head(predict(demo_fit))  # Pointwise posterior predictive
#>    response     time   predict    error       Q2.5     Q97.5
#> 1 32.842651 68.35820 30.207150 4.124578  22.144905 38.439454
#> 2 -1.160003 87.29038 -3.210645 4.151474 -11.321679  4.760017
#> 3 27.564248 69.01173 30.518943 4.165452  22.479043 38.801177
#> 4 10.062971 11.59361 10.414788 4.078595   2.277202 18.328893
#> 5 14.056859 19.50091 10.347444 4.292330   2.282246 18.331710
#> 6 18.292640 46.12009 18.469592 4.078643  10.340106 26.560291
head(predict(demo_fit, probs = c(0.1, 0.5, 0.9)))  # Median and 80% posterior predictive interval.
#>    response     time   predict    error       Q10       Q50       Q90
#> 1 32.842651 68.35820 30.292511 4.177484 25.003207 30.301147 35.591733
#> 2 -1.160003 87.29038 -3.293718 4.043098 -8.512640 -3.292723  1.936534
#> 3 27.564248 69.01173 30.589537 4.167861 25.342076 30.649014 35.948693
#> 4 10.062971 11.59361 10.285350 4.140679  5.086742 10.298795 15.515312
#> 5 14.056859 19.50091 10.382461 4.005179  5.090210 10.301096 15.517589
#> 6 18.292640 46.12009 18.539889 4.098412 13.190974 18.464812 23.725677
head(predict(demo_fit, prior = TRUE))  # Prior predictive
#>    response     time   predict    error      Q2.5    Q97.5
#> 1 32.842651 68.35820 11.171082 32.66496 -49.74378 70.57951
#> 2 -1.160003 87.29038  9.864536 32.13784 -52.55783 71.93426
#> 3 27.564248 69.01173 10.111668 30.81516 -49.85891 70.60924
#> 4 10.062971 11.59361 10.282018 33.31249 -49.26736 68.24167
#> 5 14.056859 19.50091 11.242114 31.90038 -49.54965 68.65070
#> 6 18.292640 46.12009  9.654039 31.85612 -49.21910 70.55472
head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#> # A tibble: 6 × 14
#>   .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3 time_3
#>    <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>  <dbl>
#> 1      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 2      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 3      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 4      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 5      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> 6      1          1     1  19.8  70.2        8.44  0.395        2.11 -0.298
#> # ℹ 5 more variables: sigma_1 <dbl>, response <dbl>, time <dbl>,
#> #   data_row <int>, .epred <dbl>
head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#>    response     time   fitted    error     Q2.5    Q97.5
#> 1 32.842651 68.35820 4.010101 0.311438 3.448848 4.656798
#> 2 -1.160003 87.29038 4.010101 0.311438 3.448848 4.656798
#> 3 27.564248 69.01173 4.010101 0.311438 3.448848 4.656798
#> 4 10.062971 11.59361 4.010101 0.311438 3.448848 4.656798
#> 5 14.056859 19.50091 4.010101 0.311438 3.448848 4.656798
#> 6 18.292640 46.12009 4.010101 0.311438 3.448848 4.656798

# Evaluate at novel data
novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
head(predict(demo_fit, newdata = novel_data, probs = c(0.025, 0.5, 0.975)))
#>   time   predict     error       Q2.5       Q50     Q97.5
#> 1   -5  10.19913  4.178428   2.277202  10.29880  18.32889
#> 2   20  10.36892  3.985927   2.284798  10.30251  18.33282
#> 3  300 -50.84247 20.327313 -90.157846 -50.79259 -12.33334

# Work with missing responses
missing_fit = mcp_example("missing", plot = FALSE)
#> NA values detected in 'y'. JAGS will treat them as latent responses and impute them during sampling.
fitted(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Expected responses for missing y
#>    y  x condition   fitted     error     Q2.5    Q97.5
#> 1 NA  8         B 35.59410 1.0808576 33.50826 37.76016
#> 2 NA 19         A 15.67323 0.8369118 14.05072 17.33612
#> 3 NA 27         A 17.18394 0.7645671 15.69881 18.68233
#> 4 NA 28         B 39.37086 0.7545884 37.87833 40.86702
#> 5 NA 29         A 17.56162 0.7634425 16.07303 19.04347
#> 6 NA 30         B 39.74854 0.7546261 38.25985 41.23888
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
#> 1 NA  8         B 35.56012 4.428752 26.945587 44.26715
#> 2 NA 19         A 15.72108 4.293560  7.125812 24.23679
#> 3 NA 27         A 17.18509 4.292231  8.656528 25.71505
#> 4 NA 28         B 39.36003 4.364408 30.842388 47.89444
#> 5 NA 29         A 17.53601 4.340269  9.033103 26.09081
#> 6 NA 30         B 39.72210 4.371321 31.218725 48.27069
# }
```
