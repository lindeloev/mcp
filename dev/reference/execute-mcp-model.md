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
  width = 0.95,
  rate = FALSE,
  prior = FALSE,
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  probs = lifecycle::deprecated(),
  conditional = TRUE,
  ...
)

# S3 method for class 'mcpfit'
fitted(
  object,
  newdata = NULL,
  summary = TRUE,
  width = 0.95,
  rate = FALSE,
  prior = FALSE,
  dpar = "epred",
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = "response",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  probs = lifecycle::deprecated(),
  ...
)

log_lik(object, ...)

# S3 method for class 'mcpfit'
log_lik(
  object,
  newdata = NULL,
  summary = FALSE,
  width = 0.95,
  prior = FALSE,
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "matrix",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  probs = lifecycle::deprecated(),
  ...
)

# S3 method for class 'mcpfit'
residuals(
  object,
  newdata = NULL,
  summary = TRUE,
  width = 0.95,
  prior = FALSE,
  group_effects = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  nsamples = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  probs = lifecycle::deprecated(),
  ...
)
```

## Arguments

- object:

  An `mcpfit` object.

- newdata:

  A `tibble` or a `data.frame` containing predictors in the model.

  - If `NULL` (default), the original data is used.

  - Grouping factors must only contain levels from the fitted data. Use
    `group_effects = FALSE` to predict without group-level effects,
    e.g., for new groups.

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

- width:

  Width of the central interval (strictly between 0 and 1), returned in
  the `.lower` and `.upper` columns. `FALSE` omits the interval. Only in
  effect when `summary == TRUE`. To get custom quantiles, summarise the
  draws from `summary = FALSE`, e.g., using
  [`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html)
  (see examples).

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

- group_effects:

  Group-level effects to include. One of:

  - `TRUE` All group-level deviations.

  - `FALSE` or `NA` No group-level deviations, i.e., population-level
    predictions. Like `re_formula = NA` in `{brms}`.

  - `"cp"` or `"predictor"`: All group-level deviations belonging to
    that part of the model.

  - Character vector of group-level parameter names, e.g.,
    `c("cp_1_id", "x_2_id")`. `*` is a wildcard, e.g., `"state*_1_id"`
    selects all levels of a factor and `"*_id"` selects all effects by
    `id`. See names using `ranef(fit)`. Can be combined, e.g.,
    `c("cp", "sigma_*")`.

- arma:

  Whether to include AR and MA effects.

  - `TRUE` Compute the GARMA residual recurrence. Requires the response
    variable in `newdata`.

  - `FALSE` Disregard AR and MA effects. For `family = gaussian()`,
    `predict()` uses only `sigma` for residuals. For posterior
    evaluation of the original data, retained JAGS imputations supply
    missing GARMA histories. In models with group-level effects, this
    currently requires all such effects (`group_effects = TRUE`).

- ndraws:

  Integer or `NULL`. Number of posterior draws to return/summarise.
  `NULL` means "all". More draws trade speed for accuracy.

- draws_format:

  One of "tidy" or "matrix". Controls the output format when
  `summary == FALSE`.

- nsamples:

  Deprecated. Use `ndraws` instead.

- samples_format:

  Deprecated. Use `draws_format` instead.

- varying:

  Deprecated. Use `group_effects` instead.

- probs:

  Deprecated. Use `width` instead.

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

- If `summary = TRUE`: A tibble with the columns of `newdata` (or
  `fit$data` if `newdata` is `NULL`), followed by summaries for each
  row:

  - `.mean`: The mean of the draws.

  - `.sd`: The SD of the draws. With posterior draws (the default), this
    is the posterior predictive SD for `type = "predict"` and the
    posterior SD of the evaluated quantity otherwise.

  - `.lower` and `.upper`: The bounds of the central interval given in
    `width`.

  With `prior = TRUE`, these are the analogous prior summaries. mcp's
  own columns are dot-prefixed to avoid collisions with data columns.

- If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes`
  `tibble` with all the posterior draws (`Nd`) evaluated at each row in
  `newdata` (`Nn`), i.e., with `Nd x Nn` rows. If there are group-level
  effects, the returned data is expanded with the relevant levels for
  each row.

  The return columns are:

  - Predictors from `newdata`, plus its response column when supplied.

  - Draw descriptors: ".chain", ".iteration", ".draw" (see the
    `posterior` and `tidybayes` packages), and `.data_row`, the row
    number in the evaluated `newdata`.

  - Draw values: one column for each parameter in the model.

  - The estimate (matching tidybayes conventions): `.epred` for expected
    responses, `.linpred` with `scale = "linear"`, the `.dpar` for other
    dpars (e.g., `.mu` or `.sigma`), `.prediction`, `.residual`, or
    `.loglik`.

- If `summary = FALSE` and `draws_format = "matrix"`: An `N_draws` X
  `nrows(newdata)` matrix with fitted/predicted values (depending on
  `type`). This format is used by `brms` and it's useful as `yrep` in
  `bayesplot::ppc_*` functions.

## Details

`fitted()` and
[`posterior_epred()`](https://lindeloev.github.io/mcp/dev/reference/posterior_epred.mcpfit.md)
evaluate the same expected responses; `predict()` and
[`posterior_predict()`](https://lindeloev.github.io/mcp/dev/reference/posterior_epred.mcpfit.md)
evaluate the same response distributions. The `posterior_*()` methods
return draws-by-observation matrices, while `fitted()` and `predict()`
summarise by default and also offer tidy draws. For binomial models, the
default response scale is counts. Use `rate = TRUE` for proportions or
`fitted(..., dpar = "mu")` for success probabilities. During migration
from v0.3.4, an omitted `rate` warns once per session per function when
counts differ from proportions. Explicit `rate` settings do not warn.

`residuals(fit)` is equivalent to
`fit$data[[mcp_columns(fit)$response]] - fitted(fit, ...)` (or
`newdata[[mcp_columns(fit)$response]] - fitted(fit, ...)`), but with
fixed arguments for `fitted`: `rate = FALSE, dpar = 'epred'`.

`log_lik()` defaults to an unsummarised draws-by-observation matrix, as
used by `loo` and other posterior workflows. Non-default `group_effects`
and `arma` settings evaluate conditional or counterfactual
log-likelihoods (e.g., omitting random effects or serial correlation);
they cannot be used in
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
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  12.9 1.02   11.0    14.9
#> 2     6.35  83.5  11.4 0.742   9.95   12.9
#> 3    28.0   60.2  25.4 0.863  23.8    27.1
#> 4    15.7   74.7  13.2 1.12   11.0    15.5
#> 5    16.2   85.9  10.9 0.736   9.45   12.3
#> 6    22.3   40.1  14.5 0.656  13.2    15.7
head(residuals(demo_fit))  # Residuals for each row of demo_fit$data
#> # A tibble: 6 × 6
#>   response  time  .mean   .sd .lower .upper
#>      <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3 -0.637 1.02  -2.68    1.28
#> 2     6.35  83.5 -5.04  0.742 -6.51   -3.60
#> 3    28.0   60.2  2.60  0.863  0.938   4.29
#> 4    15.7   74.7  2.48  1.12   0.215   4.65
#> 5    16.2   85.9  5.31  0.736  3.93    6.76
#> 6    22.3   40.1  7.82  0.656  6.65    9.16
log_lik(demo_fit)[1:3, 1:3]  # Log-likelihood at each demo_fit$data
#>              1         2         3
#> [1,] -2.508965 -3.555050 -2.610885
#> [2,] -2.313553 -3.096477 -2.895883
#> [3,] -2.520140 -3.468677 -2.504247

# All of the above take a range of arguments. E.g.,:
# \donttest{
head(predict(demo_fit))  # Pointwise posterior predictive
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  12.7  3.90   4.99   20.7
#> 2     6.35  83.5  11.4  3.91   3.63   19.1
#> 3    28.0   60.2  25.3  4.04  17.6    33.3
#> 4    15.7   74.7  13.3  4.01   5.27   21.1
#> 5    16.2   85.9  10.8  4.11   3.14   18.6
#> 6    22.3   40.1  14.6  4.03   6.77   22.2
head(predict(demo_fit, width = 0.8))  # 80% posterior predictive interval
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  12.8  4.16   7.75   18.0
#> 2     6.35  83.5  11.5  4.15   6.35   16.4
#> 3    28.0   60.2  25.3  3.95  20.4    30.5
#> 4    15.7   74.7  13.3  4.08   8.05   18.4
#> 5    16.2   85.9  11.0  4.17   5.86   15.9
#> 6    22.3   40.1  14.4  3.77   9.48   19.5
head(predict(demo_fit, prior = TRUE))  # Prior predictive
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  14.0  28.1  -47.3   70.7
#> 2     6.35  83.5  14.4  30.4  -50.8   75.1
#> 3    28.0   60.2  13.7  26.2  -39.8   66.2
#> 4    15.7   74.7  13.9  28.4  -46.3   69.6
#> 5    16.2   85.9  14.4  29.5  -52.1   76.2
#> 6    22.3   40.1  14.2  22.5  -30.0   58.3
head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#> # A tibble: 6 × 14
#>   .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3 time_3
#>    <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>  <dbl>
#> 1      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 2      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 3      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 4      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 5      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> 6      1          1     1  28.6  71.6        8.78  0.508        16.8 -0.372
#> # ℹ 5 more variables: sigma_1 <dbl>, response <dbl>, time <dbl>,
#> #   .data_row <int>, .epred <dbl>
head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  3.87 0.280   3.37   4.44
#> 2     6.35  83.5  3.87 0.280   3.37   4.44
#> 3    28.0   60.2  3.87 0.280   3.37   4.44
#> 4    15.7   74.7  3.87 0.280   3.37   4.44
#> 5    16.2   85.9  3.87 0.280   3.37   4.44
#> 6    22.3   40.1  3.87 0.280   3.37   4.44

# Custom quantiles from the draws, here for fitted().
# Note that for predict(), this adds more Monte Carlo error than `width`.
fitted(demo_fit, summary = FALSE) |>
  dplyr::group_by(.data_row, time) |>
  dplyr::reframe(quantile = c(0.1, 0.5, 0.9), .epred = quantile(.epred, quantile)) |>
  head()
#> # A tibble: 6 × 4
#>   .data_row  time quantile .epred
#>       <int> <dbl>    <dbl>  <dbl>
#> 1         1  76.3      0.1   11.6
#> 2         1  76.3      0.5   12.8
#> 3         1  76.3      0.9   14.2
#> 4         2  83.5      0.1   10.4
#> 5         2  83.5      0.5   11.4
#> 6         2  83.5      0.9   12.3

# Evaluate at novel data
novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
head(predict(demo_fit, newdata = novel_data))
#> # A tibble: 3 × 5
#>    time  .mean   .sd .lower .upper
#>   <dbl>  <dbl> <dbl>  <dbl>  <dbl>
#> 1    -5   9.78  3.89   2.22  17.7 
#> 2    20  10.1   3.82   2.22  17.7 
#> 3   300 -33.5  17.9  -67.8    2.91

# Work with missing responses
missing_fit = mcp_example("missing", plot = FALSE)
#> NA values detected in 'y'. JAGS will treat them as latent responses and impute them during sampling.
fitted(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Expected responses for missing y
#> # A tibble: 6 × 7
#>       y     x state .mean   .sd .lower .upper
#>   <dbl> <int> <fct> <dbl> <dbl>  <dbl>  <dbl>
#> 1    NA     8 B      35.4 1.09    33.2   37.5
#> 2    NA    19 A      15.5 0.831   13.9   17.1
#> 3    NA    27 A      17.1 0.758   15.6   18.6
#> 4    NA    28 B      39.4 0.758   37.9   40.9
#> 5    NA    29 A      17.5 0.758   16.0   19.0
#> 6    NA    30 B      39.8 0.759   38.3   41.3
fitted(missing_fit, summary = FALSE) |> dplyr::filter(is.na(y)) |> head()  # Same, but draws
#> # A tibble: 6 × 14
#>   .chain .iteration .draw  cp_1 Intercept_1   x_1 stateB_1    x_2 sigma_1     y
#>    <int>      <int> <int> <dbl>       <dbl> <dbl>    <dbl>  <dbl>   <dbl> <dbl>
#> 1      1          1     1  57.9        11.1 0.192     22.0 -0.360    4.13    NA
#> 2      1          1     1  57.9        11.1 0.192     22.0 -0.360    4.13    NA
#> 3      1          1     1  57.9        11.1 0.192     22.0 -0.360    4.13    NA
#> 4      1          1     1  57.9        11.1 0.192     22.0 -0.360    4.13    NA
#> 5      1          1     1  57.9        11.1 0.192     22.0 -0.360    4.13    NA
#> 6      1          1     1  57.9        11.1 0.192     22.0 -0.360    4.13    NA
#> # ℹ 4 more variables: x <int>, state <fct>, .data_row <int>, .epred <dbl>
predict(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Posterior predictive for missing y
#> # A tibble: 6 × 7
#>       y     x state .mean   .sd .lower .upper
#>   <dbl> <int> <fct> <dbl> <dbl>  <dbl>  <dbl>
#> 1    NA     8 B      35.4  4.41  26.7    44.1
#> 2    NA    19 A      15.5  4.39   6.95   24.1
#> 3    NA    27 A      17.1  4.39   8.56   25.7
#> 4    NA    28 B      39.4  4.32  30.8    47.9
#> 5    NA    29 A      17.5  4.36   8.96   26.1
#> 6    NA    30 B      39.7  4.38  31.2    48.3
# }
```
