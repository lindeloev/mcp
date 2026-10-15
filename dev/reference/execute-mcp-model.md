# Fitted and predicted values of `mcp` models fits

Evaluate the model on data, either summarised (per data-row) or per
draw. You can use draws from the prior (`prior = TRUE`), select a
distributional parameter with `dpar`, and choose the response or
linear-predictor scale with `scale` where applicable. See [Fits and
predictions](https://lindeloev.github.io/mcp/articles/predict.html) for
worked examples, including forecasts with future change points.

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

  Logical. Summarise the draws for each row in `newdata` (`TRUE`) or
  return the draws (`FALSE`)?

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
    `c("cp_1_id", "x_2_id")`. **\[experimental\]** `*` is a wildcard,
    e.g., `"state*_1_id"` selects all levels of a factor and `"*_id"`
    selects all effects by `id`. See names using `ranef(fit)`. Can be
    combined, e.g., `c("cp", "sigma_*")`.

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

  What distributional parameter to evaluate. E.g.,

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
    linear trends are modeled.

## Value

- If `summary = TRUE`: A tibble with the columns of `newdata` (or
  `fit$data` if `newdata` is `NULL`), followed by summaries for each
  row:

  - `.mean`: The mean of the draws.

  - `.sd`: The SD of the draws. With posterior draws (the default), this
    is the posterior predictive SD for `predict()` and the posterior SD
    of the evaluated quantity otherwise.

  - `.lower` and `.upper`: The lower and upper ends of the central
    interval given in `width`.

  With `prior = TRUE`, these are the analogous prior summaries. mcp's
  own columns are dot-prefixed to avoid collisions with data columns.

- If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes`
  `tibble` with all the posterior draws (`Nd`) evaluated at each row in
  `newdata` (`Nn`), i.e., with `Nd x Nn` rows, ordered by `.data_row`
  and then `.draw`. If there are group-level effects, the returned data
  is expanded with the relevant levels for each row.

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
  `nrows(newdata)` matrix with fitted/predicted values (depending on the
  function). This format is used by `brms` and it's useful as `yrep` in
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
log-likelihoods (e.g., without group-level effects or serial
correlation); they cannot be used in
[`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
or
[`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
because estimating information criteria for reduced models requires
refitting.

Missing responses in the original data remain missing in the response
column. `fitted()` returns their expected responses, while `predict()`
generates fresh posterior predictive response draws (including at
missing rows). In GARMA models, retained JAGS imputations supply the
history used for later fitted and predicted rows. [Read
more](https://lindeloev.github.io/mcp/articles/missing.html).

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
#> 1    12.2   76.3  12.9 1.01   10.9    14.8
#> 2     6.35  83.5  11.4 0.717   9.96   12.8
#> 3    28.0   60.2  25.4 0.873  23.7    27.3
#> 4    15.7   74.7  13.2 1.11   11.0    15.4
#> 5    16.2   85.9  10.9 0.711   9.51   12.3
#> 6    22.3   40.1  14.5 0.621  13.2    15.6
head(residuals(demo_fit))  # Residuals for each row of demo_fit$data
#> # A tibble: 6 × 6
#>   response  time  .mean   .sd .lower .upper
#>      <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3 -0.645 1.01  -2.54    1.38
#> 2     6.35  83.5 -5.05  0.717 -6.43   -3.61
#> 3    28.0   60.2  2.61  0.873  0.784   4.33
#> 4    15.7   74.7  2.47  1.11   0.327   4.71
#> 5    16.2   85.9  5.30  0.711  3.89    6.70
#> 6    22.3   40.1  7.83  0.621  6.68    9.14
log_lik(demo_fit)[1:3, 1:3]  # Log-likelihood at each demo_fit$data
#>              1         2         3
#> [1,] -2.193497 -3.449599 -2.445390
#> [2,] -2.317434 -3.373387 -2.342121
#> [3,] -2.311202 -3.090324 -2.419474

# All of the above take a range of arguments. E.g.,:
# \donttest{
head(predict(demo_fit))  # Pointwise posterior predictive
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  12.7  3.93   5.01   20.7
#> 2     6.35  83.5  11.5  3.85   3.65   19.1
#> 3    28.0   60.2  25.3  4.06  17.6    33.2
#> 4    15.7   74.7  13.3  3.92   5.29   21.1
#> 5    16.2   85.9  10.9  4.14   3.17   18.6
#> 6    22.3   40.1  14.6  4.03   6.79   22.2
head(predict(demo_fit, width = 0.8))  # 80% posterior predictive interval
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  12.8  4.03   7.77   18.0
#> 2     6.35  83.5  11.5  4.12   6.37   16.4
#> 3    28.0   60.2  25.3  3.97  20.4    30.5
#> 4    15.7   74.7  13.3  4.02   8.07   18.4
#> 5    16.2   85.9  11.1  4.09   5.88   15.9
#> 6    22.3   40.1  14.4  3.74   9.49   19.5
head(predict(demo_fit, prior = TRUE))  # Prior predictive
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  13.5  29.4  -42.7   76.7
#> 2     6.35  83.5  12.9  31.0  -46.6   81.2
#> 3    28.0   60.2  13.6  26.8  -37.9   68.0
#> 4    15.7   74.7  13.4  28.5  -41.0   76.2
#> 5    16.2   85.9  12.8  31.1  -48.3   80.8
#> 6    22.3   40.1  14.2  23.5  -31.2   55.4
head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#> # A tibble: 6 × 14
#>   response  time .data_row .chain .iteration .draw  cp_1  cp_2 Intercept_1
#>      <dbl> <dbl>     <int>  <int>      <int> <int> <dbl> <dbl>       <dbl>
#> 1     12.2  76.3         1      1          1     1  31.9  71.0        8.95
#> 2     12.2  76.3         1      1          2     2  31.7  71.4        9.87
#> 3     12.2  76.3         1      1          3     3  32.1  70.5       11.4 
#> 4     12.2  76.3         1      1          4     4  33.3  72.8       10.1 
#> 5     12.2  76.3         1      1          5     5  31.7  69.9        9.97
#> 6     12.2  76.3         1      1          6     6  33.3  69.9       10.0 
#> # ℹ 5 more variables: time_2 <dbl>, Intercept_3 <dbl>, time_3 <dbl>,
#> #   sigma_1 <dbl>, .epred <dbl>
head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  3.86 0.293   3.34   4.54
#> 2     6.35  83.5  3.86 0.293   3.34   4.54
#> 3    28.0   60.2  3.86 0.293   3.34   4.54
#> 4    15.7   74.7  3.86 0.293   3.34   4.54
#> 5    16.2   85.9  3.86 0.293   3.34   4.54
#> 6    22.3   40.1  3.86 0.293   3.34   4.54

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
#> 4         2  83.5      0.1   10.5
#> 5         2  83.5      0.5   11.4
#> 6         2  83.5      0.9   12.3

# Evaluate at novel data
novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
head(predict(demo_fit, newdata = novel_data))
#> # A tibble: 3 × 5
#>    time  .mean   .sd .lower .upper
#>   <dbl>  <dbl> <dbl>  <dbl>  <dbl>
#> 1    -5   9.80  3.89   2.25  17.7 
#> 2    20  10.1   3.81   2.25  17.7 
#> 3   300 -33.4  18.4  -69.9    2.39

# Work with missing responses
missing_fit = mcp_example("missing", plot = FALSE)
#> NA values detected in 'y'. JAGS will treat them as latent responses and impute them during sampling.
fitted(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Expected responses for missing y
#> # A tibble: 6 × 7
#>       y     x state .mean   .sd .lower .upper
#>   <dbl> <int> <fct> <dbl> <dbl>  <dbl>  <dbl>
#> 1    NA     8 B      35.4 1.09    33.2   37.5
#> 2    NA    19 A      15.5 0.833   13.8   17.1
#> 3    NA    27 A      17.1 0.759   15.6   18.6
#> 4    NA    28 B      39.4 0.758   37.9   40.9
#> 5    NA    29 A      17.5 0.759   16.0   19.0
#> 6    NA    30 B      39.8 0.760   38.3   41.3
fitted(missing_fit, summary = FALSE) |> dplyr::filter(is.na(y)) |> head()  # Same, but draws
#> # A tibble: 6 × 14
#>       y     x state .data_row .chain .iteration .draw  cp_1 Intercept_1   x_1
#>   <dbl> <int> <fct>     <int>  <int>      <int> <int> <dbl>       <dbl> <dbl>
#> 1    NA     8 B             8      1          1     1  54.1        12.9 0.202
#> 2    NA     8 B             8      1          2     2  54.1        11.8 0.204
#> 3    NA     8 B             8      1          3     3  57.7        11.7 0.183
#> 4    NA     8 B             8      1          4     4  61.1        13.0 0.177
#> 5    NA     8 B             8      1          5     5  57.9        13.7 0.160
#> 6    NA     8 B             8      1          6     6  58.0        11.0 0.242
#> # ℹ 4 more variables: stateB_1 <dbl>, x_2 <dbl>, sigma_1 <dbl>, .epred <dbl>
predict(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Posterior predictive for missing y
#> # A tibble: 6 × 7
#>       y     x state .mean   .sd .lower .upper
#>   <dbl> <int> <fct> <dbl> <dbl>  <dbl>  <dbl>
#> 1    NA     8 B      35.4  4.40  26.7    44.1
#> 2    NA    19 A      15.5  4.36   6.93   24.1
#> 3    NA    27 A      17.0  4.38   8.55   25.6
#> 4    NA    28 B      39.4  4.31  30.8    47.9
#> 5    NA    29 A      17.5  4.34   8.95   26.0
#> 6    NA    30 B      39.7  4.38  31.2    48.3
# }
```
