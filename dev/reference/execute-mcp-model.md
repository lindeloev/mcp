# Fitted and Predicted Values of `mcpfit` Objects

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
    `arma = "replicate"`, `predict()` and
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
  itself (e.g., success probability) and are unaffected by `rate`. Other
  families do not support `rate = TRUE`. For rates per unit of exposure
  in models with `offset(log(exposure))`, set `exposure = 1` in
  `newdata`.

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

  How to include AR and MA effects.

  - `TRUE`: Condition on the observed response history, i.e.,
    one-step-ahead evaluation of the GARMA residual recurrence. Requires
    the response variable in `newdata`.

  - `FALSE`: Disregard AR and MA effects. For `family = gaussian()`,
    `predict()` uses only `sigma` for residuals.

  - `"replicate"`: Generate fresh response series recursively instead of
    conditioning on the observed responses. Only for predictions:
    `predict()`,
    [`posterior_predict()`](https://lindeloev.github.io/mcp/dev/reference/posterior_epred.mcpfit.md),
    and
    [`pp_check()`](https://lindeloev.github.io/mcp/dev/reference/pp_check.md)
    (where it is the default). Applies equally to prior and posterior
    draws.

  For posterior evaluation of the original data, retained JAGS
  imputations supply missing GARMA histories. In models with group-level
  effects, this currently requires all such effects
  (`group_effects = TRUE`).

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

- ...:

  Must be empty. Reserved for future use.

- dpar:

  What distributional parameter to evaluate. E.g.,

  - `"epred"` (default): Expected response from the full model (or
    `NULL` for compatibility with brms etc.).

  - `"mu"`: The conditional mean (or success probability per trial for
    binomial/bernoulli models), on the link or response scale.

  - `"sigma"`: The standard deviation of the residuals. With
    [`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()`, it is the SD of
    the innovations, i.e., the part of the residuals not explained by
    the AR/MA terms.

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

  For `predict()` with built-in families, these are the exact moments
  and quantiles of the posterior predictive distribution given the
  evaluated draws (the mixture over draws), so they do not vary with the
  seed. Two exceptions: `ndraws` selects draws at random, and
  `arma = "replicate"` conditions on randomly generated response
  histories, i.e., a Monte Carlo approximation.

  With `prior = TRUE`, these are the analogous prior summaries. mcp's
  own columns are dot-prefixed to avoid collisions with data columns.

- If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes`
  `tibble` with all the posterior draws (`Nd`) evaluated at each row in
  `newdata` (`Nn`), i.e., with `Nd x Nn` rows, ordered by `.row` and
  then `.draw`. If there are group-level effects, the returned data is
  expanded with the relevant levels for each row.

  The return columns are:

  - Predictors from `newdata`, plus its response column when supplied.

  - Draw descriptors: ".chain", ".iteration", ".draw" (see the
    `posterior` and `tidybayes` packages), and `.row`, the row number in
    the evaluated `newdata`.

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
missing rows). In AR/MA models, retained JAGS imputations supply the
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
#> 1    12.2   76.3  12.9 0.964  11.0    14.8
#> 2     6.35  83.5  11.4 0.700  10.1    12.8
#> 3    28.0   60.2  25.4 0.874  23.7    27.1
#> 4    15.7   74.7  13.2 1.06   11.1    15.3
#> 5    16.2   85.9  10.9 0.696   9.65   12.3
#> 6    22.3   40.1  14.5 0.639  13.2    15.6
head(residuals(demo_fit))  # Residuals for each row of demo_fit$data
#> # A tibble: 6 × 6
#>   response  time  .mean   .sd .lower .upper
#>      <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3 -0.655 0.964 -2.58    1.23
#> 2     6.35  83.5 -5.06  0.700 -6.48   -3.74
#> 3    28.0   60.2  2.64  0.874  0.971   4.39
#> 4    15.7   74.7  2.46  1.06   0.348   4.57
#> 5    16.2   85.9  5.29  0.696  3.88    6.56
#> 6    22.3   40.1  7.85  0.639  6.70    9.08
log_lik(demo_fit)[1:3, 1:3]  # Log-likelihood at each demo_fit$data
#>              1         2         3
#> [1,] -2.503263 -2.748959 -2.493721
#> [2,] -2.358787 -3.159376 -2.418801
#> [3,] -2.308629 -3.061485 -2.470527

# All of the above take a range of arguments. E.g.,:
# \donttest{
head(predict(demo_fit))  # Pointwise posterior predictive
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  12.9  4.00   5.02   20.8
#> 2     6.35  83.5  11.4  3.95   3.65   19.2
#> 3    28.0   60.2  25.4  3.98  17.6    33.2
#> 4    15.7   74.7  13.2  4.03   5.31   21.1
#> 5    16.2   85.9  10.9  3.95   3.16   18.7
#> 6    22.3   40.1  14.5  3.94   6.74   22.2
head(predict(demo_fit, width = 0.8))  # 80% posterior predictive interval
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  12.9  4.00   7.78   18.0
#> 2     6.35  83.5  11.4  3.95   6.36   16.5
#> 3    28.0   60.2  25.4  3.98  20.3    30.5
#> 4    15.7   74.7  13.2  4.03   8.08   18.4
#> 5    16.2   85.9  10.9  3.95   5.87   16.0
#> 6    22.3   40.1  14.5  3.94   9.44   19.5
head(predict(demo_fit, prior = TRUE))  # Prior predictive
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  13.5  29.3  -42.7   76.7
#> 2     6.35  83.5  13.4  30.8  -46.6   81.2
#> 3    28.0   60.2  13.5  26.4  -37.9   68.0
#> 4    15.7   74.7  13.7  29.2  -41.0   76.2
#> 5    16.2   85.9  13.0  31.2  -48.3   80.8
#> 6    22.3   40.1  13.7  22.5  -31.2   55.4
head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#> # A tibble: 6 × 14
#>   response  time  .row .chain .iteration .draw  cp_1  cp_2 Intercept_1 time_2 Intercept_3  time_3 sigma_1 .epred
#>      <dbl> <dbl> <int>  <int>      <int> <int> <dbl> <dbl>       <dbl>  <dbl>       <dbl>   <dbl>   <dbl>  <dbl>
#> 1     12.2  76.3     1      1          1     1  32.5  69.6       10.2   0.582        10.8 -0.0577    4.49   10.4
#> 2     12.2  76.3     1      1          2     2  29.4  72.3        8.48  0.582        13.9 -0.197     4.14   13.1
#> 3     12.2  76.3     1      1          3     3  29.0  69.6        9.49  0.521        14.0 -0.199     3.99   12.7
#> 4     12.2  76.3     1      1          4     4  32.8  69.7       11.0   0.521        14.8 -0.277     4.05   12.9
#> 5     12.2  76.3     1      1          5     5  33.9  69.5       10.2   0.609        15.2 -0.276     4.16   13.3
#> 6     12.2  76.3     1      1          6     6  30.4  70.2        9.94  0.539        15.0 -0.226     4.04   13.6
head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#> # A tibble: 6 × 6
#>   response  time .mean   .sd .lower .upper
#>      <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1    12.2   76.3  3.88 0.277   3.39   4.51
#> 2     6.35  83.5  3.88 0.277   3.39   4.51
#> 3    28.0   60.2  3.88 0.277   3.39   4.51
#> 4    15.7   74.7  3.88 0.277   3.39   4.51
#> 5    16.2   85.9  3.88 0.277   3.39   4.51
#> 6    22.3   40.1  3.88 0.277   3.39   4.51

# Custom quantiles from the draws, here for fitted().
# Note that for predict(), this adds more Monte Carlo error than `width`.
fitted(demo_fit, summary = FALSE) |>
  dplyr::group_by(.row, time) |>
  dplyr::reframe(quantile = c(0.1, 0.5, 0.9), .epred = quantile(.epred, quantile)) |>
  head()
#> # A tibble: 6 × 4
#>    .row  time quantile .epred
#>   <int> <dbl>    <dbl>  <dbl>
#> 1     1  76.3      0.1   11.7
#> 2     1  76.3      0.5   12.9
#> 3     1  76.3      0.9   14.2
#> 4     2  83.5      0.1   10.5
#> 5     2  83.5      0.5   11.4
#> 6     2  83.5      0.9   12.3

# Evaluate at novel data
novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
head(predict(demo_fit, newdata = novel_data))
#> # A tibble: 3 × 5
#>    time  .mean   .sd .lower .upper
#>   <dbl>  <dbl> <dbl>  <dbl>  <dbl>
#> 1    -5   9.92  3.94   2.18  17.7 
#> 2    20   9.92  3.94   2.18  17.7 
#> 3   300 -33.3  17.4  -68.1    1.75

# Work with missing responses
missing_fit = mcp_example("missing", plot = FALSE)
#> NA values detected in 'y'. JAGS will treat them as latent responses and impute them during sampling.
fitted(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Expected responses for missing y
#> # A tibble: 6 × 7
#>       y     x state .mean   .sd .lower .upper
#>   <dbl> <int> <fct> <dbl> <dbl>  <dbl>  <dbl>
#> 1    NA     8 B      35.4 1.10    33.2   37.6
#> 2    NA    19 A      15.5 0.824   13.9   17.2
#> 3    NA    27 A      17.1 0.757   15.6   18.6
#> 4    NA    28 B      39.4 0.774   37.9   40.9
#> 5    NA    29 A      17.5 0.758   16.1   19.0
#> 6    NA    30 B      39.8 0.775   38.3   41.3
fitted(missing_fit, summary = FALSE) |> dplyr::filter(is.na(y)) |> head()  # Same, but draws
#> # A tibble: 6 × 14
#>       y     x state  .row .chain .iteration .draw  cp_1 Intercept_1   x_1 stateB_1    x_2 sigma_1 .epred
#>   <dbl> <int> <fct> <int>  <int>      <int> <int> <dbl>       <dbl> <dbl>    <dbl>  <dbl>   <dbl>  <dbl>
#> 1    NA     8 B         8      1          1     1  57.9        10.3 0.223     22.2 -0.419    4.33   34.4
#> 2    NA     8 B         8      1          2     2  58.5        11.8 0.187     22.3 -0.438    4.11   35.6
#> 3    NA     8 B         8      1          3     3  59.1        11.5 0.194     23.5 -0.473    4.53   36.5
#> 4    NA     8 B         8      1          4     4  60.9        11.7 0.168     22.9 -0.414    3.81   35.9
#> 5    NA     8 B         8      1          5     5  63.8        13.7 0.120     23.4 -0.463    4.53   38.0
#> 6    NA     8 B         8      1          6     6  63.5        12.1 0.160     21.7 -0.462    4.56   35.1
predict(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Posterior predictive for missing y
#> # A tibble: 6 × 7
#>       y     x state .mean   .sd .lower .upper
#>   <dbl> <int> <fct> <dbl> <dbl>  <dbl>  <dbl>
#> 1    NA     8 B      35.4  4.42  26.7    44.1
#> 2    NA    19 A      15.5  4.36   6.98   24.1
#> 3    NA    27 A      17.1  4.34   8.60   25.7
#> 4    NA    28 B      39.4  4.35  30.8    47.9
#> 5    NA    29 A      17.5  4.34   8.99   26.1
#> 6    NA    30 B      39.8  4.35  31.2    48.3
# }
```
