# Posterior prediction draws for `mcpfit` objects

Methods for the `{rstantools}` posterior-prediction generics. They
return a draws-by-observation matrix and enable `{tidybayes}` workflows
such as
[`add_epred_draws()`](https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html),
[`add_predicted_draws()`](https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html),
and
[`add_linpred_draws()`](https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html).
These methods and workflows require the suggested package
`{rstantools}`.

## Usage

``` r
posterior_epred.mcpfit(
  object,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  dpar = NULL,
  seed = NULL,
  rate = FALSE,
  ...
)

posterior_predict.mcpfit(
  object,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  seed = NULL,
  rate = FALSE,
  conditional = TRUE,
  ...
)

posterior_linpred.mcpfit(
  object,
  transform = FALSE,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  dpar = NULL,
  seed = NULL,
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
    [`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
    [`residuals()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
    [`log_lik()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
    and
    [`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
    condition on the response history by default, so `newdata` must
    include the response. For
    [`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
    and
    [`residuals()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
    missing response histories are supported only in the original fitted
    data, using retained posterior imputations as histories. Predictions
    are fresh response draws, including at missing rows. With
    `conditional = FALSE`,
    [`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
    and
    [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
    generate fresh response series recursively, so their `newdata` need
    only contain predictors and required response auxiliaries.
    [`log_lik()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
    is unavailable when a missing response enters a later observed
    history.

  - For models with `y | weights()`: Require the weights column except
    for
    [`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
    and
    [`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).

- draws, ndraws:

  Number of posterior draws to return. `draws` follows the
  `{rstantools}` convention; `ndraws` is the mcp spelling. Supply at
  most one.

- re.form, re_formula:

  Group-level effects to include. `NULL` includes all effects and `NA`
  excludes them.

- dpar:

  Distributional parameter for `posterior_epred()` and
  `posterior_linpred()`; `NULL` uses the expected response.

- seed:

  Optional integer seed for draw selection and posterior prediction.

- rate:

  Logical scalar. For binomial models, return counts (`rate = FALSE`,
  the default for
  [`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  and
  [`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md))
  or the observed or expected success proportion (`rate = TRUE`).
  Predictions and count-scale fitted values require a trials column in
  `newdata`. Distributional parameters such as `dpar = "mu"` evaluate
  the parameter itself (e.g., success probability) and are unaffected by
  `rate`.

- ...:

  Must be empty. Reserved for future use.

- conditional:

  Logical. For AR/MA models, condition on observed response histories
  (`TRUE`, default) or generate fresh histories recursively (`FALSE`).
  Applies equally to prior and posterior draws. Predictive checks with
  [`pp_check()`](https://lindeloev.github.io/mcp/reference/pp_check.md)
  generate fresh histories.

- transform:

  For `posterior_linpred()`, return the inverse-link transformed
  expected response instead of the linear predictor.

## Value

A numeric `N_draws` by `nrow(newdata)` matrix.

## Details

For GARMA models, `posterior_predict()` conditions on the observed
response history, just like
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).
Use `conditional = FALSE` in either method to generate fresh response
histories recursively. Missing responses in the original data are filled
with retained imputations only to supply histories. Conditional
`posterior_predict()` draws from the response distributions whose means
`posterior_epred()` returns, including at missing rows. These methods
require posterior draws. For prior prediction, use
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
with `prior = TRUE`; `conditional` selects the same behavior.

For binomial models, `posterior_epred()` and `posterior_predict()` (and
corresponding `{tidybayes}` workflows such as
[`add_epred_draws()`](https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html))
follow `{brms}` and `{rstantools}` conventions by returning values on
the outcome count scale (`rate = FALSE`), i.e., expected counts \\E\[Y\]
= n\mu\\ and simulated counts in \\\\0, \dots, n\\\\, matching
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
and
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).
Use `rate = TRUE` for proportions. To obtain the success probability
parameter \\\mu\\ on the \\\[0, 1\]\\ scale regardless of trial counts,
pass `dpar = "mu"`.

## See also

[`fitted.mcpfit()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
[`predict.mcpfit()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
