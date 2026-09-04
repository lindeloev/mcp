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

  Optional data frame at which to evaluate the model. For GARMA
  `posterior_predict()`, only predictors and required response
  auxiliaries are needed: each response series is generated recursively
  without conditioning on an observed response column.

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

- ...:

  Must be empty. Reserved for future use.

- transform:

  For `posterior_linpred()`, return the inverse-link transformed
  expected response instead of the linear predictor.

## Value

A numeric `N_draws` by `nrow(newdata)` matrix.

## Details

For GARMA models, `posterior_predict()` generates each replicated
response series recursively. It does not condition later predictions on
the observed response history, unlike
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
and posterior
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).
These methods require posterior draws. For prior prediction, use
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
with `prior = TRUE`, which also generates fresh series.

For binomial models, `posterior_epred()` and `posterior_predict()` (and
corresponding `{tidybayes}` workflows such as
[`add_epred_draws()`](https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html))
follow `{brms}` and `{rstantools}` conventions by returning values on
the outcome count scale (`rate = FALSE`), i.e., expected counts \\E\[Y\]
= n\mu\\ and simulated counts in \\\\0, \dots, n\\\\. In contrast,
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
and
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
default to proportions (`rate = TRUE`). To obtain the success
probability parameter \\\mu\\ on the \\\[0, 1\]\\ scale regardless of
trial counts, pass `dpar = "mu"`.

## See also

[`fitted.mcpfit()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
[`predict.mcpfit()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
