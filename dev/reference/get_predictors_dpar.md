# Get predictors for one distributional parameter

This function extracts an `par_x`-less design matrix. `par_x` will be
relative to the segment onset, so it will be multiplied in in the
formula (`jags_code` and `fit$simulate()`).

## Usage

``` r
get_predictors_dpar(
  data,
  form_rhs,
  segment,
  dpar,
  par_x,
  order = NULL,
  check_rank = TRUE,
  design_id = NULL
)

get_predictor_tables(model, data, family, par_x, check_rank = TRUE)

get_predictors(model, data, family, par_x, check_rank = TRUE)
```

## Arguments

- data:

  Table-like data in long format (data.frame, tibble, data.table, etc.)
  with syntactic column names. Missing values in the response variable
  are imputed using the posterior predictive.
  [`fitted.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  or
  [`predict.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  details how to see the imputed values.

- form_rhs:

  The full predictor formula of a segment, including one or several
  distributional terms.

- segment:

  Integer. The segment number

- dpar:

  A distributional parameter or an `ar`/`ma` component.

- par_x:

  String (default: `NULL` which is auto-detect).

- order:

  Applies to `dpar %in% c("ar", "ma")`.

- check_rank:

  Boolean. Whether to stop on rank deficiency.

- model:

  A list of formulas - one for each segment. The many examples on the
  [mcp website](https://lindeloev.github.io/mcp/). But briefly:

  The first formula has the format `response ~ predictors` while the
  following formulas have the format `response ~ cp ~ predictors`. Here,
  `cp` names the change-point part of the formula rather than a literal
  variable. The response and change-point parts can be omitted
  (`cp ~ predictor` assumes the same response; `~ predictor` assumes an
  intercept-only change point). Terms normally carry into later segments
  until redefined (see details).

  The following terms can be modeled:

  - *Regular formulas:* e.g., `~ 1 + x`. [Read
    more](https://lindeloev.github.io/mcp/articles/formulas.html).

  - *Extended formulas:* e.g., `~ x:group + I(x^2) + exp(z)`. [Read
    more](https://lindeloev.github.io/mcp/articles/formulas.html).
    R-side bases such as [`scale()`](https://rdrr.io/r/base/scale.html),
    [`poly()`](https://rdrr.io/r/stats/poly.html), and
    [`splines::ns()`](https://rdrr.io/r/splines/ns.html) are evaluated
    before sampling, and their fitted scaling or basis is reused for
    `newdata`.

  - *Group-level effects (random effects):* e.g., `~ 1 + (1 | id)` for a
    group-level intercept, or `~ 1 + (factor || id)` for independent
    intercept and factor-contrast deviations. [Read
    more](https://lindeloev.github.io/mcp/articles/group_effects.html).

  - *Gaussian residual standard deviation:* e.g., `~sigma(1)` for a
    simple standard-deviation change or `~sigma(1 + x + group)` for more
    advanced structures. Explicit
    [`sigma()`](https://rdrr.io/r/stats/sigma.html) formulas model
    log-SD, while the implicit constant `sigma_1` in a model without
    [`sigma()`](https://rdrr.io/r/stats/sigma.html) remains on the
    response scale. [Read
    more](https://lindeloev.github.io/mcp/articles/dpar.html)

  - *Time-series residuals:* use `ar(p)` and `ma(q)` separately or
    together, e.g., `~ 1 + ar(1) + ma(1)`. Both accept an optional
    regression formula and observation `boundary`. GARMA terms support
    Gaussian, binomial, Poisson, and negative-binomial families with
    their default links. They define a finite conditional recurrence on
    the link scale. For `~ ar(p) + ma(q)`, the model is:

    \$\$\eta_t = b_t + \sum\_{j=1}^{p} \phi\_{j,t}
    \left\[g(y^\*\_{t-j}) - b\_{t-j}\right\] + \sum\_{k=1}^{q}
    \theta\_{k,t} \left\[g(y^\*\_{t-k}) - \eta\_{t-k}\right\]\$\$

    where \\b_t\\ is the linear predictor from the segment formulas,
    \\\phi\_{j,t}\\ is the lag-\\j\\ autoregressive (AR) coefficient at
    time \\t\\, \\\theta\_{k,t}\\ is the lag-\\k\\ moving-average (MA)
    coefficient at time \\t\\, \\g(\cdot)\\ is the link function,
    \\y^\*\_t\\ is the observation, and \\\eta_t\\ is the resulting full
    linear predictor including serial dependence. Note some
    implications:

    - For an N-order component, the last N values *before* the segment
      onset are input to the first \\\eta_t\\ in the segment.

    - AR coefficients are not jointly constrained to stationarity; nor
      MA coefficients to invertibility. [Read
      more](https://lindeloev.github.io/mcp/articles/arma.html)

  - *Weights:* `y | weights(w) ~ ...` specifies observation
    log-likelihood weights.

  - *Binomial:* use `successes | trials(total) ~ ...` with
    `family = binomial()`.

- family:

  A supported family:
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`bernoulli()`](https://lindeloev.github.io/mcp/dev/reference/bernoulli.md),
  [`poisson()`](https://rdrr.io/r/stats/family.html), or
  [`negbinomial()`](https://lindeloev.github.io/mcp/dev/reference/negbinomial.md),
  with a supported link function; e.g., `gaussian(link = "log")`.

## Value

A tibble with one row per model parameter and the columns

- `dpar`: character.

- `segment`: the segment number (positive integer).

- `matrix_name`: original column name from the model matrix. Used to
  diagnose collisions after parameter names are converted for JAGS.

- `display_name`: user-facing parameter name used in summary functions.

- `code_name`: parameter name used in JAGS and internally in mcp.

- `par_type`: One of "Intercept", "dummy", or "slope". Used for setting
  priors and for change point indicator func.

- `order`: positive integer or NA. Only relevant for `ar` and `ma`.

- `explicit`: whether the distributional parameter was supplied in the
  formula.

- `design_id`: key of the fitted component formula that produced the
  row.

- `design_col`: column occupied by the row in that component's model
  matrix.

- `matrix_data`: column of the design matrix less the `par_x` term.

## Functions

- `get_predictor_tables()`: Apply `get_predictors_segment` to all
  segments of a model.

- `get_predictors()`: Return only the population predictor table.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
