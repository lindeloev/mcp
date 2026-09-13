# Get predictors for one distributional parameter

This function extracts a `par_x`-less design matrix. `par_x` will be
relative to the segment onset, so it will be multiplied in the formula
(`jags_code` and `fit$simulate()`).

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
  [`fitted.mcpfit`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  or
  [`predict.mcpfit`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
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

  Logical scalar. Whether to stop on rank deficiency.

- model:

  A list of formulas - one for each segment. The general format is
  `response ~ cp ~ predictors` (e.g., `y ~ 1 ~ 1 + x`), except the first
  segment has no change point and uses `response ~ predictors`. The
  response and change-point parts can be omitted (`cp ~ predictor`
  assumes the same response; `~ predictor` assumes an intercept-only
  change point). Non-\$x\$ terms persist into later segments until
  replaced or removed by an intercept reset (see details). See examples
  on the [mcp website](https://lindeloev.github.io/mcp/).

  **1. Response (segment 1 only):**

  - `y ~ ...`: Standard continuous or count response (Gaussian, Poisson,
    Bernoulli).

  - `successes | trials(total) ~ ...`: Binomial response
    (`family = binomial()`).

  - `y | weights(w) ~ ...`: Observation log-likelihood weights
    (multiplies each observation's log-likelihood contribution by
    `w > 0`; affects posterior inference and
    [`log_lik()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
    but not predictions).

  - `y | trials(total) + weights(w) ~ ...`: Combine response auxiliaries
    using `+`.

  **2. Change-point modeling (`cp`, segments 2+):**

  - `~ 1 ~ ...` (or omitted, e.g., `~ x`): Population-level change point
    (default).

  - `1 + (1 | id) ~ ...`: Group-level change-point deviations around the
    population change point. [Read
    more](https://lindeloev.github.io/mcp/articles/group_effects.html).

  **3. Regression formula (all segments):** [Read
  more](https://lindeloev.github.io/mcp/articles/formulas.html)

  - `~ 1 + x`: Disjoined slope with a new segment intercept.

  - `~ 0 + x`: Joined slope (no intercept; continuous from previous
    segment).

  - `~ 1`: Plateau (intercept only, no slope).

  - `~ x:group + I(x^2) + exp(z)`: Extended terms, interactions, and
    R-side bases ([`scale()`](https://rdrr.io/r/base/scale.html),
    [`poly()`](https://rdrr.io/r/stats/poly.html),
    [`splines::ns()`](https://rdrr.io/r/splines/ns.html)). Bases are
    evaluated before sampling and reused for `newdata`.

  - `~ 1 + (1 | id)`: Group-level intercepts (or `(1 + x || id)` for
    independent slopes and intercepts). [Read
    more](https://lindeloev.github.io/mcp/articles/group_effects.html).

  - `~ sigma(1 + x)`: Distributional parameters on the link scale (e.g.,
    log residual SD). [Read
    more](https://lindeloev.github.io/mcp/articles/dpar.html).

  - `~ ar(1) + ma(1)`: Autoregressive and moving-average time-series
    residuals on the link scale (accepts regression formulas,
    `series = id`, and `boundary`; use `ar(0)` or `ma(0)` to turn off in
    later segments). [Read
    more](https://lindeloev.github.io/mcp/articles/arma.html).

- family:

  A supported family:
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`bernoulli()`](https://lindeloev.github.io/mcp/reference/bernoulli.md),
  [`poisson()`](https://rdrr.io/r/stats/family.html), or
  [`negbinomial()`](https://lindeloev.github.io/mcp/reference/negbinomial.md),
  with a supported link function; e.g., `gaussian(link = "log")`.

## Value

A tibble with one row per model parameter and the columns

- `dpar`: character.

- `segment`: the segment number (positive integer).

- `matrix_name`: original column name from the model matrix. Used to
  diagnose collisions after parameter names are converted for JAGS.

- `display_name`: user-facing parameter name used in summary functions.

- `code_name`: parameter name used in JAGS and internally in mcp.

- `term_key`: identifier for the formula term which generated the
  coefficient. Multi-column terms share one key.

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
