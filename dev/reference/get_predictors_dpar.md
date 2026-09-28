# Get predictors for one distributional parameter

This function extracts a `par_x`-less design matrix. x-terms are
measured from the change point, so `par_x` is multiplied in the formula
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

  Logical scalar. Whether to stop on rank deficiency.

- model:

  A list of formulas, one for each segment, in the format
  `model = list(response ~ predictors, cp ~ predictors)`. See
  [mcp-formula](https://lindeloev.github.io/mcp/dev/reference/mcp-formula.md)
  for how segments connect and the full formula syntax.

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
