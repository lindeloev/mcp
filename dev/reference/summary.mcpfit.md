# Summarise mcpfit objects

Summarise parameter estimates and model diagnostics.

## Usage

``` r
# S3 method for class 'mcpfit'
summary(
  object,
  width = 0.95,
  digits = 2,
  prior = FALSE,
  verbose = FALSE,
  diagnostics = NULL,
  ...
)

# S3 method for class 'mcpsummary'
print(x, digits = attr(x, "mcpsummary")$digits, ...)

# S3 method for class 'mcpfit'
fixef(object, width = 0.95, prior = FALSE, verbose = FALSE, dpar = "mu", ...)

# S3 method for class 'mcpfit'
ranef(object, width = 0.95, prior = FALSE, verbose = FALSE, ...)

# S3 method for class 'mcpfit'
print(x, ...)
```

## Arguments

- object:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object.

- width:

  Float. The width of the central posterior interval (between 0 and 1).

- digits:

  Non-negative integer. Number of significant digits used when printing
  the summary. Defaults to 2. The returned tibble retains the unrounded
  values.

- prior:

  Logical. Summarise prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)?

- verbose:

  Logical. Include the `segment` and `dpar` columns. Defaults to `FALSE`
  for a compact summary.

- diagnostics:

  Named list of diagnostic warning thresholds. Available elements are
  `rhat = 1.01`, `ess_bulk = 400`, `ess_tail = 400`, `ar = 0.10`, and
  `ma = 0.10`. An empty list uses these defaults; a partial list
  overrides only the supplied values. Set an element to `NULL` to
  disable that diagnostic, or use `FALSE` to disable all configurable
  diagnostic warnings. In `summary.mcpfit()`, `NULL` inherits the
  settings used to fit the model, while a list or `FALSE` overrides the
  diagnostic footer.

- ...:

  Must be empty. Reserved for future use.

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object, or an `mcpsummary` returned by `summary()`.

- dpar:

  Distributional parameter(s) whose regression coefficients to return,
  e.g., `"mu"`, `"sigma"`, or `"ar1"` (see
  [mcp-formula](https://lindeloev.github.io/mcp/dev/reference/mcp-formula.md)).
  For modeled distributional parameters such as
  [`sigma()`](https://rdrr.io/r/stats/sigma.html), these coefficients
  are on the link scale.

## Value

An `mcpsummary` tibble with parameter estimates and MCMC diagnostics. It
prints with model information, separate tables for change points and
population-level parameters, and convergence warnings. Rows are ordered
by change point first, then `mu`, then the other distributional
parameters, then AR/MA lags - each ascending by segment. Note: The
change point distributions are often not unimodal and symmetric so the
intervals can be deceiving. Plot them using `plot_pars(fit)`.

Columns:

- `segment` (only with `verbose = TRUE`) is the segment the parameter
  belongs to. For change points, it is the segment that the change point
  starts, e.g., 2 for `cp_1`.

- `dpar` (only with `verbose = TRUE`) is the distributional parameter
  (e.g., `"mu"`, `"sigma"`, or `"ar1"`; see
  [mcp-formula](https://lindeloev.github.io/mcp/dev/reference/mcp-formula.md))
  as accepted by the `dpar` argument of `fixef.mcpfit()`,
  [`fitted.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  etc. It is `NA` for change points. See also
  [`mcp_pars()`](https://lindeloev.github.io/mcp/dev/reference/mcp_pars.md).

- `mean` is the posterior mean

- `sd` is the posterior standard deviation.

- `lower` and `upper` are the lower and upper ends of the central
  posterior interval given in `width`.

- `rhat` is the rank-normalized split-Rhat convergence diagnostic.

- `ess_bulk` and `ess_tail` are the bulk and tail effective sample
  sizes. Low effective sample sizes are also obvious as poor mixing in
  trace plots (see `plot_pars(fit)`). Read how to deal with such
  problems in [Diagnosing and fixing
  problems](https://lindeloev.github.io/mcp/articles/tips.html#convergence).

Group-level change-point deviations (`cp_i_id`) follow a hierarchical
normal distribution around the population-level change point before
truncation. The group-specific change points are truncated to the
observed range of x and ordered, so the population-level change point is
generally close to their mean, but not exactly. Predictor group-level
deviations (such as `Intercept_1_id`) use standard hierarchical
zero-mean priors, without these bounds. [Read
more](https://lindeloev.github.io/mcp/articles/group_effects.html).

For simulated data, the summary contains two additional columns so that
it is easy to inspect whether the model can recover the parameters. Run
simulation and summary multiple times to get a sense of the robustness.

- `sim` is the value used to generate the data.

- `match` is `"OK"` if `sim` is contained in the central posterior
  interval (`lower` to `upper`). For change points starting a disjoined
  segment, `sim` only needs to be between the same two adjacent observed
  `x` as the interval, since the likelihood is identical anywhere
  between them.

## Functions

- `print(mcpsummary)`: Print an `mcpsummary` returned by `summary()`:
  model information, change point and population-level parameters, and
  diagnostics.

- `fixef(mcpfit)`: Population-level coefficients (fixed effects) of
  `mcpfit`.

- `ranef(mcpfit)`: Group-level deviations (random effects) of `mcpfit`.
  Change-point deviations are relative to their population-level change
  point; `cp_i_<group>_sd` (e.g., `cp_1_id_sd`) is the SD of their
  normal distribution before truncation.

- `print(mcpfit)`: Print the posterior summary of an
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object, or the prior summary with `prior = TRUE`. Without the
  requested draws, only the model is shown.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# Typical usage
summary(demo_fit)
#> Family: gaussian
#> Links: mu = identity; sigma = identity
#> Draws: 1000 (2 chains with iter = 500)
#> Segments:
#>   1: response ~ 1
#>   2: ~ 0 + time
#>   3: ~ 1 + time
#> 
#> Change point parameters:
#>    parameter  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  cp_1        31.64 1.691 28.23 34.883 1.00     1012      945 30.0    OK
#>  cp_2        71.10 0.991 69.48 72.770 1.01      869      676 70.0    OK
#> 
#> Population-level parameters:
#>    parameter  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.97 0.665  8.68 11.245 1.00     1032      747 10.0    OK
#>  time_2       0.54 0.048  0.45  0.643 1.00     1087      970  0.5    OK
#>  Intercept_3 13.97 1.386 11.11 16.647 1.00      993      823 15.0    OK
#>  time_3      -0.21 0.084 -0.37 -0.048 1.00      900      740 -0.3    OK
#>  sigma_1      3.86 0.293  3.34  4.540 1.00      931      921  3.5    OK
summary(demo_fit, width = 0.8, digits = 4)  # Set interval width
#> Family: gaussian
#> Links: mu = identity; sigma = identity
#> Draws: 1000 (2 chains with iter = 500)
#> Segments:
#>   1: response ~ 1
#>   2: ~ 0 + time
#>   3: ~ 1 + time
#> 
#> Change point parameters:
#>    parameter    mean      sd   lower   upper   rhat ess_bulk ess_tail  sim match
#>  cp_1        31.6369 1.69122 29.4401 33.7774 1.0028     1012      945 30.0    OK
#>  cp_2        71.0983 0.99102 69.7097 72.4884 1.0065      869      676 70.0    OK
#> 
#> Population-level parameters:
#>    parameter    mean      sd   lower   upper   rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.9675 0.66463  9.1113 10.8454 1.0003     1032      747 10.0    OK
#>  time_2       0.5436 0.04754  0.4826  0.6026 1.0032     1087      970  0.5    OK
#>  Intercept_3 13.9661 1.38591 12.2422 15.7274 0.9994      993      823 15.0    OK
#>  time_3      -0.2062 0.08367 -0.3078 -0.1008 1.0003      900      740 -0.3    OK
#>  sigma_1      3.8615 0.29314  3.5054  4.2163 1.0027      931      921  3.5      

# Get the results as a tibble
results = summary(demo_fit)
results$mean
#> [1] 31.6369306 71.0983320  9.9674657  0.5435759 13.9660544 -0.2062269  3.8615146

# Group-level deviations (random effects)
# ranef(my_fit)

# Summarise prior
summary(demo_fit, prior = TRUE)
#> Family: gaussian
#> Links: mu = identity; sigma = identity
#> Draws: 1000 (2 chains with iter = 500)
#> Segments:
#>   1: response ~ 1
#>   2: ~ 0 + time
#>   3: ~ 1 + time
#> 
#> Change point parameters:
#>    parameter    mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        34.6020 23.52   1.99  85.8 1.00     1004      882 30.0    OK
#>  cp_2        65.9423 23.39  16.38  98.1 1.00      975      916 70.0    OK
#> 
#> Population-level parameters:
#>    parameter    mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 13.7421 17.74 -22.83  48.1 1.00      786      822 10.0    OK
#>  time_2       0.0068  0.62  -1.23   1.3 1.00      859      956  0.5    OK
#>  Intercept_3 13.4314 17.44 -19.85  47.6 1.00      920      834 15.0    OK
#>  time_3      -0.0109  0.63  -1.25   1.2 1.00     1031      980 -0.3    OK
#>  sigma_1      6.5183  9.00   0.13  25.8 1.00      912     1028  3.5    OK
```
