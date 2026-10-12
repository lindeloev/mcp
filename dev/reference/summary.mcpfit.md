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

  Distributional parameter(s) whose regression coefficients to return.
  For modeled distributional parameters such as
  [`sigma()`](https://rdrr.io/r/stats/sigma.html), these coefficients
  are on the link scale.

## Value

An `mcpsummary` tibble with parameter estimates and MCMC diagnostics. It
prints with model information, separate tables for change points and
population-level parameters, and convergence warnings. Rows are ordered
by change point first, then `mu`, then the other distributional
parameters, then `ar`/`ma` components - each ascending by segment. Note:
The change point distributions are often not unimodal and symmetric so
the intervals can be deceiving. Plot them using `plot_pars(fit)`.

With `verbose = TRUE`:

- `segment` is the segment the parameter belongs to.

- `dpar` is the distributional parameter (`"cp"`, `"mu"`, `"sigma"`,
  `"ar"`, `"ma"`, etc.) the parameter belongs to. For AR/MA terms, the
  lag order is encoded in `parameter`, e.g. `ar2_1`.

- `mean` is the posterior mean

- `sd` is the posterior standard deviation.

- `lower` and `upper` are the bounds of the central posterior interval
  given in `width`.

- `rhat` is the rank-normalized split-Rhat convergence diagnostic.

- `ess_bulk` and `ess_tail` are the bulk and tail effective sample
  sizes. Low effective sample sizes are also obvious as poor mixing in
  trace plots (see `plot_pars(fit)`). Read how to deal with such
  problems [here](https://lindeloev.github.io/mcp/articles/tips.html)

Group-level change-point deviations (`cp_i_id`) follow a hierarchical
normal distribution around the population change point before
truncation. Realized locations are truncated to remain in range and
ordered, so the population change point is generally close to the mean
of the realized group change points, but not exactly. Predictor
group-level effects (such as `Intercept_1_id`) use standard hierarchical
zero-mean priors, without change-point constraints.

For simulated data, the summary contains two additional columns so that
it is easy to inspect whether the model can recover the parameters. Run
simulation and summary multiple times to get a sense of the robustness.

- `sim` is the value used to generate the data.

- `match` is `"OK"` if `sim` is contained in the central posterior
  interval (`lower` to `upper`).

## Functions

- `print(mcpsummary)`: Print an `mcpsummary` returned by `summary()`:
  model information, change point and population-level parameters, and
  diagnostics.

- `fixef(mcpfit)`: Population-level fixed effects (regression
  coefficients) of `mcpfit`.

- `ranef(mcpfit)`: Group-level deviations (random effects) of `mcpfit`.
  Change-point deviations are relative to their population change point;
  `cp_i_<group>_sd` (e.g., `cp_1_id_sd`) is the SD of their normal
  distribution before truncation.

- `print(mcpfit)`: Print the posterior summary of an
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object.

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
#>  cp_1        31.59 1.744 28.24 34.914 1.00      950      907 30.0    OK
#>  cp_2        71.14 1.015 69.48 72.788 1.00     1055     1029 70.0    OK
#> 
#> Population-level parameters:
#>    parameter  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.95 0.665  8.49 11.223 1.00      894      973 10.0    OK
#>  time_2       0.54 0.047  0.46  0.634 1.00      979      776  0.5    OK
#>  Intercept_3 13.94 1.367 11.30 16.700 1.00      951      824 15.0    OK
#>  time_3      -0.21 0.082 -0.36 -0.044 1.00      869      837 -0.3    OK
#>  sigma_1      3.87 0.280  3.37  4.442 1.00      890      945  3.5    OK
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
#>    parameter    mean      sd   lower    upper   rhat ess_bulk ess_tail  sim
#>  cp_1        31.5935 1.74401 29.3232 33.82035 1.0005      950      907 30.0
#>  cp_2        71.1442 1.01458 69.7466 72.53458 0.9991     1055     1029 70.0
#>  match
#>     OK
#>     OK
#> 
#> Population-level parameters:
#>    parameter    mean      sd   lower    upper   rhat ess_bulk ess_tail  sim
#>  Intercept_1  9.9465 0.66456  9.0751 10.77530 1.0033      894      973 10.0
#>  time_2       0.5437 0.04660  0.4836  0.60348 0.9993      979      776  0.5
#>  Intercept_3 13.9424 1.36672 12.1993 15.70817 0.9998      951      824 15.0
#>  time_3      -0.2066 0.08208 -0.3136 -0.09996 0.9996      869      837 -0.3
#>  sigma_1      3.8675 0.28029  3.5116  4.23835 0.9997      890      945  3.5
#>  match
#>     OK
#>     OK
#>     OK
#>     OK
#>       

# Get the results as a tibble
results = summary(demo_fit)
results$mean
#> [1] 31.5935403 71.1441814  9.9464536  0.5437021 13.9424489 -0.2065739  3.8674841

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
#>    parameter   mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        33.161 22.78   2.10  82.9 1.00      901     1074 30.0    OK
#>  cp_2        64.742 23.34  15.21  98.3 1.00      904      927 70.0    OK
#> 
#> Population-level parameters:
#>    parameter   mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 13.671 18.03 -21.28  50.1 1.00     1039      982 10.0    OK
#>  time_2      -0.020  0.61  -1.21   1.2 1.00      946     1017  0.5    OK
#>  Intercept_3 14.016 17.63 -20.72  47.1 1.00      903      911 15.0    OK
#>  time_3       0.012  0.62  -1.21   1.2 1.00      884      826 -0.3    OK
#>  sigma_1      5.906  6.61   0.15  22.5 1.00      901      784  3.5    OK
```
