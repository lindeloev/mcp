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
  the summary. Defaults to 2. The invisibly returned data frame retains
  the unrounded values.

- prior:

  Logical. Summarise prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)?

- verbose:

  Logical. Include the `segment` and `dpar` columns. Defaults to `FALSE`
  for a compact, v0.3.4-compatible summary.

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

- dpar:

  Distributional parameter(s) whose regression coefficients to return.
  For modeled distributional parameters such as
  [`sigma()`](https://rdrr.io/r/stats/sigma.html), these coefficients
  are on the link scale.

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object.

## Value

A data frame with parameter estimates and MCMC diagnostics. Rows are
ordered by change point first, then `mu`, then the other distributional
parameters, then `ar`/`ma` components - each ascending by segment. OBS:
The change point distributions are often not unimodal and symmetric so
the intervals can be deceiving. Plot them using `plot_pars(fit)`.

With `verbose = TRUE`:

- `segment` is the segment the parameter belongs to.

- `dpar` is the distributional parameter (`"cp"`, `"mu"`, `"sigma"`,
  `"ar"`, `"ma"`, etc.) the parameter belongs to. For AR/MA terms, the
  lag order is encoded in `variable`, e.g. `ar2_1`.

- `mean` is the posterior mean

- `sd` is the posterior standard deviation.

- `lower` and `upper` are the bounds of the central posterior interval
  given in `width`.

- `rhat` is the rank-normalized split-Rhat convergence diagnostic.

- `ess_bulk` and `ess_tail` are the bulk and tail effective sample
  sizes. Low effective sample sizes are also obvious as poor mixing in
  trace plots (see `plot_pars(fit)`). Read how to deal with such
  problems [here](https://lindeloev.github.io/mcp/articles/tips.html)

Group-level change-point deviations (`cp_i_id`) follow a standard
hierarchical normal distribution around the population change point.
Their realized locations are truncated to remain in range and ordered.
Predictor group-level effects (such as `Intercept_1_id`) also use
standard hierarchical zero-mean priors, without change-point
constraints.

For simulated data, the summary contains two additional columns so that
it is easy to inspect whether the model can recover the parameters. Run
simulation and summary multiple times to get a sense of the robustness.

- `sim` is the value used to generate the data.

- `match` is `"OK"` if `sim` is contained in the central posterior
  interval (`lower` to `upper`).

## Functions

- `fixef(mcpfit)`: Population-level fixed effects (regression
  coefficients) of `mcpfit`.

- `ranef(mcpfit)`: Group-level deviations (random effects) of `mcpfit`.
  Change-point deviations are relative to their population change point;
  `cp_i_sd` is the scale of their latent normal distribution.

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
#> Iterations: 500 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>     variable  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  cp_1        31.46 1.729 28.06 34.998 1.00      726      951 30.0    OK
#>  cp_2        71.12 1.016 69.44 72.765 1.00      867      958 70.0    OK
#> 
#> Population-level parameters:
#>     variable  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 10.05 0.648  8.78 11.324 1.00      992      921 10.0    OK
#>  time_2       0.53 0.046  0.44  0.619 1.00      898      858  0.5    OK
#>  Intercept_3 17.46 1.175 15.38 19.965 1.00      799      941 20.0      
#>  time_3      -0.10 0.073 -0.26  0.019 1.00      836      860 -0.3      
#>  sigma_1      3.89 0.276  3.38  4.466 1.00     1043     1121  3.5    OK
summary(demo_fit, width = 0.8, digits = 4)  # Set interval width
#> Family: gaussian
#> Links: mu = identity; sigma = identity
#> Iterations: 500 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>     variable    mean      sd   lower   upper   rhat ess_bulk ess_tail  sim
#>  cp_1        31.4599 1.72856 29.4131 33.5661 1.0024      726      951 30.0
#>  cp_2        71.1172 1.01571 69.6973 72.5154 0.9987      867      958 70.0
#>  match
#>     OK
#>     OK
#> 
#> Population-level parameters:
#>     variable    mean      sd   lower   upper   rhat ess_bulk ess_tail  sim
#>  Intercept_1 10.0470 0.64793  9.1995 10.8528 0.9986      992      921 10.0
#>  time_2       0.5303 0.04610  0.4708  0.5878 1.0006      898      858  0.5
#>  Intercept_3 17.4634 1.17476 15.9842 19.0219 1.0001      799      941 20.0
#>  time_3      -0.1021 0.07325 -0.2060 -0.0142 0.9998      836      860 -0.3
#>  sigma_1      3.8934 0.27572  3.5711  4.2535 0.9985     1043     1121  3.5
#>  match
#>     OK
#>     OK
#>       
#>       
#>       

# Get the results as a data frame
results = summary(demo_fit)
#> Family: gaussian
#> Links: mu = identity; sigma = identity
#> Iterations: 500 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>     variable  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  cp_1        31.46 1.729 28.06 34.998 1.00      726      951 30.0    OK
#>  cp_2        71.12 1.016 69.44 72.765 1.00      867      958 70.0    OK
#> 
#> Population-level parameters:
#>     variable  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 10.05 0.648  8.78 11.324 1.00      992      921 10.0    OK
#>  time_2       0.53 0.046  0.44  0.619 1.00      898      858  0.5    OK
#>  Intercept_3 17.46 1.175 15.38 19.965 1.00      799      941 20.0      
#>  time_3      -0.10 0.073 -0.26  0.019 1.00      836      860 -0.3      
#>  sigma_1      3.89 0.276  3.38  4.466 1.00     1043     1121  3.5    OK

# Group-level deviations (random effects)
# ranef(my_fit)

# Summarise prior
summary(demo_fit, prior = TRUE)
#> Family: gaussian
#> Links: mu = identity; sigma = identity
#> Iterations: 500 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>     variable     mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        33.82751 23.81  2.01 85.25 1.00      992      948 30.0    OK
#>  cp_2        65.97862 23.39 17.59 98.37 1.00      809     1067 70.0    OK
#> 
#> Population-level parameters:
#>     variable     mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 14.16572 10.68 -6.80 32.94 1.00      927      667 10.0    OK
#>  time_2      -0.00596  0.12 -0.22  0.16 1.00     1042      945  0.5      
#>  Intercept_3 13.83889 10.09 -4.97 31.09 1.00      891      879 20.0    OK
#>  time_3       0.00022  0.10 -0.18  0.20 1.00      945      854 -0.3      
#>  sigma_1      6.51879  8.85  0.28 26.16 1.00     1045     1026  3.5    OK
```
