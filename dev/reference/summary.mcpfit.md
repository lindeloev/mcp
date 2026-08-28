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

  Currently ignored.

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
  lag order is encoded in `name`, e.g. `ar2_1`.

- `mean` is the posterior mean

- `sd` is the posterior standard deviation.

- `lower` and `upper` are the bounds of the central posterior interval
  given in `width`.

- `rhat` is the rank-normalized split-Rhat convergence diagnostic.

- `ess_bulk` and `ess_tail` are the bulk and tail effective sample
  sizes. Low effective sample sizes are also obvious as poor mixing in
  trace plots (see `plot_pars(fit)`). Read how to deal with such
  problems [here](https://lindeloev.github.io/mcp/articles/tips.html)

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

- `print(mcpfit)`: Print the posterior summary of an
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# Typical usage
summary(demo_fit)
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  cp_1        30.45 3.606 23.28 37.885 1.01      135      249 30.0    OK
#>  cp_2        69.76 0.288 69.29 70.250 1.00     1424     1529 70.0    OK
#> 
#> Population-level parameters:
#>         name  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 10.30 0.717  8.87 11.670 1.01      389      560 10.0    OK
#>  time_2       0.53 0.060  0.43  0.666 1.01      148      278  0.5    OK
#>  Intercept_3  0.62 1.558 -2.38  3.723 1.01      193      347  0.0    OK
#>  time_3      -0.22 0.092 -0.40 -0.049 1.01      182      388 -0.2    OK
#>  sigma_1      4.01 0.311  3.45  4.657 1.00      795      854  4.0    OK
summary(demo_fit, width = 0.8, digits = 4)  # Set interval width
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name    mean      sd   lower   upper   rhat ess_bulk ess_tail  sim
#>  cp_1        30.4479 3.60618 25.9129 34.8948 1.0123      135      249 30.0
#>  cp_2        69.7583 0.28765 69.3651 70.1560 0.9998     1424     1529 70.0
#>  match
#>     OK
#>     OK
#> 
#> Population-level parameters:
#>         name    mean      sd   lower   upper   rhat ess_bulk ess_tail  sim
#>  Intercept_1 10.2999 0.71676  9.3981 11.2243 1.0053      389      560 10.0
#>  time_2       0.5323 0.05991  0.4584  0.6112 1.0105      148      278  0.5
#>  Intercept_3  0.6178 1.55829 -1.4041  2.6248 1.0087      193      347  0.0
#>  time_3      -0.2230 0.09218 -0.3389 -0.1041 1.0105      182      388 -0.2
#>  sigma_1      4.0101 0.31144  3.6464  4.4134 1.0008      795      854  4.0
#>  match
#>     OK
#>     OK
#>     OK
#>     OK
#>     OK

# Get the results as a data frame
results = summary(demo_fit)
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  cp_1        30.45 3.606 23.28 37.885 1.01      135      249 30.0    OK
#>  cp_2        69.76 0.288 69.29 70.250 1.00     1424     1529 70.0    OK
#> 
#> Population-level parameters:
#>         name  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 10.30 0.717  8.87 11.670 1.01      389      560 10.0    OK
#>  time_2       0.53 0.060  0.43  0.666 1.01      148      278  0.5    OK
#>  Intercept_3  0.62 1.558 -2.38  3.723 1.01      193      347  0.0    OK
#>  time_3      -0.22 0.092 -0.40 -0.049 1.01      182      388 -0.2    OK
#>  sigma_1      4.01 0.311  3.45  4.657 1.00      795      854  4.0    OK

# Group-level deviations (random effects)
# ranef(my_fit)

# Summarise prior
summary(demo_fit, prior = TRUE)
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name     mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        37.70780 25.43   4.37 92.62 1.00     1807     1681 30.0    OK
#>  cp_2        60.81057 24.43  13.75 96.93 1.00     1813     1728 70.0    OK
#> 
#> Population-level parameters:
#>         name     mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 10.16128 25.31 -28.41 51.34 1.00     1822     1887 10.0    OK
#>  time_2       0.00017  0.21  -0.38  0.41 1.00     1943     1562  0.5      
#>  Intercept_3 10.62569 21.84 -29.30 51.43 1.00     1961     1844  0.0    OK
#>  time_3      -0.00952  0.28  -0.49  0.46 1.00     2196     1800 -0.2    OK
#>  sigma_1     14.12281 16.49   0.45 53.18 1.00     1910     1965  4.0    OK
```
