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
#>         name  mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        24.01 5.016 14.36 33.47    1       92      177 30.0    OK
#>  cp_2        69.91 0.341 69.35 70.48    1     1231     1409 70.0    OK
#> 
#> Population-level parameters:
#>         name  mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.07 0.891  7.35 10.81    1      142      255 10.0    OK
#>  time_2       0.40 0.055  0.31  0.53    1      104      204  0.5    OK
#>  Intercept_3  2.49 1.247 -0.11  4.76    1      205      341  0.0    OK
#>  time_3      -0.28 0.068 -0.41 -0.15    1      201      351 -0.2    OK
#>  sigma_1      3.67 0.272  3.17  4.22    1     1015     1202  4.0    OK
#> 
#> Warning: 5 parameters show poor convergence (rhat > 1.01 or ess_bulk < 400 or ess_tail < 400).
summary(demo_fit, width = 0.8, digits = 4)  # Set interval width
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 2 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name    mean      sd   lower   upper  rhat ess_bulk ess_tail  sim match
#>  cp_1        24.0124 5.01586 17.6723 30.5623 1.019       92      177 30.0    OK
#>  cp_2        69.9116 0.34082 69.4413 70.3846 1.002     1231     1409 70.0    OK
#> 
#> Population-level parameters:
#>         name    mean      sd   lower   upper  rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.0716 0.89137  7.9229 10.1885 1.018      142      255 10.0    OK
#>  time_2       0.4041 0.05485  0.3417  0.4728 1.012      104      204  0.5      
#>  Intercept_3  2.4904 1.24740  0.8677  4.0253 1.007      205      341  0.0      
#>  time_3      -0.2843 0.06790 -0.3659 -0.1917 1.008      201      351 -0.2    OK
#>  sigma_1      3.6675 0.27184  3.3248  4.0329 1.003     1015     1202  4.0    OK
#> 
#> Warning: 5 parameters show poor convergence (rhat > 1.01 or ess_bulk < 400 or ess_tail < 400).

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
#>         name  mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        24.01 5.016 14.36 33.47    1       92      177 30.0    OK
#>  cp_2        69.91 0.341 69.35 70.48    1     1231     1409 70.0    OK
#> 
#> Population-level parameters:
#>         name  mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.07 0.891  7.35 10.81    1      142      255 10.0    OK
#>  time_2       0.40 0.055  0.31  0.53    1      104      204  0.5    OK
#>  Intercept_3  2.49 1.247 -0.11  4.76    1      205      341  0.0    OK
#>  time_3      -0.28 0.068 -0.41 -0.15    1      201      351 -0.2    OK
#>  sigma_1      3.67 0.272  3.17  4.22    1     1015     1202  4.0    OK
#> 
#> Warning: 5 parameters show poor convergence (rhat > 1.01 or ess_bulk < 400 or ess_tail < 400).

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
#>  cp_1        36.03396 26.53   1.26  93.3    1     1807     1681 30.0    OK
#>  cp_2        60.13504 25.49  11.04  97.8    1     1813     1728 70.0    OK
#> 
#> Population-level parameters:
#>         name     mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.29112 22.15 -24.46  45.3    1     1822     1887 10.0    OK
#>  time_2       0.00042  0.52  -0.96   1.0    1     1943     1562  0.5    OK
#>  Intercept_3  9.69748 19.11 -25.24  45.4    1     1961     1844  0.0    OK
#>  time_3      -0.02395  0.70  -1.24   1.2    1     2196     1800 -0.2    OK
#>  sigma_1     12.35746 14.43   0.39  46.5    1     1910     1965  4.0    OK
```
