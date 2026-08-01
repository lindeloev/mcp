# Summarise mcpfit objects

Summarise parameter estimates and model diagnostics.

## Usage

``` r
# S3 method for class 'mcpfit'
summary(object, width = 0.95, digits = 2, prior = FALSE, verbose = FALSE, ...)

# S3 method for class 'mcpfit'
fixef(object, width = 0.95, prior = FALSE, verbose = FALSE, ...)

# S3 method for class 'mcpfit'
ranef(object, width = 0.95, prior = FALSE, verbose = FALSE, ...)

# S3 method for class 'mcpfit'
print(x, ...)
```

## Arguments

- object:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object.

- width:

  Float. The width of the central posterior interval (between 0 and 1).

- digits:

  a non-null value for digits specifies the minimum number of
  significant digits to be printed in values. The default, NULL, uses
  getOption("digits"). (For the interpretation for complex numbers see
  signif.) Non-integer values will be rounded down, and only values
  greater than or equal to 1 and no greater than 22 are accepted.

- prior:

  TRUE/FALSE. Summarise prior instead of posterior?

- verbose:

  Logical. Include the `segment` and `dpar` columns. Defaults to `FALSE`
  for a compact, v0.3.4-compatible summary.

- ...:

  Currently ignored

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
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
  `"ar1"`, `"ma1"`, etc.) the parameter belongs to.

- `mean` is the posterior mean

- `lower` and `upper` are the bounds of the central posterior interval
  given in `width`.

- `Rhat` is the rank-normalized split-Rhat convergence diagnostic.

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

- `fixef(mcpfit)`: Fixed (population-level) effects of `mcpfit`.

- `ranef(mcpfit)`: Random (varying) effects of `mcpfit`.

- `print(mcpfit)`: Print the posterior summary of an
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# Typical usage
summary(demo_fit)
#> Family: gaussian(link = 'identity')
#> Iterations: 3000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Population-level parameters:
#>         name match  sim  mean lower upper Rhat ess_bulk ess_tail
#>         cp_1       30.0 10.96   8.8  13.4  2.7        4       19
#>         cp_2       70.0 15.63  14.4  16.3  3.1        3       12
#>  Intercept_1    OK 10.0  8.35   5.5  12.8  2.5        4       22
#>       time_2    OK  0.5  0.15  -2.2   1.5  2.8        3       11
#>  Intercept_3        0.0 11.71  10.3  14.1  3.4        3       11
#>       time_3    OK -0.2 -0.31  -2.8   1.3  3.0        3       14
#>      sigma_1        4.0  9.17   5.7  11.7  3.1        3       12
#> 
#> Warning: 7 parameters show poor convergence (Rhat > 1.01 or ESS < 400).
summary(demo_fit, width = 0.8, digits = 4)  # Set interval width
#> Family: gaussian(link = 'identity')
#> Iterations: 3000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Population-level parameters:
#>         name match  sim    mean  lower  upper  Rhat ess_bulk ess_tail
#>         cp_1       30.0 10.9614  8.784 13.381 2.655        4       19
#>         cp_2       70.0 15.6302 14.440 16.250 3.054        3       12
#>  Intercept_1    OK 10.0  8.3516  5.539 12.843 2.540        4       22
#>       time_2    OK  0.5  0.1523 -2.142  1.501 2.826        3       11
#>  Intercept_3        0.0 11.7093 10.282 14.095 3.404        3       11
#>       time_3    OK -0.2 -0.3131 -2.749  1.299 2.951        3       14
#>      sigma_1        4.0  9.1750  5.749 11.717 3.106        3       12
#> 
#> Warning: 7 parameters show poor convergence (Rhat > 1.01 or ESS < 400).

# Get the results as a data frame
results = summary(demo_fit)
#> Family: gaussian(link = 'identity')
#> Iterations: 3000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Population-level parameters:
#>         name match  sim  mean lower upper Rhat ess_bulk ess_tail
#>         cp_1       30.0 10.96   8.8  13.4  2.7        4       19
#>         cp_2       70.0 15.63  14.4  16.3  3.1        3       12
#>  Intercept_1    OK 10.0  8.35   5.5  12.8  2.5        4       22
#>       time_2    OK  0.5  0.15  -2.2   1.5  2.8        3       11
#>  Intercept_3        0.0 11.71  10.3  14.1  3.4        3       11
#>       time_3    OK -0.2 -0.31  -2.8   1.3  3.0        3       14
#>      sigma_1        4.0  9.17   5.7  11.7  3.1        3       12
#> 
#> Warning: 7 parameters show poor convergence (Rhat > 1.01 or ESS < 400).

# Varying (random) effects
# ranef(my_fit)

# Summarise prior
summary(demo_fit, prior = TRUE)
#> Family: gaussian(link = 'identity')
#> Iterations: 3000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Population-level parameters:
#>         name match  sim   mean  lower upper Rhat ess_bulk ess_tail
#>         cp_1    OK 30.0 37.799   4.23  91.5    1     2673     2481
#>         cp_2    OK 70.0 60.947  14.47  97.0    1     2459     2898
#>  Intercept_1    OK 10.0 10.601 -28.21  49.9    1     2811     2902
#>       time_2    OK  0.5 -0.014  -1.18   1.2    1     2638     2750
#>  Intercept_3    OK  0.0 10.405 -30.56  50.6    1     2872     2863
#>       time_3    OK -0.2  0.002  -1.30   1.3    1     2716     2883
#>      sigma_1    OK  4.0 14.225   0.42  51.8    1     2935     2990
```
