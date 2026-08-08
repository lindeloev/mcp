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
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
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
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# Typical usage
summary(demo_fit)
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Population-level parameters:
#>         name match  sim  mean lower upper Rhat ess_bulk ess_tail
#>         cp_1    OK 30.0 23.21 10.08 32.58    1       95      147
#>         cp_2    OK 70.0 69.93 69.35 70.48    1     1618     2087
#>  Intercept_1    OK 10.0  8.92  6.66 10.62    1      143      200
#>       time_2    OK  0.5  0.40  0.30  0.51    1      115      339
#>  Intercept_3    OK  0.0  2.26 -0.17  4.64    1      299      451
#>       time_3    OK -0.2 -0.27 -0.40 -0.14    1      308      440
#>      sigma_1    OK  4.0  3.68  3.20  4.30    1     1478     1496
#> 
#> Warning: 5 parameters show poor convergence (Rhat > 1.01 or ESS < 400).
summary(demo_fit, width = 0.8, digits = 4)  # Set interval width
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Population-level parameters:
#>         name match  sim    mean   lower   upper   Rhat ess_bulk ess_tail
#>         cp_1    OK 30.0 23.2088 16.4267 30.1189 1.0438       95      147
#>         cp_2    OK 70.0 69.9283 69.4451 70.3969 1.0008     1618     2087
#>  Intercept_1    OK 10.0  8.9217  7.7206 10.0853 1.0324      143      200
#>       time_2        0.5  0.3991  0.3334  0.4707 1.0357      115      339
#>  Intercept_3        0.0  2.2600  0.6608  3.8577 1.0129      299      451
#>       time_3    OK -0.2 -0.2707 -0.3589 -0.1821 1.0162      308      440
#>      sigma_1    OK  4.0  3.6818  3.3514  4.0441 0.9999     1478     1496
#> 
#> Warning: 5 parameters show poor convergence (Rhat > 1.01 or ESS < 400).

# Get the results as a data frame
results = summary(demo_fit)
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Population-level parameters:
#>         name match  sim  mean lower upper Rhat ess_bulk ess_tail
#>         cp_1    OK 30.0 23.21 10.08 32.58    1       95      147
#>         cp_2    OK 70.0 69.93 69.35 70.48    1     1618     2087
#>  Intercept_1    OK 10.0  8.92  6.66 10.62    1      143      200
#>       time_2    OK  0.5  0.40  0.30  0.51    1      115      339
#>  Intercept_3    OK  0.0  2.26 -0.17  4.64    1      299      451
#>       time_3    OK -0.2 -0.27 -0.40 -0.14    1      308      440
#>      sigma_1    OK  4.0  3.68  3.20  4.30    1     1478     1496
#> 
#> Warning: 5 parameters show poor convergence (Rhat > 1.01 or ESS < 400).

# Varying (random) effects
# ranef(my_fit)

# Summarise prior
summary(demo_fit, prior = TRUE)
#> Family: gaussian(link = 'identity')
#> Iterations: 1000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Population-level parameters:
#>         name match  sim    mean  lower upper Rhat ess_bulk ess_tail
#>         cp_1    OK 30.0 36.4153   1.33  93.0    1     2658     2473
#>         cp_2    OK 70.0 60.9086  11.83  97.8    1     2904     2675
#>  Intercept_1    OK 10.0  9.2724 -25.89  43.8    1     3013     2649
#>       time_2    OK  0.5  0.0115  -0.98   1.1    1     2859     2989
#>  Intercept_3    OK  0.0  9.9600 -27.24  46.1    1     2934     2875
#>       time_3    OK -0.2  0.0066  -1.12   1.1    1     2873     2905
#>      sigma_1    OK  4.0 12.5162   0.39  45.9    1     3059     2988
```
