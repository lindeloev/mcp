# Prior and Posterior Covariance and Central Intervals for `mcpfit` Objects

Summarise the joint and marginal uncertainty of population-level model
parameters using posterior or prior draws.

## Usage

``` r
# S3 method for class 'mcpfit'
vcov(object, correlation = FALSE, pars = NULL, dpar = "mu", prior = FALSE, ...)

# S3 method for class 'mcpfit'
confint(object, parm, level = 0.95, prior = FALSE, ...)
```

## Arguments

- object:

  An `mcpfit` object.

- correlation:

  Return the correlation matrix instead of the covariance matrix?

- pars:

  Optional names of population-level parameters to extract, or `"all"`
  for all population-level parameters.

- dpar:

  Distributional parameter(s) to select when `pars = NULL`.

- prior:

  Logical. Use prior draws instead of posterior draws?

- ...:

  Must be empty. Reserved for future use.

- parm:

  Optional names or positions of population-level parameters to include
  in the intervals.

- level:

  Width of the central interval.

## Value

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns a covariance or
correlation matrix. [`confint()`](https://rdrr.io/r/stats/confint.html)
returns a two-column matrix of central intervals.

## Examples

``` r
# Posterior covariance of the primary-response coefficients, matching fixef().
vcov(demo_fit)
#>               Intercept_1        time_2   Intercept_3        time_3
#> Intercept_1  0.4416343541 -0.0072905718  0.0304178993 -0.0002884934
#> time_2      -0.0072905718  0.0021712914 -0.0004292144 -0.0001262524
#> Intercept_3  0.0304178993 -0.0004292144  1.8679145063 -0.0929213047
#> time_3      -0.0002884934 -0.0001262524 -0.0929213047  0.0067363802

# Central posterior intervals for all population-level parameters, or a selection.
confint(demo_fit)
#>                  2.5 %     97.5 %
#> cp_1        28.2437124 34.9139556
#> cp_2        69.4788981 72.7884742
#> Intercept_1  8.4927383 11.2226971
#> time_2       0.4591914  0.6340292
#> Intercept_3 11.2996947 16.6996581
#> time_3      -0.3618760 -0.0443138
#> sigma_1      3.3661781  4.4420735
confint(demo_fit, parm = "cp_1")
#>         2.5 %   97.5 %
#> cp_1 28.24371 34.91396
confint(demo_fit, parm = c("Intercept_1", "time_2"), level = 0.8)
#>                 10 %       90 %
#> Intercept_1 9.075140 10.7753049
#> time_2      0.483584  0.6034823
confint(demo_fit, prior = TRUE)
#>                  2.5 %    97.5 %
#> cp_1          2.100901 82.850839
#> cp_2         15.205694 98.320851
#> Intercept_1 -21.277525 50.113926
#> time_2       -1.206786  1.182860
#> Intercept_3 -20.720652 47.144012
#> time_3       -1.210176  1.247309
#> sigma_1       0.154470 22.495286

# Include change points, residual SDs, group SDs, and AR/MA parameters.
vcov(demo_fit, pars = "all")
#>                     cp_1          cp_2   Intercept_1        time_2
#> cp_1         3.041576158 -1.471207e-02  0.4877127091  4.569380e-02
#> cp_2        -0.014712065  1.029368e+00 -0.0070031106  5.243563e-05
#> Intercept_1  0.487712709 -7.003111e-03  0.4416343541 -7.290572e-03
#> time_2       0.045693804  5.243563e-05 -0.0072905718  2.171291e-03
#> Intercept_3  0.022820566 -1.683385e-01  0.0304178993 -4.292144e-04
#> time_3      -0.001816144 -2.897673e-03 -0.0002884934 -1.262524e-04
#> sigma_1      0.001773736 -1.118579e-02  0.0031484330  6.143777e-05
#>               Intercept_3        time_3       sigma_1
#> cp_1         0.0228205659 -0.0018161441  1.773736e-03
#> cp_2        -0.1683384888 -0.0028976730 -1.118579e-02
#> Intercept_1  0.0304178993 -0.0002884934  3.148433e-03
#> time_2      -0.0004292144 -0.0001262524  6.143777e-05
#> Intercept_3  1.8679145063 -0.0929213047 -1.723748e-02
#> time_3      -0.0929213047  0.0067363802  9.624823e-04
#> sigma_1     -0.0172374806  0.0009624823  7.856290e-02

# Inspect posterior parameter correlations across the full population model.
# Useful to quickly check identifiability (high correlation). Inspecting
# `bayesplot::mcmc_pairs(as_draws(demo_fit))` is better, though.
vcov(demo_fit, pars = "all", correlation = TRUE)
#>                     cp_1         cp_2  Intercept_1       time_2  Intercept_3
#> cp_1         1.000000000 -0.008314554  0.420807334  0.562275175  0.009574113
#> cp_2        -0.008314554  1.000000000 -0.010386619  0.001109129 -0.121400251
#> Intercept_1  0.420807334 -0.010386619  1.000000000 -0.235435045  0.033490320
#> time_2       0.562275175  0.001109129 -0.235435045  1.000000000 -0.006739641
#> Intercept_3  0.009574113 -0.121400251  0.033490320 -0.006739641  1.000000000
#> time_3      -0.012687839 -0.034797716 -0.005289214 -0.033011667 -0.828368365
#> sigma_1      0.003628531 -0.039334393  0.016902630  0.004704003 -0.044997308
#>                   time_3      sigma_1
#> cp_1        -0.012687839  0.003628531
#> cp_2        -0.034797716 -0.039334393
#> Intercept_1 -0.005289214  0.016902630
#> time_2      -0.033011667  0.004704003
#> Intercept_3 -0.828368365 -0.044997308
#> time_3       1.000000000  0.041837978
#> sigma_1      0.041837978  1.000000000
```
