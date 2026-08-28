# Posterior Covariance and Central Intervals for `mcpfit` Objects

Summarise the joint and marginal posterior uncertainty of
population-level model parameters.

## Usage

``` r
# S3 method for class 'mcpfit'
vcov(object, correlation = FALSE, pars = NULL, dpar = "mu", ...)

# S3 method for class 'mcpfit'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An `mcpfit` object.

- correlation:

  Return the posterior correlation matrix instead of the covariance
  matrix?

- pars:

  Optional names of population-level parameters to extract, or `"all"`
  for all population-level parameters.

- dpar:

  Distributional parameter(s) to select when `pars = NULL`.

- ...:

  Currently unused.

- parm:

  Optional names or positions of population-level parameters to include
  in the intervals.

- level:

  Width of the central posterior interval.

## Value

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns a posterior
covariance or correlation matrix.
[`confint()`](https://rdrr.io/r/stats/confint.html) returns a two-column
matrix of central posterior intervals.

## Examples

``` r
# Posterior covariance of the primary-response coefficients, matching fixef().
vcov(demo_fit)
#>               Intercept_1        time_2  Intercept_3        time_3
#> Intercept_1  0.5137397261  0.0076750258 -0.032563602  0.0004612561
#> time_2       0.0076750258  0.0035887962  0.007317392 -0.0005024094
#> Intercept_3 -0.0325636015  0.0073173924  2.428266606 -0.1270849180
#> time_3       0.0004612561 -0.0005024094 -0.127084918  0.0084964259

# Central posterior intervals for all population-level parameters, or a selection.
confint(demo_fit)
#>                  2.5 %      97.5 %
#> cp_1        23.2818627 37.88476201
#> cp_2        69.2948696 70.25002448
#> Intercept_1  8.8747825 11.66981091
#> time_2       0.4273044  0.66607328
#> Intercept_3 -2.3834621  3.72290985
#> time_3      -0.4021069 -0.04872518
#> sigma_1      3.4488478  4.65679839
confint(demo_fit, parm = "cp_1")
#>         2.5 %   97.5 %
#> cp_1 23.28186 37.88476
confint(demo_fit, parm = c("Intercept_1", "time_2"), level = 0.8)
#>                  10 %       90 %
#> Intercept_1 9.3981254 11.2242586
#> time_2      0.4584488  0.6111839

# Include change points, residual SDs, group SDs, and AR/MA parameters.
vcov(demo_fit, pars = "all")
#>                    cp_1          cp_2   Intercept_1        time_2  Intercept_3
#> cp_1        13.00453024 -0.0078184799  1.3685958360  0.1817068751  0.335986912
#> cp_2        -0.00781848  0.0827428878 -0.0009984594 -0.0002535124  0.003969372
#> Intercept_1  1.36859584 -0.0009984594  0.5137397261  0.0076750258 -0.032563602
#> time_2       0.18170688 -0.0002535124  0.0076750258  0.0035887962  0.007317392
#> Intercept_3  0.33598691  0.0039693719 -0.0325636015  0.0073173924  2.428266606
#> time_3      -0.02591721 -0.0015059827  0.0004612561 -0.0005024094 -0.127084918
#> sigma_1     -0.02475858 -0.0017215457  0.0041714465 -0.0004574783 -0.022651261
#>                    time_3       sigma_1
#> cp_1        -0.0259172125 -0.0247585782
#> cp_2        -0.0015059827 -0.0017215457
#> Intercept_1  0.0004612561  0.0041714465
#> time_2      -0.0005024094 -0.0004574783
#> Intercept_3 -0.1270849180 -0.0226512615
#> time_3       0.0084964259  0.0021237310
#> sigma_1      0.0021237310  0.0969936201

# Inspect posterior parameter correlations across the full population model.
# Useful to quickly check identifiability (high correlation). Inspecting
# `bayesplot::mcmc_pairs(as_draws(demo_fit))` is better, though.
vcov(demo_fit, pars = "all", correlation = TRUE)
#>                     cp_1         cp_2  Intercept_1      time_2  Intercept_3
#> cp_1         1.000000000 -0.007537193  0.529488223  0.84110387  0.059789764
#> cp_2        -0.007537193  1.000000000 -0.004842766 -0.01471158  0.008855402
#> Intercept_1  0.529488223 -0.004842766  1.000000000  0.17874499 -0.029154975
#> time_2       0.841103865 -0.014711583  0.178744987  1.00000000  0.078385144
#> Intercept_3  0.059789764  0.008855402 -0.029154975  0.07838514  1.000000000
#> time_3      -0.077969153 -0.056798452  0.006981559 -0.09098405 -0.884764506
#> sigma_1     -0.022044831 -0.019216824  0.018687169 -0.02452023 -0.046673742
#>                   time_3     sigma_1
#> cp_1        -0.077969153 -0.02204483
#> cp_2        -0.056798452 -0.01921682
#> Intercept_1  0.006981559  0.01868717
#> time_2      -0.090984050 -0.02452023
#> Intercept_3 -0.884764506 -0.04667374
#> time_3       1.000000000  0.07397923
#> sigma_1      0.073979225  1.00000000
```
