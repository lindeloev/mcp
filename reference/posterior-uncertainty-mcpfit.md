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
#>              Intercept_1        time_2   Intercept_3       time_3
#> Intercept_1  0.419813848 -0.0070498819  0.0254177297 -0.001137126
#> time_2      -0.007049882  0.0021254541 -0.0004574416  0.000112729
#> Intercept_3  0.025417730 -0.0004574416  1.3800726485 -0.067992377
#> time_3      -0.001137126  0.0001127290 -0.0679923771  0.005365290

# Central posterior intervals for all population-level parameters, or a selection.
confint(demo_fit)
#>                  2.5 %      97.5 %
#> cp_1        28.0610620 34.99837379
#> cp_2        69.4405568 72.76525323
#> Intercept_1  8.7825010 11.32416084
#> time_2       0.4411370  0.61909520
#> Intercept_3 15.3786096 19.96529691
#> time_3      -0.2593416  0.01851412
#> sigma_1      3.3817469  4.46576333
confint(demo_fit, parm = "cp_1")
#>         2.5 %   97.5 %
#> cp_1 28.06106 34.99837
confint(demo_fit, parm = c("Intercept_1", "time_2"), level = 0.8)
#>                  10 %       90 %
#> Intercept_1 9.1995326 10.8527692
#> time_2      0.4708233  0.5877652
confint(demo_fit, prior = TRUE)
#>                  2.5 %     97.5 %
#> cp_1         2.0148120 85.2529521
#> cp_2        17.5944576 98.3678642
#> Intercept_1 -6.7964331 32.9447081
#> time_2      -0.2200525  0.1636851
#> Intercept_3 -4.9707893 31.0920523
#> time_3      -0.1827619  0.2032179
#> sigma_1      0.2832097 26.1626792

# Include change points, residual SDs, group SDs, and AR/MA parameters.
vcov(demo_fit, pars = "all")
#>                     cp_1          cp_2  Intercept_1        time_2   Intercept_3
#> cp_1         2.987920575 -0.0015724881  0.486313050  0.0405205419 -0.0089416257
#> cp_2        -0.001572488  1.0316609631 -0.001714171  0.0007583457 -0.0884850424
#> Intercept_1  0.486313050 -0.0017141706  0.419813848 -0.0070498819  0.0254177297
#> time_2       0.040520542  0.0007583457 -0.007049882  0.0021254541 -0.0004574416
#> Intercept_3 -0.008941626 -0.0884850424  0.025417730 -0.0004574416  1.3800726485
#> time_3       0.002976195  0.0005563016 -0.001137126  0.0001127290 -0.0679923771
#> sigma_1     -0.018688319  0.0073972559  0.003711368 -0.0004579957 -0.0297783267
#>                    time_3       sigma_1
#> cp_1         0.0029761949 -0.0186883193
#> cp_2         0.0005563016  0.0073972559
#> Intercept_1 -0.0011371257  0.0037113677
#> time_2       0.0001127290 -0.0004579957
#> Intercept_3 -0.0679923771 -0.0297783267
#> time_3       0.0053652904  0.0019448920
#> sigma_1      0.0019448920  0.0760206125

# Inspect posterior parameter correlations across the full population model.
# Useful to quickly check identifiability (high correlation). Inspecting
# `bayesplot::mcmc_pairs(as_draws(demo_fit))` is better, though.
vcov(demo_fit, pars = "all", correlation = TRUE)
#>                      cp_1          cp_2  Intercept_1       time_2  Intercept_3
#> cp_1         1.0000000000 -0.0008956418  0.434213197  0.508469442 -0.004403327
#> cp_2        -0.0008956418  1.0000000000 -0.002604697  0.016194699 -0.074156699
#> Intercept_1  0.4342131975 -0.0026046967  1.000000000 -0.236008476  0.033393151
#> time_2       0.5084694420  0.0161946988 -0.236008476  1.000000000 -0.008446151
#> Intercept_3 -0.0044033275 -0.0741566991  0.033393151 -0.008446151  1.000000000
#> time_3       0.0235060809  0.0074773064 -0.023959816  0.033382067 -0.790155619
#> sigma_1     -0.0392120944  0.0264141392  0.020774947 -0.036030467 -0.091935583
#>                   time_3     sigma_1
#> cp_1         0.023506081 -0.03921209
#> cp_2         0.007477306  0.02641414
#> Intercept_1 -0.023959816  0.02077495
#> time_2       0.033382067 -0.03603047
#> Intercept_3 -0.790155619 -0.09193558
#> time_3       1.000000000  0.09630153
#> sigma_1      0.096301533  1.00000000
```
