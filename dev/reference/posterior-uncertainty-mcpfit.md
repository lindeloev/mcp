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

  Optional names of population-level parameters to extract, or
  `"population"` for all population-level parameters.

- dpar:

  Distributional parameter(s) to select when `pars = NULL`, e.g.,
  `"mu"`, `"sigma"`, or `"ar1"`.

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

## Details

As in `lme4`, [`vcov()`](https://rdrr.io/r/stats/vcov.html) covers the
coefficients of
[fixef()](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
(by default those of `mu`), while
[`confint()`](https://rdrr.io/r/stats/confint.html) covers all
population-level parameters, including change points.

## Examples

``` r
# Posterior covariance of the primary-response coefficients, matching fixef().
vcov(demo_fit)
#>              Intercept_1        time_2  Intercept_3        time_3
#> Intercept_1  0.441730996 -0.0076944262 -0.008694904  0.0010577664
#> time_2      -0.007694426  0.0022602415 -0.002418374  0.0002482322
#> Intercept_3 -0.008694904 -0.0024183736  1.920748300 -0.0976086443
#> time_3       0.001057766  0.0002482322 -0.097608644  0.0070004302

# Central posterior intervals for all population-level parameters, or a selection.
confint(demo_fit)
#>                  2.5 %     97.5 %
#> cp_1        28.2323892 34.8833815
#> cp_2        69.4796848 72.7703641
#> Intercept_1  8.6811979 11.2449667
#> time_2       0.4539101  0.6431522
#> Intercept_3 11.1116538 16.6473670
#> time_3      -0.3729092 -0.0477629
#> sigma_1      3.3443922  4.5399852
confint(demo_fit, parm = "cp_1")
#>         2.5 %   97.5 %
#> cp_1 28.23239 34.88338
confint(demo_fit, parm = c("Intercept_1", "time_2"), level = 0.8)
#>                  10 %       90 %
#> Intercept_1 9.1112810 10.8454130
#> time_2      0.4825704  0.6025843
confint(demo_fit, prior = TRUE)
#>                   2.5 %    97.5 %
#> cp_1          1.9860251 85.773990
#> cp_2         16.3799745 98.109020
#> Intercept_1 -22.8308241 48.063439
#> time_2       -1.2259134  1.275298
#> Intercept_3 -19.8533241 47.565743
#> time_3       -1.2460525  1.243472
#> sigma_1       0.1277497 25.805782

# Include change points, residual SDs, group-level SDs, and AR/MA parameters.
vcov(demo_fit, pars = "population")
#>                     cp_1         cp_2  Intercept_1        time_2  Intercept_3        time_3       sigma_1
#> cp_1         2.860236828 -0.023006037  0.479259027  0.0450919199 -0.113947879  7.548661e-03  1.694093e-02
#> cp_2        -0.023006037  0.982126142  0.026726008 -0.0018991431 -0.277033079  4.117322e-03 -3.602333e-03
#> Intercept_1  0.479259027  0.026726008  0.441730996 -0.0076944262 -0.008694904  1.057766e-03  8.129028e-03
#> time_2       0.045091920 -0.001899143 -0.007694426  0.0022602415 -0.002418374  2.482322e-04  3.602231e-04
#> Intercept_3 -0.113947879 -0.277033079 -0.008694904 -0.0024183736  1.920748300 -9.760864e-02 -1.340728e-02
#> time_3       0.007548661  0.004117322  0.001057766  0.0002482322 -0.097608644  7.000430e-03 -2.727354e-05
#> sigma_1      0.016940926 -0.003602333  0.008129028  0.0003602231 -0.013407276 -2.727354e-05  8.593321e-02

# Inspect posterior parameter correlations across the full population model.
# Useful to quickly check identifiability (high correlation). Inspecting
# `bayesplot::mcmc_pairs(as_draws(demo_fit))` is better, though.
vcov(demo_fit, pars = "population", correlation = TRUE)
#>                    cp_1        cp_2  Intercept_1      time_2  Intercept_3       time_3      sigma_1
#> cp_1         1.00000000 -0.01372642  0.426373660  0.56081588 -0.048614970  0.053346578  0.034170837
#> cp_2        -0.01372642  1.00000000  0.040576214 -0.04030850 -0.201703189  0.049655673 -0.012399945
#> Intercept_1  0.42637366  0.04057621  1.000000000 -0.24351172 -0.009439535  0.019021663  0.041723350
#> time_2       0.56081588 -0.04030850 -0.243511725  1.00000000 -0.036703775  0.062404847  0.025847209
#> Intercept_3 -0.04861497 -0.20170319 -0.009439535 -0.03670378  1.000000000 -0.841764654 -0.033000824
#> time_3       0.05334658  0.04965567  0.019021663  0.06240485 -0.841764654  1.000000000 -0.001111984
#> sigma_1      0.03417084 -0.01239994  0.041723350  0.02584721 -0.033000824 -0.001111984  1.000000000
```
