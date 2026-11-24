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
#> Intercept_1  0.432601772 -0.0072880609 -0.025366661  0.0018197182
#> time_2      -0.007288061  0.0021951341  0.002298114 -0.0001429837
#> Intercept_3 -0.025366661  0.0022981145  1.708083869 -0.0851944535
#> time_3       0.001819718 -0.0001429837 -0.085194454  0.0061900325

# Central posterior intervals for all population-level parameters, or a selection.
confint(demo_fit)
#>                  2.5 %      97.5 %
#> cp_1        28.4436593 35.12834318
#> cp_2        69.4651798 72.76760453
#> Intercept_1  8.6041580 11.22531492
#> time_2       0.4519378  0.63916695
#> Intercept_3 11.3599241 16.71961788
#> time_3      -0.3606648 -0.04875881
#> sigma_1      3.3921200  4.51374442
confint(demo_fit, parm = "cp_1")
#>         2.5 %   97.5 %
#> cp_1 28.44366 35.12834
confint(demo_fit, parm = c("Intercept_1", "time_2"), level = 0.8)
#>                  10 %      90 %
#> Intercept_1 9.0800790 10.762627
#> time_2      0.4820554  0.604029
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
#>                     cp_1          cp_2   Intercept_1        time_2  Intercept_3        time_3       sigma_1
#> cp_1         2.915379907 -0.0014760645  0.4843835253  0.0439545558  0.088473453 -2.815580e-03 -3.435019e-03
#> cp_2        -0.001476065  1.0066137713 -0.0162539711  0.0008672107 -0.234489443  5.894607e-04 -2.355474e-03
#> Intercept_1  0.484383525 -0.0162539711  0.4326017716 -0.0072880609 -0.025366661  1.819718e-03  3.708866e-04
#> time_2       0.043954556  0.0008672107 -0.0072880609  0.0021951341  0.002298114 -1.429837e-04  2.700886e-04
#> Intercept_3  0.088473453 -0.2344894431 -0.0253666614  0.0022981145  1.708083869 -8.519445e-02 -9.734890e-03
#> time_3      -0.002815580  0.0005894607  0.0018197182 -0.0001429837 -0.085194454  6.190032e-03  7.044311e-05
#> sigma_1     -0.003435019 -0.0023554739  0.0003708866  0.0002700886 -0.009734890  7.044311e-05  7.685546e-02

# Inspect posterior parameter correlations across the full population model.
# Useful to quickly check identifiability (high correlation). Inspecting
# `bayesplot::mcmc_pairs(as_draws(demo_fit))` is better, though.
vcov(demo_fit, pars = "population", correlation = TRUE)
#>                     cp_1         cp_2  Intercept_1      time_2 Intercept_3       time_3      sigma_1
#> cp_1         1.000000000 -0.000861641  0.431318055  0.54944698  0.03964704 -0.020959154 -0.007256784
#> cp_2        -0.000861641  1.000000000 -0.024631105  0.01844857 -0.17882872  0.007467530 -0.008468559
#> Intercept_1  0.431318055 -0.024631105  1.000000000 -0.23650336 -0.02950968  0.035165235  0.002034041
#> time_2       0.549446978  0.018448567 -0.236503360  1.00000000  0.03753069 -0.038789094  0.020794012
#> Intercept_3  0.039647043 -0.178828719 -0.029509685  0.03753069  1.00000000 -0.828533860 -0.026868231
#> time_3      -0.020959154  0.007467530  0.035165235 -0.03878909 -0.82853386  1.000000000  0.003229644
#> sigma_1     -0.007256784 -0.008468559  0.002034041  0.02079401 -0.02686823  0.003229644  1.000000000
```
