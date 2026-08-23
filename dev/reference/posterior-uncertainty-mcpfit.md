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
#>              Intercept_1        time_2  Intercept_3        time_3
#> Intercept_1  0.794544888  0.0151052438  0.005135117 -0.0002433870
#> time_2       0.015105244  0.0030089358 -0.004954620  0.0002716614
#> Intercept_3  0.005135117 -0.0049546204  1.556002444 -0.0728512322
#> time_3      -0.000243387  0.0002716614 -0.072851232  0.0046103775

# Central posterior intervals for all population-level parameters, or a selection.
confint(demo_fit)
#>                  2.5 %     97.5 %
#> cp_1        14.3594277 33.4693818
#> cp_2        69.3476085 70.4764874
#> Intercept_1  7.3470487 10.8055330
#> time_2       0.3118097  0.5260004
#> Intercept_3 -0.1117829  4.7570735
#> time_3      -0.4094870 -0.1464992
#> sigma_1      3.1746954  4.2236295
confint(demo_fit, parm = "cp_1")
#>         2.5 %   97.5 %
#> cp_1 14.35943 33.46938
confint(demo_fit, parm = c("Intercept_1", "time_2"), level = 0.8)
#>                  10 %       90 %
#> Intercept_1 7.9228753 10.1885454
#> time_2      0.3417202  0.4728081

# Include change points, residual SDs, group SDs, and AR/MA parameters.
vcov(demo_fit, pars = "all")
#>                     cp_1          cp_2   Intercept_1        time_2  Intercept_3
#> cp_1        25.158871662  6.406286e-02  3.0619949618  0.2295109170 -0.196183242
#> cp_2         0.064062858  1.161571e-01  0.0077159033  0.0006905697 -0.035459392
#> Intercept_1  3.061994962  7.715903e-03  0.7945448880  0.0151052438  0.005135117
#> time_2       0.229510917  6.905697e-04  0.0151052438  0.0030089358 -0.004954620
#> Intercept_3 -0.196183242 -3.545939e-02  0.0051351173 -0.0049546204  1.556002444
#> time_3       0.010912540  7.938059e-05 -0.0002433870  0.0002716614 -0.072851232
#> sigma_1      0.005186872 -2.287294e-03 -0.0009893186  0.0001566627 -0.009764095
#>                    time_3       sigma_1
#> cp_1         1.091254e-02  0.0051868718
#> cp_2         7.938059e-05 -0.0022872943
#> Intercept_1 -2.433870e-04 -0.0009893186
#> time_2       2.716614e-04  0.0001566627
#> Intercept_3 -7.285123e-02 -0.0097640946
#> time_3       4.610378e-03  0.0005880054
#> sigma_1      5.880054e-04  0.0738989038

# Inspect posterior parameter correlations across the full population model.
# Useful to quickly check identifiability (high correlation). Inspecting
# `bayesplot::mcmc_pairs(as_draws(demo_fit))` is better, though.
vcov(demo_fit, pars = "all", correlation = TRUE)
#>                     cp_1        cp_2  Intercept_1      time_2  Intercept_3
#> cp_1         1.000000000  0.03747467  0.684856643  0.83416375 -0.031355317
#> cp_2         0.037474672  1.00000000  0.025398301  0.03693841 -0.083407146
#> Intercept_1  0.684856643  0.02539830  1.000000000  0.30893141  0.004618341
#> time_2       0.834163746  0.03693841  0.308931408  1.00000000 -0.072410041
#> Intercept_3 -0.031355317 -0.08340715  0.004618341 -0.07241004  1.000000000
#> time_3       0.032041439  0.00343023 -0.004021332  0.07293789 -0.860128861
#> sigma_1      0.003804004 -0.02468767 -0.004082799  0.01050607 -0.028794387
#>                   time_3      sigma_1
#> cp_1         0.032041439  0.003804004
#> cp_2         0.003430230 -0.024687674
#> Intercept_1 -0.004021332 -0.004082799
#> time_2       0.072937894  0.010506068
#> Intercept_3 -0.860128861 -0.028794387
#> time_3       1.000000000  0.031856207
#> sigma_1      0.031856207  1.000000000
```
