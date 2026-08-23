# Summarise priors used by an mcp model

**\[experimental\]**

## Usage

``` r
prior_summary(fit, verbose = FALSE)
```

## Arguments

- fit:

  An `mcpfit` object.

- verbose:

  Logical. Include rule, description, source, and kind.

## Value

A tibble with one row per model parameter, ordered and labeled the same
way as
[`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md):
change points first, then `mu`, then the other distributional
parameters, then `ar`/`ma` components - each with `segment` and `dpar`
columns.

## Details

Shows the effective, resolved prior distributions on the familiar
SD/scale parameterization rather than JAGS precision. Symbolic
expressions in bounds (e.g. `min(x)` or `max(x)`) may be retained in
symbolic form in compact output, while `verbose = TRUE` displays the
underlying rule, description, source, and kind.

## Examples

``` r
prior_summary(demo_fit)  # Show the effective priors and bounds
#> # A tibble: 7 × 5
#>   parameter   segment dpar  prior                                         bounds
#>   <chr>         <int> <chr> <chr>                                         <chr> 
#> 1 cp_1              2 cp    student_t(df = 1, location = 0.02388966, sca… [min(…
#> 2 cp_2              3 cp    student_t(df = 1, location = 0.02388966, sca… [cp_1…
#> 3 Intercept_1       1 mu    student_t(df = 3, location = 9.5, scale = 11… none  
#> 4 time_2            2 mu    student_t(df = 3, location = 0, scale = 0.33… none  
#> 5 Intercept_3       3 mu    student_t(df = 3, location = 9.5, scale = 11… none  
#> 6 time_3            3 mu    student_t(df = 3, location = 0, scale = 0.33… none  
#> 7 sigma_1           1 sigma student_t(df = 3, location = 0, scale = 11.2) [0, I…
prior_summary(demo_fit, verbose = TRUE)  # Include their rules and sources
#> # A tibble: 7 × 9
#>   parameter   segment dpar  prior          bounds rule  description source kind 
#>   <chr>         <int> <chr> <chr>          <chr>  <chr> <chr>       <chr>  <chr>
#> 1 cp_1              2 cp    student_t(df … [min(… stud… Regularizi… defau… dist…
#> 2 cp_2              3 cp    student_t(df … [cp_1… stud… Regularizi… defau… dist…
#> 3 Intercept_1       1 mu    student_t(df … none   stud… Robustly c… defau… dist…
#> 4 time_2            2 mu    student_t(df … none   stud… Regularizi… defau… dist…
#> 5 Intercept_3       3 mu    student_t(df … none   stud… Robustly c… defau… dist…
#> 6 time_3            3 mu    student_t(df … none   stud… Regularizi… defau… dist…
#> 7 sigma_1           1 sigma student_t(df … [0, I… stud… Positive r… defau… dist…
```
