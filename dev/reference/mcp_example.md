# Run example models

Run example models

## Usage

``` r
mcp_example(name, sample = "post", plot = TRUE)

mcp_example_data(name)
```

## Arguments

- name:

  Name of the example. One of:

  - `"demo"`: Two change points between intercepts and joined/disjoined
    slopes.

  - `"intercepts"`: An intercept-only change point.

  - `"multiple"`: Multiple regression with categorical predictors and
    interactions.

  - `"binomial"`: Binomial with two change points. Much like `"demo"` on
    a logit scale.

  - `"group_mu"`: Group-level predictor deviations (random
    intercepts/slopes) across a change point.

  - `"group_cp"`: Group-level change-point deviations (random effects).

  - `"missing"`: Missing data imputation (NAs in response variable y).

  - `"quadratic"`: A change point to a quadratic segment where there is
    no data.

  - `"ar"`: One change point in autoregressive residuals (the `ar1`
    dpar)

  - `"sigma"`: A change in "sigma" dpar, including a slope on sigma.

- sample:

  One of

  - `"post"`: Sample the posterior.

  - `"prior"`: Sample only the prior. Plots, summaries, etc. will use
    the prior. This is useful for prior predictive checks.

  - `"both"`: Sample both prior and posterior. Plots, summaries, etc.
    will default to using the posterior. The prior only has effect when
    doing Savage-Dickey density ratios in
    [`hypothesis`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md).

  - `"none"` or `FALSE`: Do not sample. Returns an mcpfit object without
    sample. This is useful if you only want to check prior strings
    (`fit$prior`), the JAGS model (`fit$jags_code`), etc.

- plot:

  Logical. Plot the fitted example? Requires sample != "none".

## Value

An `mcpfit`, enriched with an `$example_code` field. It contains the
code to reproduce the data and the fit.

## Functions

- `mcp_example_data()`: Conveniently get simulated data only.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# \donttest{
fit = mcp_example("multiple")

print(fit$example_code) # See how the data was simulated
#> # Define model
#> model = list(
#>   y ~ 1 + x:group + z,
#>   ~ 1 + x + group,
#>   ~ 0 + I(x^2)
#> )
#> 
#> # Simulate data
#> set.seed(42)
#> data = data.frame(
#>   x = 1:120,
#>   group = rep(c('A', 'B', 'C', 'D'), 30),
#>   z = rnorm(120, mean = 1:120, sd = 25),
#>   y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
#> )
#> empty = mcp(model, data, sample = FALSE, par_x = 'x')
#> data$y = empty$simulate(empty, data,
#>   cp_1 = 69.5,
#>   cp_2 = 99.5,
#> 
#>   Intercept_1 = 10,
#>   z_1 = 0.2,
#>   xgroupA_1 = -0.75,
#>   xgroupB_1 = -0.25,
#>   xgroupC_1 = 0.25,
#>   xgroupD_1 = 0.75,
#> 
#>   Intercept_2 = 10,
#>   x_2 = -1,
#>   groupB_2 = 15,
#>   groupC_2 = 30,
#>   groupD_2 = 45,
#> 
#>   xE2_3 = 0.2,
#> 
#>   sigma_1 = 5
#> )
#> 
#> # Run sampling
#> fit = mcp(model, data, par_x = 'x', iter = 10000, sample = sample, seed = 42)
#> 
#> # Illustrative plot
#> if (plot) {
#>   set.seed(42)
#>   print(plot(fit, color_by = 'group') + ggplot2::labs(title = 'plot(fit, color_by = "group")'))
#> }

# Without sampling
empty = mcp_example("binomial", sample = "none", plot = FALSE)
print(empty)
#> Family: binomial
#> Links: mu = logit
#> Segments:
#>   1: y | trials(N) ~ 1
#>   2: y | trials(N) ~ 1 ~ 0 + x
#>   3: y | trials(N) ~ 1 ~ 1 + x
#> 
#> No draws. Nothing to summarise.
print(empty$example_code)
#> # Define model
#> model = list(
#>   y | trials(N) ~ 1,  # constant success probability
#>   ~ 0 + x,  # joined changing success probability
#>   ~ 1 + x  # disjoined changing success probability
#> )
#> 
#> # Simulate data
#> set.seed(42)
#> data = data.frame(
#>   x = 1:100,
#>   N = base::sample(15, 25, replace=TRUE),
#>   y = 0.  # Numeric placeholder that is valid for every sampled trial count.
#> )
#> empty = mcp(model, data, family = binomial(), sample = FALSE)
#> data$y = empty$simulate(empty, data,
#>   cp_1 = 29.5,
#>   cp_2 = 69.5,
#>   Intercept_1 = 1.5,
#>   Intercept_3 = -1,
#>   x_2 = -0.15,
#>   x_3 = 0.05
#> )
#> 
#> # Run sampling
#> fit = mcp(model, data, family = binomial(), iter = 4000, sample = sample, seed = 42)
#> 
#> # Illustrative plot
#> if (plot) {
#>   set.seed(42)
#>   print(plot(fit, q_fit = TRUE) + ggplot2::labs(title = 'plot(fit, q_fit = TRUE)'))
#> }

# Now sample this model
fit2 = mcp(empty$model, empty$data, family = empty$family, warmup = 2000, iter = 6000, seed = 42)
plot(fit2)

# }
```
