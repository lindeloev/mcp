# Get the JAGS Model, Data, and Initial Values of an mcp Model

**\[experimental\]**

## Usage

``` r
as_jags(fit, prior = FALSE)
```

## Arguments

- fit:

  An
  [mcpfit](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object, e.g., from `mcp(..., sample = FALSE)`.

- prior:

  Logical. Return data and initial values for sampling the prior instead
  of the posterior, as in `mcp(..., sample = "prior")`.

## Value

An `mcpjags` object, which prints the JAGS code. It is a list with:

- `code`: String. The JAGS model.

- `data`: List with the data.

- `inits`: List with initial values for each chain, including JAGS
  random-number seeds. These are the initial values that
  [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) used
  to sample `fit`. The seeds are random for fits with `sample = FALSE`
  and `seed = NULL`.

- `variable.names`: Character vector. The nodes that
  [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md)
  monitors, i.e., `mcp_pars(fit)$parameter`.

- `modules`: Character vector. The JAGS modules that
  [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) loads,
  e.g., `"glm"`.

- `levels`: Named list with the levels of each grouping variable. JAGS
  indexes them as 1, 2, etc. in this order, e.g., `x_1_id[2]` is the
  deviation for `levels$id[2]`.

## Details

Returns everything needed to sample an `mcp` model directly in JAGS.
This gives more flexibility to modify the JAGS code or to monitor other
nodes. The downside is that `mcp` methods
([`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
[`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md),
[`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
etc.) cannot use draws sampled this way. Use
[`mcpfamily()`](https://lindeloev.github.io/mcp/dev/reference/mcpfamily.md)
to fit other likelihoods in `mcp`.

## Examples

``` r
data = mcp_example_data("demo")
fit = mcp(list(response ~ 1, ~ 0 + time, ~ 1 + time), data, sample = FALSE)
jags = as_jags(fit)
jags  # Show the JAGS code. Also in jags$code
#> model {
#>   # Priors for change points
#>   cp_0 = min(time)
#>   cp_1 ~ dunif(cp_0, cp_3)  # Flat Dirichlet on the segment lengths
#>   cp_2 ~ dunif(cp_0, cp_3)  # Flat Dirichlet on the segment lengths
#>   cp_3 = max(time)
#>   cp_order_ ~ dbern(step(cp_2 - cp_1))  # Ordered change points
#> 
#>   # Priors for population-level effects
#>   Intercept_1 ~ dnorm(13.825120048974153, 1/(17.382500430759464)^2)  # Mean intercept (rstanarm default)
#>   time_2_end_ ~ dnorm(Intercept_1, 1/(0.62123709509255609*(cp_2-cp_1))^2)  # Autoscaled mean coefficient (rstanarm default); sampled as the value at the segment end
#>   time_2 = (time_2_end_ - (Intercept_1)) / (cp_2 - cp_1)
#>   Intercept_3_end_ ~ dnorm(13.825120048974153 + time_3 * (cp_3 - cp_2), 1/(17.382500430759464)^2)  # Mean intercept (rstanarm default); sampled as the value at the segment end
#>   Intercept_3 = Intercept_3_end_ - (time_3 * (cp_3 - cp_2))
#>   time_3 ~ dnorm(0, 1/(0.62123709509255609)^2)  # Autoscaled mean coefficient (rstanarm default)
#>   sigma_1 ~ dt(0, 1/(5.7)^2, 3) T(0, )  # Positive residual SD calibrated on the response scale
#> 
#>   # Model and likelihood
#>   for (i_ in 1:length(time)) {
#>     # time measured from each segment's change point
#>     x_local_1_[i_] = min(time[i_], cp_1)
#>     x_local_2_[i_] = min(time[i_], cp_2) - cp_1
#>     x_local_3_[i_] = min(time[i_], cp_3) - cp_2
#>     
#>     # Formula for mu
#>     link_mu_[i_] =
#>       (time[i_] >= cp_0) * (time[i_] < cp_2) * inprod(rhs_matrix_[i_, c(1)], c(Intercept_1)) * 1 + 
#>       (time[i_] >= cp_1) * (time[i_] < cp_2) * inprod(rhs_matrix_[i_, c(2)], c(time_2)) * x_local_2_[i_] + 
#>       (time[i_] >= cp_2) * inprod(rhs_matrix_[i_, c(3)], c(Intercept_3)) * 1 + 
#>       (time[i_] >= cp_2) * inprod(rhs_matrix_[i_, c(4)], c(time_3)) * x_local_3_[i_]
#>     
#>     # Formula for sigma
#>     link_sigma_[i_] =
#>       (time[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(5)], c(sigma_1)) * 1
#> 
#>     # Likelihood for family = gaussian()
#>     mu_[i_] = link_mu_[i_]
#>     sigma_[i_] = max(1e-05, link_sigma_[i_])
#>     response[i_] ~ dnorm(mu_[i_], 1 / sigma_[i_]^2)
#>   }
#> }
#> 
#> # Also contains $data, $inits (3 chains), $variable.names, $modules, and $levels. See ?as_jags

# \donttest{
# After making changes to jags$code or other elements, sample it using rjags
for (module in jags$modules)
  rjags::load.module(module)
#> module glm loaded
model = rjags::jags.model(textConnection(jags$code), jags$data, jags$inits, n.chains = length(jags$inits))
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 101
#>    Unobserved stochastic nodes: 7
#>    Total graph size: 2150
#> 
#> Initializing model
#> 
draws = rjags::coda.samples(model, jags$variable.names, n.iter = 1000)
summary(draws)  # This calls coda:::summary.mcmc.list() - not mcp:::summary.mcpfit()
#> 
#> Iterations = 1001:2000
#> Thinning interval = 1 
#> Number of chains = 3 
#> Sample size per chain = 1000 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                Mean      SD  Naive SE Time-series SE
#> Intercept_1  9.9291 0.64874 0.0118442       0.014390
#> Intercept_3 13.9308 1.34393 0.0245367       0.045865
#> cp_1        31.6979 1.75381 0.0320200       0.051989
#> cp_2        71.1182 1.00442 0.0183382       0.030091
#> sigma_1      3.8741 0.29155 0.0053230       0.007309
#> time_2       0.5465 0.04690 0.0008564       0.001134
#> time_3      -0.2022 0.08118 0.0014821       0.003915
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%     25%     50%     75%    97.5%
#> Intercept_1  8.6643  9.4843  9.9311 10.3633 11.21971
#> Intercept_3 11.3221 13.0147 13.9202 14.8497 16.53749
#> cp_1        28.4005 30.5060 31.6739 32.8211 35.11139
#> cp_2        69.4631 70.2414 71.1240 71.9814 72.78092
#> sigma_1      3.3495  3.6782  3.8540  4.0538  4.49282
#> time_2       0.4592  0.5144  0.5461  0.5772  0.64008
#> time_3      -0.3539 -0.2587 -0.2032 -0.1475 -0.03601
#> 
# }
```
