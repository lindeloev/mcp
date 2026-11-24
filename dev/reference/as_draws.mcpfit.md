# Extract MCMC Draws from `mcpfit` Objects

Extract posterior or prior draws using posterior, coda, or the optional
tidybayes package's S3 generics. Base R's
[`as.matrix()`](https://rdrr.io/r/base/matrix.html),
[`as.array()`](https://rdrr.io/r/base/array.html), and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
plain draws, like in rstanarm and brms, so bayesplot functions work
directly on the fit.

## Usage

``` r
# S3 method for class 'mcpfit'
as_draws(x, prior = FALSE, ...)

as_draws(x, ...)

as_draws_df(x, ...)

as_draws_array(x, ...)

as_draws_matrix(x, ...)

as_draws_rvars(x, ...)

# S3 method for class 'mcpfit'
as.matrix(x, prior = FALSE, ...)

# S3 method for class 'mcpfit'
as.array(x, prior = FALSE, ...)

# S3 method for class 'mcpfit'
as.data.frame(x, row.names = NULL, optional = FALSE, prior = FALSE, ...)
```

## Arguments

- x:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object.

- prior:

  Logical. Extract prior draws (`TRUE`) instead of posterior draws
  (`FALSE`)? Errors if the requested draws are unavailable.

- ...:

  Passed to posterior or tidybayes format conversion functions.

- row.names, optional:

  Unused. For compatibility with
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

## Value

A posterior `draws` object or a coda `mcmc.list` object.
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) returns a
draws-by-parameter matrix,
[`as.array()`](https://rdrr.io/r/base/array.html) an
iteration-by-chain-by-parameter array, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) a
draws-by-parameter data frame. The data frame has no chain information,
so use `as_draws_df()` for chain-aware diagnostics such as
[`posterior::rhat()`](https://mc-stan.org/posterior/reference/rhat.html).

## Examples

``` r
# Default posterior draws, with one row per iteration and chain
draws = as_draws(demo_fit)  # Return a posterior::draws object
head(as_draws_df(demo_fit))  # Convert draws to a data frame
#> # A draws_df: 6 iterations, 1 chains, and 7 variables
#>   cp_1 cp_2 Intercept_1 time_2 Intercept_3 time_3 sigma_1
#> 1   32   70        10.2   0.58          11 -0.058     4.5
#> 2   29   72         8.5   0.58          14 -0.197     4.1
#> 3   29   70         9.5   0.52          14 -0.199     4.0
#> 4   33   70        11.0   0.52          15 -0.277     4.1
#> 5   34   70        10.2   0.61          15 -0.276     4.2
#> 6   30   70         9.9   0.54          15 -0.226     4.0
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

# Other posterior formats are useful in different downstream packages
as_draws_matrix(demo_fit)[1:3, 1:3]  # Matrix of draws by parameter
#> # A draws_matrix: 3 iterations, 1 chains, and 3 variables
#>     variable
#> draw cp_1 cp_2 Intercept_1
#>    1   32   70        10.2
#>    2   29   72         8.5
#>    3   29   70         9.5
as_draws_array(demo_fit)[1:2, , 1:2]  # Iteration-by-chain-by-parameter array
#> # A draws_array: 2 iterations, 2 chains, and 2 variables
#> , , variable = cp_1
#> 
#>          chain
#> iteration  1  2
#>         1 32 32
#>         2 29 30
#> 
#> , , variable = cp_2
#> 
#>          chain
#> iteration  1  2
#>         1 70 71
#>         2 72 69
#> 
as_draws_rvars(demo_fit)[c("cp_1", "cp_2")]  # Random-variable representation
#> # A draws_rvars: 500 iterations, 2 chains, and 2 variables
#> $cp_1: rvar<500,2>[1] mean ± sd:
#> [1] 32 ± 1.7 
#> 
#> $cp_2: rvar<500,2>[1] mean ± sd:
#> [1] 71 ± 1 
#> 

# Base R formats, e.g., for bayesplot
as.matrix(demo_fit)[1:3, 1:3]  # Plain matrix of draws by parameter
#>     variable
#> draw     cp_1     cp_2 Intercept_1
#>    1 32.48487 69.61262   10.200329
#>    2 29.41653 72.31654    8.483119
#>    3 29.01040 69.64645    9.492935
bayesplot::mcmc_trace(demo_fit, pars = c("cp_1", "cp_2"))


# mcp also supports the coda and tidybayes conventions
head(coda::as.mcmc(demo_fit)[[1]])  # First chain as a coda mcmc object
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1 
#> End = 7 
#> Thinning interval = 1 
#>          cp_1     cp_2 Intercept_1    time_2 Intercept_3       time_3  sigma_1
#> [1,] 32.48487 69.61262   10.200329 0.5821892    10.79263 -0.057669450 4.486435
#> [2,] 29.41653 72.31654    8.483119 0.5820990    13.86030 -0.197292607 4.135718
#> [3,] 29.01040 69.64645    9.492935 0.5209982    14.03013 -0.199308936 3.986916
#> [4,] 32.84258 69.67823   10.969884 0.5213880    14.78298 -0.276674271 4.051947
#> [5,] 33.87082 69.52184   10.191440 0.6092685    15.15793 -0.276271007 4.163503
#> [6,] 30.41488 70.16472    9.937738 0.5385329    14.98238 -0.226361690 4.036436
#> [7,] 29.69202 69.44907    9.475047 0.4710464    11.16928 -0.005983112 3.700439
if (requireNamespace("tidybayes", quietly = TRUE))
  head(tidybayes::tidy_draws(demo_fit))
#> # A tibble: 6 × 10
#>    cp_1  cp_2 Intercept_1 time_2 Intercept_3  time_3 sigma_1 .chain .iteration .draw
#>   <dbl> <dbl>       <dbl>  <dbl>       <dbl>   <dbl>   <dbl>  <int>      <int> <int>
#> 1  32.5  69.6       10.2   0.582        10.8 -0.0577    4.49      1          1     1
#> 2  29.4  72.3        8.48  0.582        13.9 -0.197     4.14      1          2     2
#> 3  29.0  69.6        9.49  0.521        14.0 -0.199     3.99      1          3     3
#> 4  32.8  69.7       11.0   0.521        14.8 -0.277     4.05      1          4     4
#> 5  33.9  69.5       10.2   0.609        15.2 -0.276     4.16      1          5     5
#> 6  30.4  70.2        9.94  0.539        15.0 -0.226     4.04      1          6     6
```
