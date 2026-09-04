# Extract MCMC Draws from `mcpfit` Objects

Extract posterior or prior draws using posterior, tidybayes, or coda S3
generics.

## Usage

``` r
# S3 method for class 'mcpfit'
as_draws(x, prior = FALSE, ...)

as_draws(x, ...)

as_draws_df(x, ...)

as_draws_array(x, ...)

as_draws_matrix(x, ...)

as_draws_rvars(x, ...)
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

## Value

A posterior `draws` object or a coda `mcmc.list` object.

## Examples

``` r
# Default posterior draws, with one row per iteration and chain
draws = as_draws(demo_fit)  # Return a posterior::draws object
head(as_draws_df(demo_fit))  # Convert draws to a data frame
#> # A draws_df: 6 iterations, 1 chains, and 7 variables
#>   Intercept_1 Intercept_3 cp_1 cp_2 sigma_1 time_2 time_3
#> 1        10.2          16   31   72     4.5   0.51  0.077
#> 2         9.9          16   31   70     3.6   0.51 -0.047
#> 3         9.4          17   31   71     3.6   0.51 -0.078
#> 4        10.8          17   32   71     3.4   0.54 -0.075
#> 5         9.3          19   33   73     3.8   0.59 -0.178
#> 6         9.3          17   30   72     4.0   0.54 -0.113
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

# Other posterior formats are useful in different downstream packages
as_draws_matrix(demo_fit)[1:3, 1:3]  # Matrix of draws by parameter
#> # A draws_matrix: 3 iterations, 1 chains, and 3 variables
#>     variable
#> draw Intercept_1 Intercept_3 cp_1
#>    1        10.2          16   31
#>    2         9.9          16   31
#>    3         9.4          17   31
as_draws_array(demo_fit)[1:2, , 1:2]  # Iteration-by-chain-by-parameter array
#> # A draws_array: 2 iterations, 2 chains, and 2 variables
#> , , variable = Intercept_1
#> 
#>          chain
#> iteration    1    2
#>         1 10.2  9.6
#>         2  9.9 10.1
#> 
#> , , variable = Intercept_3
#> 
#>          chain
#> iteration  1  2
#>         1 16 19
#>         2 16 16
#> 
as_draws_rvars(demo_fit)[c("cp_1", "cp_2")]  # Random-variable representation
#> # A draws_rvars: 500 iterations, 2 chains, and 2 variables
#> $cp_1: rvar<500,2>[1] mean ± sd:
#> [1] 31 ± 1.7 
#> 
#> $cp_2: rvar<500,2>[1] mean ± sd:
#> [1] 71 ± 1 
#> 

# mcp also supports the coda and tidybayes conventions
head(coda::as.mcmc(demo_fit)[[1]])  # First chain as a coda mcmc object
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1 
#> End = 7 
#> Thinning interval = 1 
#>      Intercept_1 Intercept_3     cp_1     cp_2  sigma_1    time_2      time_3
#> [1,]   10.191209    15.69228 30.84290 72.12381 4.459326 0.5126160  0.07732289
#> [2,]    9.865942    16.35426 30.53939 69.99791 3.571110 0.5118716 -0.04739203
#> [3,]    9.380313    16.53452 30.63572 71.06167 3.615166 0.5141171 -0.07802441
#> [4,]   10.789719    16.97402 31.93417 70.89540 3.437456 0.5390244 -0.07489507
#> [5,]    9.323356    18.69243 33.12817 72.57656 3.843492 0.5946020 -0.17761014
#> [6,]    9.305478    16.87526 30.01061 72.06966 3.968703 0.5380811 -0.11323448
#> [7,]   10.021265    15.98717 32.82439 72.41413 4.029903 0.5880111 -0.04399810
head(tidybayes::tidy_draws(demo_fit))  # Tidybayes-compatible draw data
#> # A draws_df: 6 iterations, 1 chains, and 7 variables
#>   Intercept_1 Intercept_3 cp_1 cp_2 sigma_1 time_2 time_3
#> 1        10.2          16   31   72     4.5   0.51  0.077
#> 2         9.9          16   31   70     3.6   0.51 -0.047
#> 3         9.4          17   31   71     3.6   0.51 -0.078
#> 4        10.8          17   32   71     3.4   0.54 -0.075
#> 5         9.3          19   33   73     3.8   0.59 -0.178
#> 6         9.3          17   30   72     4.0   0.54 -0.113
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}
```
