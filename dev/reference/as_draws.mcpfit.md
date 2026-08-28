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
  (`FALSE`)?

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
#> 1         8.4       2.111   20   70     4.0   0.39 -0.298
#> 2         8.9       0.875   20   69     3.7   0.40 -0.205
#> 3         9.7      -0.905   20   70     3.9   0.40 -0.105
#> 4         9.6      -1.210   21   69     4.9   0.42 -0.095
#> 5         9.4      -0.168   22   70     4.5   0.40 -0.199
#> 6        11.2      -0.015   24   70     4.5   0.40 -0.195
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

# Other posterior formats are useful in different downstream packages
as_draws_matrix(demo_fit)[1:3, 1:3]  # Matrix of draws by parameter
#> # A draws_matrix: 3 iterations, 1 chains, and 3 variables
#>     variable
#> draw Intercept_1 Intercept_3 cp_1
#>    1         8.4        2.11   20
#>    2         8.9        0.88   20
#>    3         9.7       -0.91   20
as_draws_array(demo_fit)[1:2, , 1:2]  # Iteration-by-chain-by-parameter array
#> # A draws_array: 2 iterations, 2 chains, and 2 variables
#> , , variable = Intercept_1
#> 
#>          chain
#> iteration   1  2
#>         1 8.4 10
#>         2 8.9 11
#> 
#> , , variable = Intercept_3
#> 
#>          chain
#> iteration    1   2
#>         1 2.11 1.7
#>         2 0.88 1.4
#> 
as_draws_rvars(demo_fit)[c("cp_1", "cp_2")]  # Random-variable representation
#> # A draws_rvars: 1000 iterations, 2 chains, and 2 variables
#> $cp_1: rvar<1000,2>[1] mean ± sd:
#> [1] 30 ± 3.6 
#> 
#> $cp_2: rvar<1000,2>[1] mean ± sd:
#> [1] 70 ± 0.29 
#> 

# mcp also supports the coda and tidybayes conventions
head(coda::as.mcmc(demo_fit)[[1]])  # First chain as a coda mcmc object
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 3001 
#> End = 3007 
#> Thinning interval = 1 
#>      Intercept_1 Intercept_3     cp_1     cp_2  sigma_1    time_2      time_3
#> [1,]    8.440762  2.11057600 19.84137 70.18666 3.996256 0.3947496 -0.29764104
#> [2,]    8.858730  0.87505090 19.55044 69.49629 3.705049 0.3998132 -0.20542023
#> [3,]    9.729883 -0.90545278 20.18012 69.94807 3.907430 0.4014265 -0.10543576
#> [4,]    9.566698 -1.20979125 20.64717 69.33957 4.853139 0.4213045 -0.09496154
#> [5,]    9.393083 -0.16835228 22.46844 70.21179 4.541024 0.4040311 -0.19902490
#> [6,]   11.202711 -0.01532687 23.86199 70.20102 4.512374 0.4011385 -0.19509195
#> [7,]   10.922353  0.46771522 23.84820 69.91665 4.587211 0.4174504 -0.12629247
head(tidybayes::tidy_draws(demo_fit))  # Tidybayes-compatible draw data
#> # A draws_df: 6 iterations, 1 chains, and 7 variables
#>   Intercept_1 Intercept_3 cp_1 cp_2 sigma_1 time_2 time_3
#> 1         8.4       2.111   20   70     4.0   0.39 -0.298
#> 2         8.9       0.875   20   69     3.7   0.40 -0.205
#> 3         9.7      -0.905   20   70     3.9   0.40 -0.105
#> 4         9.6      -1.210   21   69     4.9   0.42 -0.095
#> 5         9.4      -0.168   22   70     4.5   0.40 -0.199
#> 6        11.2      -0.015   24   70     4.5   0.40 -0.195
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}
```
