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
#> 1        10.1         4.7   24   70     4.1   0.36  -0.33
#> 2        10.1         4.2   25   70     3.8   0.37  -0.33
#> 3        10.1         2.5   25   70     3.8   0.39  -0.26
#> 4         9.9         1.6   25   70     3.8   0.37  -0.26
#> 5         9.5         3.2   25   70     3.5   0.39  -0.36
#> 6         9.1         3.0   25   70     3.5   0.42  -0.25
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

# Other posterior formats are useful in different downstream packages
as_draws_matrix(demo_fit)[1:3, 1:3]  # Matrix of draws by parameter
#> # A draws_matrix: 3 iterations, 1 chains, and 3 variables
#>     variable
#> draw Intercept_1 Intercept_3 cp_1
#>    1          10         4.7   24
#>    2          10         4.2   25
#>    3          10         2.5   25
as_draws_array(demo_fit)[1:2, , 1:2]  # Iteration-by-chain-by-parameter array
#> # A draws_array: 2 iterations, 2 chains, and 2 variables
#> , , variable = Intercept_1
#> 
#>          chain
#> iteration  1   2
#>         1 10 8.7
#>         2 10 8.7
#> 
#> , , variable = Intercept_3
#> 
#>          chain
#> iteration   1   2
#>         1 4.7 2.2
#>         2 4.2 2.2
#> 
as_draws_rvars(demo_fit)[c("cp_1", "cp_2")]  # Random-variable representation
#> # A draws_rvars: 1000 iterations, 2 chains, and 2 variables
#> $cp_1: rvar<1000,2>[1] mean ± sd:
#> [1] 24 ± 5 
#> 
#> $cp_2: rvar<1000,2>[1] mean ± sd:
#> [1] 70 ± 0.34 
#> 

# mcp also supports the coda and tidybayes conventions
head(coda::as.mcmc(demo_fit)[[1]])  # First chain as a coda mcmc object
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 3001 
#> End = 3007 
#> Thinning interval = 1 
#>      Intercept_1 Intercept_3     cp_1     cp_2  sigma_1    time_2     time_3
#> [1,]   10.122805    4.733476 23.62971 70.13133 4.079635 0.3602996 -0.3309533
#> [2,]   10.073020    4.211827 25.26966 70.03862 3.839004 0.3735152 -0.3343376
#> [3,]   10.146927    2.519549 25.07811 69.54561 3.846829 0.3870755 -0.2639302
#> [4,]    9.917544    1.576157 24.69675 69.56728 3.794215 0.3747316 -0.2649606
#> [5,]    9.548849    3.185002 25.32573 69.73441 3.520012 0.3934620 -0.3631516
#> [6,]    9.063937    2.963639 25.00800 70.31785 3.529331 0.4162720 -0.2506363
#> [7,]    9.349596    2.441074 25.60126 69.35883 3.434344 0.3968397 -0.2671036
head(tidybayes::tidy_draws(demo_fit))  # Tidybayes-compatible draw data
#> # A draws_df: 6 iterations, 1 chains, and 7 variables
#>   Intercept_1 Intercept_3 cp_1 cp_2 sigma_1 time_2 time_3
#> 1        10.1         4.7   24   70     4.1   0.36  -0.33
#> 2        10.1         4.2   25   70     3.8   0.37  -0.33
#> 3        10.1         2.5   25   70     3.8   0.39  -0.26
#> 4         9.9         1.6   25   70     3.8   0.37  -0.26
#> 5         9.5         3.2   25   70     3.5   0.39  -0.36
#> 6         9.1         3.0   25   70     3.5   0.42  -0.25
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}
```
