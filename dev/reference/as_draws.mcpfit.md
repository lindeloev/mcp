# Extract MCMC Draws from `mcpfit` Objects

Extract posterior or prior draws using posterior, coda, or the optional
tidybayes package's S3 generics.

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
#> 1         8.8          17   29   72     3.7   0.51  -0.37
#> 2         9.6          13   35   72     4.0   0.56  -0.13
#> 3         9.9          16   31   73     4.2   0.54  -0.28
#> 4         9.6          14   29   72     3.8   0.47  -0.18
#> 5        10.2          14   30   72     3.8   0.50  -0.21
#> 6        10.5          14   35   70     4.3   0.55  -0.24
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

# Other posterior formats are useful in different downstream packages
as_draws_matrix(demo_fit)[1:3, 1:3]  # Matrix of draws by parameter
#> # A draws_matrix: 3 iterations, 1 chains, and 3 variables
#>     variable
#> draw Intercept_1 Intercept_3 cp_1
#>    1         8.8          17   29
#>    2         9.6          13   35
#>    3         9.9          16   31
as_draws_array(demo_fit)[1:2, , 1:2]  # Iteration-by-chain-by-parameter array
#> # A draws_array: 2 iterations, 2 chains, and 2 variables
#> , , variable = Intercept_1
#> 
#>          chain
#> iteration   1    2
#>         1 8.8 10.2
#>         2 9.6  9.1
#> 
#> , , variable = Intercept_3
#> 
#>          chain
#> iteration  1  2
#>         1 17 16
#>         2 13 15
#> 
as_draws_rvars(demo_fit)[c("cp_1", "cp_2")]  # Random-variable representation
#> # A draws_rvars: 500 iterations, 2 chains, and 2 variables
#> $cp_1: rvar<500,2>[1] mean ± sd:
#> [1] 32 ± 1.7 
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
#>      Intercept_1 Intercept_3     cp_1     cp_2  sigma_1    time_2     time_3
#> [1,]    8.784091    16.77079 28.62198 71.64352 3.677526 0.5075254 -0.3716232
#> [2,]    9.616996    12.86030 35.23290 71.99687 4.032746 0.5642833 -0.1268837
#> [3,]    9.855172    15.70700 30.65926 72.62149 4.201375 0.5383178 -0.2829757
#> [4,]    9.636267    13.63040 29.41700 72.01096 3.803651 0.4724210 -0.1752059
#> [5,]   10.237401    14.31229 30.46828 72.24478 3.804450 0.5016851 -0.2106269
#> [6,]   10.533919    14.21816 34.84658 70.38014 4.314083 0.5508531 -0.2415054
#> [7,]   11.096123    16.06382 31.28970 71.73028 4.119481 0.5326566 -0.4110181
if (requireNamespace("tidybayes", quietly = TRUE))
  head(tidybayes::tidy_draws(demo_fit))
#> # A tibble: 6 × 10
#>   Intercept_1 Intercept_3  cp_1  cp_2 sigma_1 time_2 time_3 .chain .iteration
#>         <dbl>       <dbl> <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <int>      <int>
#> 1        8.78        16.8  28.6  71.6    3.68  0.508 -0.372      1          1
#> 2        9.62        12.9  35.2  72.0    4.03  0.564 -0.127      1          2
#> 3        9.86        15.7  30.7  72.6    4.20  0.538 -0.283      1          3
#> 4        9.64        13.6  29.4  72.0    3.80  0.472 -0.175      1          4
#> 5       10.2         14.3  30.5  72.2    3.80  0.502 -0.211      1          5
#> 6       10.5         14.2  34.8  70.4    4.31  0.551 -0.242      1          6
#> # ℹ 1 more variable: .draw <int>
```
