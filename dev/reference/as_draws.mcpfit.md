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
#> 1   29   72         8.8   0.51          17  -0.37     3.7
#> 2   35   72         9.6   0.56          13  -0.13     4.0
#> 3   31   73         9.9   0.54          16  -0.28     4.2
#> 4   29   72         9.6   0.47          14  -0.18     3.8
#> 5   30   72        10.2   0.50          14  -0.21     3.8
#> 6   35   70        10.5   0.55          14  -0.24     4.3
#> # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

# Other posterior formats are useful in different downstream packages
as_draws_matrix(demo_fit)[1:3, 1:3]  # Matrix of draws by parameter
#> # A draws_matrix: 3 iterations, 1 chains, and 3 variables
#>     variable
#> draw cp_1 cp_2 Intercept_1
#>    1   29   72         8.8
#>    2   35   72         9.6
#>    3   31   73         9.9
as_draws_array(demo_fit)[1:2, , 1:2]  # Iteration-by-chain-by-parameter array
#> # A draws_array: 2 iterations, 2 chains, and 2 variables
#> , , variable = cp_1
#> 
#>          chain
#> iteration  1  2
#>         1 29 33
#>         2 35 31
#> 
#> , , variable = cp_2
#> 
#>          chain
#> iteration  1  2
#>         1 72 71
#>         2 72 71
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
#>    1 28.62198 71.64352    8.784091
#>    2 35.23290 71.99687    9.616996
#>    3 30.65926 72.62149    9.855172
bayesplot::mcmc_trace(demo_fit, pars = c("cp_1", "cp_2"))


# mcp also supports the coda and tidybayes conventions
head(coda::as.mcmc(demo_fit)[[1]])  # First chain as a coda mcmc object
#> Markov Chain Monte Carlo (MCMC) output:
#> Start = 1 
#> End = 7 
#> Thinning interval = 1 
#>          cp_1     cp_2 Intercept_1    time_2 Intercept_3     time_3  sigma_1
#> [1,] 28.62198 71.64352    8.784091 0.5075254    16.77079 -0.3716232 3.677526
#> [2,] 35.23290 71.99687    9.616996 0.5642833    12.86030 -0.1268837 4.032746
#> [3,] 30.65926 72.62149    9.855172 0.5383178    15.70700 -0.2829757 4.201375
#> [4,] 29.41700 72.01096    9.636267 0.4724210    13.63040 -0.1752059 3.803651
#> [5,] 30.46828 72.24478   10.237401 0.5016851    14.31229 -0.2106269 3.804450
#> [6,] 34.84658 70.38014   10.533919 0.5508531    14.21816 -0.2415054 4.314083
#> [7,] 31.28970 71.73028   11.096123 0.5326566    16.06382 -0.4110181 4.119481
if (requireNamespace("tidybayes", quietly = TRUE))
  head(tidybayes::tidy_draws(demo_fit))
#> # A tibble: 6 × 10
#>    cp_1  cp_2 Intercept_1 time_2 Intercept_3 time_3 sigma_1 .chain .iteration
#>   <dbl> <dbl>       <dbl>  <dbl>       <dbl>  <dbl>   <dbl>  <int>      <int>
#> 1  28.6  71.6        8.78  0.508        16.8 -0.372    3.68      1          1
#> 2  35.2  72.0        9.62  0.564        12.9 -0.127    4.03      1          2
#> 3  30.7  72.6        9.86  0.538        15.7 -0.283    4.20      1          3
#> 4  29.4  72.0        9.64  0.472        13.6 -0.175    3.80      1          4
#> 5  30.5  72.2       10.2   0.502        14.3 -0.211    3.80      1          5
#> 6  34.8  70.4       10.5   0.551        14.2 -0.242    4.31      1          6
#> # ℹ 1 more variable: .draw <int>
```
