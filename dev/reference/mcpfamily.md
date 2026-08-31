# Create or Test Objects of Class "mcpfamily"

**\[experimental\]**

## Usage

``` r
mcpfamily(x)

# S3 method for class 'mcpfamily'
format(x, ...)

# S3 method for class 'mcpfamily'
print(x, ...)

is.mcpfamily(x)
```

## Arguments

- x:

  A family object, e.g., `binomial(link = "identity")`.

- ...:

  Must be empty. Reserved for future use.

## Details

Converts standard R family objects into `mcpfamily` objects used
internally by `mcp`. Supported family and link combinations include:

- `gaussian(link = "identity")` or `gaussian(link = "log")`

- `binomial(link = "logit")`, `binomial(link = "probit")`, or
  `binomial(link = "identity")`

- `bernoulli(link = "logit")`, `bernoulli(link = "probit")`, or
  `bernoulli(link = "identity")`

- `poisson(link = "log")` or `poisson(link = "identity")`

- `negbinomial(link = "log", link_shape = "log")`

Note: `mcpfamily` objects are shipped with mcp - there is not (yet)
support for user-supplied families.

## Methods (by generic)

- `format(mcpfamily)`: Format an `mcpfamily` object.

- `print(mcpfamily)`: Print an `mcpfamily` object.

## Functions

- `is.mcpfamily()`: Checks whether x is an `mcpfamily`.

## See also

[`family`](https://rdrr.io/r/stats/family.html)

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# mcp() converts supported standard family objects automatically
data = data.frame(time = 1:6, y = exp(seq(0, 1, length.out = 6)))
my_normal = mcpfamily(stats::gaussian(link = "log"))
fit = mcp(list(y ~ 1), data, family = my_normal, par_x = "time", sample = FALSE)
family(fit)  # Show the mcp family retained in the fit
#> Family: gaussian
#> Links: mu = log; sigma = identity

# The converted object can also be inspected directly
mcpfamily(stats::binomial())$dpars  # Show its distributional parameters
#> [1] "mu"
```
