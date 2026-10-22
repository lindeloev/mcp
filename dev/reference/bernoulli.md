# Bernoulli Family for mcp

Binary responses (0 or 1). See [Binomial and
Bernoulli](https://lindeloev.github.io/mcp/articles/binomial.html) for
worked examples.

## Usage

``` r
bernoulli(link = "logit")
```

## Arguments

- link:

  Link function.

## Examples

``` r
# Fit a binary-response model with a probit link
data = data.frame(time = 1:6, y = c(0, 0, 0, 1, 1, 1))
fit = mcp(list(y ~ 1), data, family = bernoulli(link = "probit"), par_x = "time", sample = FALSE)
mcp_pars(fit)$parameter  # Show the parameters of the fitted Bernoulli model
#> [1] "Intercept_1"
```
