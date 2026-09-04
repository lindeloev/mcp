# Families and link functions

To do GLM, you specify a family and link function in the conventional
way:

``` r

fit = mcp(model, data = df, family = gaussian("log"))
fit = mcp(model, data = df, family = binomial("identity"))
```

Currently supported families and link functions include:

- Default: `gaussian(link = "identity")`. Alternatives:
  `gaussian(link = "log")`
- Default: `binomial(link = "logit")`. Alternatives:
  `binomial(link = "probit")`, `binomial(link = "identity")`
- Default: `bernoulli(link = "logit")`. Alternatives:
  `bernoulli(link = "probit")`, `bernoulli(link = "identity")`
- Default: `poisson(link = "log")`. Alternatives:
  `poisson(link = "identity")`
- Default: `negbinomial(link = "log", link_shape = "log")`

See the “GLM” menu above for more details on using GLM with `mcp`. The
hardest part about adding new families/links is specifying default
priors. Read more in [the article on mcp
priors](https://lindeloev.github.io/mcp/articles/priors.md).

## On-default link functions

Some link functions are *default* in GLMs for good reasons: they are
computationally convenient and interpretable. A non-default link can
imply impossible values, at which point `mcp` errors. For example, a
`bernoulli("identity")` model `prob ~ 1 + x` places a slope directly on
`P(y = TRUE)`, which can easily fall below 0 or exceed 1. Informative
[priors](https://lindeloev.github.io/mcp/articles/priors.md), for
example through truncation, can prevent the sampler from visiting
invalid values.

In short: think carefully and proceed at your own risk.
