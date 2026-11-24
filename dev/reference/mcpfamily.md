# Define a Family for mcp

**\[experimental\]**

## Usage

``` r
mcpfamily(
  family,
  links,
  lower = NULL,
  jags,
  r = list(),
  default_prior = NULL,
  response = list(),
  garma = NULL,
  glm_family = NULL,
  floor = NULL
)
```

## Arguments

- family:

  String. Name the family for printing, e.g., `"student"`.

- links:

  Named character vector with the link function of each distributional
  parameter (dpar), e.g., `c(mu = "identity", sigma = "log")`. Supported
  links: `"identity"`, `"log"`, `"logit"`, and `"probit"`. The dpar `mu`
  is required.

- lower:

  Named numeric vector. The exclusive lower bound of dpars, e.g.,
  `c(sigma = 0)`. Priors on constant dpars (e.g., `sigma_1`) are
  truncated to it.

- jags:

  Named list with the JAGS likelihood:

  - `likelihood`: A function taking a `context` list and returning JAGS
    code for observation `i_` (a string or a character vector of lines).
    See examples. `context` contains:

    - `y`: The response node, e.g., `"y[i_]"`.

    - `dpar(name)`: A dpar on the response scale, e.g.,
      `context$dpar("sigma")` returns `"sigma_[i_]"`.

    - `aux(name, default)`: A response auxiliary column, e.g.,
      `context$aux("weights", "1")`. Returns `default` if the model does
      not include it.

    - `local(name)`: A node for intermediate values, e.g.,
      `context$local("prob")` returns `"prob_[i_]"`.

    - `zeros_trick(neg_log_density)`: For likelihoods that JAGS does not
      have. Takes the negative log-density of `y` as a string. Add a
      constant to keep it non-negative. Applies weights but does not
      support missing responses. See the [Gumbel
      example](https://lindeloev.github.io/mcp/articles/custom_families.html#likelihoods-that-jags-does-not-have).

  - `glm_module`: When to load the JAGS glm module: `"never"` (default),
    `"group_effects"`, or `"always"`.

- r:

  Named list of R functions. They take `dpars`, a list with a
  rows-by-draws matrix for each dpar on its response scale, and `data`,
  a list with the response auxiliary columns, e.g., `data$trials`. Like
  [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html), they must
  recycle `data` over the draws (columns). See examples and [Custom
  families
  article](https://lindeloev.github.io/mcp/articles/custom_families.html)
  for applied examples.

  - `log_lik(y, dpars, data)`: The log-density of each response,
    including normalizing constants. When supplied, `mcp` checks
    `jags$likelihood` against it after sampling. Required by
    [`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
    [`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
    and
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).
    For weighted models, the observation weights are applied by `mcp`.

  - `epred(dpars, data, rate = FALSE)`: The expected response. Used by
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md),
    etc.

  - `simulate(n, dpars, data, rate = FALSE)`: Random responses. Used by
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`pp_check()`](https://lindeloev.github.io/mcp/dev/reference/pp_check.md),
    `fit$simulate()`, etc.

  - `variance(dpars, data, rate = FALSE)`: The variance of the response.
    Used for exact
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    SDs. If missing, the SD of the simulated responses is used.

  - `cdf(q, dpars, data)` and `quantile(p, dpars, data)`: Used for exact
    predictive intervals. If missing, the quantiles of the simulated
    responses are used, which is less precise.

- default_prior:

  `NULL` or a data frame with columns `dpar`, `par_type` (`"Intercept"`,
  `"dummy"`, or `"slope"`), and `prior`, e.g.,
  `"dnorm(0, 2.5 / predictor_scale())"`. Each parameter gets the row
  matching its `dpar` and `par_type` in
  [`mcp_pars()`](https://lindeloev.github.io/mcp/dev/reference/mcp_pars.md).
  Optional columns:

  - `group_sd_prior`: Prior for the SD of group-level effects.

  - `description`: Shown in `prior_summary(fit, verbose = TRUE)`.

  - `condition`: one of

    - `"always"`: (default)

    - `"included"`: for dpar terms explicitly included in formula, e.g.,
      `sigma_Intercept_1` in `sigma(1)`.

    - `"not_included"` for dpar terms not included in any formula, e.g.,
      `sigma_1` implicit when
      [`sigma()`](https://rdrr.io/r/stats/sigma.html) term is absent.

  Priors use the syntax of `mcp(..., prior)` (see
  [mcp-priors](https://lindeloev.github.io/mcp/dev/reference/mcp-priors.md)).
  `.y` and `.x` are the response and the change-point variable, e.g.,
  `"dnorm(mean(.y), 2.5 * sd(.y))"`. `predictor_scale()` is the SD of
  each coefficient's predictor, for autoscaling.
  [`prior_summary()`](https://lindeloev.github.io/mcp/dev/reference/prior_summary.md)
  shows the result: try `prior_summary(demo_fit, verbose = TRUE)`, where
  `rule` is the default `prior` with these filled in, and `prior` is the
  resolved prior.

- response:

  Named list describing the response:

  - `auxiliary`: Named list of response auxiliaries. Built-in families
    use `trials()` and
    [`weights()`](https://rdrr.io/r/stats/weights.html); custom families
    can define additional names, used like `y | name(column)`. For
    example,
    `list(trials = list(required = TRUE, operations = c("epred", "log_lik", "simulate")))`.
    `operations` lists where the column is needed; allowed values are
    `"epred"`, `"log_lik"`, `"simulate"`, and `"garma"`.
    [`weights()`](https://rdrr.io/r/stats/weights.html) is always
    available. Use it in `jags$likelihood` to support weights.

  - `validate(y, data, data_columns)`: Function that raises an error for
    invalid responses.

  - `is_discrete`: Logical. Integer-valued responses? Defaults to
    `FALSE`.

  - `has_rate`, `observed(y, data, rate)`, `probability(rate)`: For
    responses per trial, as in
    [`binomial()`](https://rdrr.io/r/stats/family.html).

- garma:

  AR/MA support for this family. `NULL` (unsupported) or a named list
  with `observed_r(y, data, threshold)` and `observed_jags(context)`:
  Matching response transformations in R and JAGS for the GARMA
  recurrence, before applying the link function, using `threshold` where
  needed. Optionally `validate_threshold(threshold, data)` which stop()s
  if invalid.

- glm_family:

  `NULL` or a stats family, e.g.,
  [`poisson()`](https://rdrr.io/r/stats/family.html), that approximates
  the model of `mu`.
  [`stats::glm.fit()`](https://rdrr.io/r/stats/glm.html) uses it to find
  initial values for change points. If `NULL`, they are equally spaced.

- floor:

  Named numeric vector. Numerical floor on the dpar for each
  observation, e.g., `c(sigma = 1e-5)`.

## Value

An `mcpfamily` object.

## Details

Define a response distribution for
[`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md), e.g.,
one that `mcp` does not ship. The built-in families, such as
[`gaussian()`](https://rdrr.io/r/stats/family.html) and
[`negbinomial()`](https://lindeloev.github.io/mcp/dev/reference/negbinomial.md),
are defined using `mcpfamily()` too. See [Custom
families](https://lindeloev.github.io/mcp/articles/custom_families.html)
for a worked example.

Only the elements needed to fit the model are required: `family`,
`links`, and `jags$likelihood`, with priors supplied as defaults or in
[`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md). Methods
that are needed later (e.g., `r$simulate` for
[`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md))
raise an error on use.

## Examples

``` r
# A Student-t family with mu, sigma, and nu
student = mcpfamily(
  family = "student",
  links = c(mu = "identity", sigma = "log", nu = "log"),
  lower = c(sigma = 0, nu = 1),
  jags = list(likelihood = function(context) {
    paste0(context$y, " ~ dt(", context$dpar("mu"), ", 1 / ", context$dpar("sigma"), "^2, ", context$dpar("nu"), ")")
  }),
  r = list(
    log_lik = function(y, dpars, data) stats::dt((y - dpars$mu) / dpars$sigma, dpars$nu, log = TRUE) - log(dpars$sigma),
    simulate = function(n, dpars, data, rate = FALSE) dpars$mu + dpars$sigma * stats::rt(n, dpars$nu)
  )
)

# Priors are required for parameters without a default prior
data = data.frame(x = 1:20, y = c(rnorm(10), rnorm(10, 3)))
prior = list(Intercept_1 = "dnorm(0, 10)", Intercept_2 = "dnorm(0, 10)", sigma_1 = "dnorm(0, 5) T(0, )", nu_1 = "dgamma(2, 0.1) T(1, )")
fit = mcp(list(y ~ 1, ~ 1), data, family = student, prior = prior, par_x = "x", sample = FALSE)
fit$family
#> Family: student
#> Links: mu = identity; sigma = identity; nu = identity
```
