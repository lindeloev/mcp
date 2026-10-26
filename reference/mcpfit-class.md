# Class `mcpfit` of Models Fitted with the mcp Package

Models fitted with the
[`mcp`](https://lindeloev.github.io/mcp/reference/mcp.md) function are
represented as an `mcpfit` object which contains the user input (model,
data, family), derived model characteristics (prior, parameter names,
and jags code), and the fit (prior and/or posterior MCMC draws).

## Details

See `methods(class = "mcpfit")` for an overview of available methods.

Components:

- `call`: The matched call to
  [`mcp()`](https://lindeloev.github.io/mcp/reference/mcp.md).

- `model`: A list of user-provided formulas.

- `data`: The user-provided data frame reduced to model-used columns.

- `family`: An `mcpfamily` object.

- `prior`: A named list of priors.

- `mcmc_post` and `mcmc_prior`:
  [`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) objects
  with posterior and prior draws, respectively. Do not access these
  directly; use as_draws(fit) or similar.

- `jags_code`: A string with JAGS code; use `cat(fit$jags_code)` to show
  it. [Read
  more](https://lindeloev.github.io/mcp/articles/implementation.html).

- `simulate`: A function to simulate data from supplied parameter
  values.

- `version`: A named list of
  [package_version](https://rdrr.io/r/base/numeric_version.html)s for
  `mcp`, `R`, `rjags`, and `jags`, e.g., `fit$version$mcp >= "0.4.0"`.
  `rjags` and `jags` are `NULL` for fits with `sample = "none"`. `NULL`
  for fits made with mcp \< 0.4.0.

- `.internal`: Information used internally by mcp.
