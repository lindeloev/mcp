# Class `mcpfit` of Models Fitted with the mcp Package

Models fitted with the
[`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) function
are represented as an `mcpfit` object which contains the user input
(model, data, family), derived model characteristics (prior, parameter
names, and jags code), and the fit (prior and/or posterior mcmc
samples).

## Details

See `methods(class = "mcpfit")` for an overview of available methods.

User-provided information (see
[`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) for more
details):

## Slots

- `model`:

  A list of formulas, making up the model. Provided by user. See
  [`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) for more
  details.

- `data`:

  A data frame. Provided by user. See
  [`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) for more
  details.

- `family`:

  An `mcpfamily` object. Provided by user. See
  [`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) for more
  details.

- `prior`:

  A named list. Provided by user. See
  [`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) for more
  details.

- `mcmc_post`:

  An [`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) object
  with posterior samples.

- `mcmc_prior`:

  An [`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) object
  with prior samples.

- `loglik`:

  An (Nchains \* Ndraws) by N-observed-responses matrix of
  log-likelihoods.

- `pars`:

  A list of character vectors of model parameter names.

- `jags_code`:

  A string with jags code. Use `cat(fit$jags_code)` to show it.

- `simulate`:

  A method to simulate and predict data.

- `.internal`:

  Information that is used internally by mcp.
