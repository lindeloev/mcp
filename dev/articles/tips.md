# Diagnosing and fixing problems

## Convergence

A common problem when using MCMC is lacking convergence between chains.
This will show up as large `rhat` values, low `ess_bulk` or `ess_tail`,
and non-converging lines in `plot_pars(fit)`. By default,
[`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) and
[`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
warn when `rhat > 1.01` or ESS \< 400.

- The first thing to try is always to make the model warm up longer and
  sample longer, e.g.: `mcp(model, data, warmup = 10000, iter = 20000)`.
  This is often needed with many (e.g., 5+) change points.

- Change points can have several modes, e.g., at two candidate jumps in
  the data. Chains may stay in one mode for thousands of iterations, so
  sampling longer may not help. Look for modes using
  `plot_pars(fit, regex_pars = "^cp_")`. If there are several, the
  posterior is genuinely uncertain between them, and more chains (e.g.,
  `mcp(..., chains = 6)`) give a better sense of how probable each one
  is.

- Center predictors other than the change-point variable, e.g.,
  `data$z = data$z - mean(data$z)`. When `z` is far from zero, its
  coefficient and the intercept are strongly correlated, which slows
  sampling of both and, in turn, the change points. After centering, the
  intercept is the value at the mean of `z`.

- It can be a sign of a deeper non-identifiability in the model. This
  will show up as strong correlations in the joint distribution of any
  pair of implicated parameters:
  `plot_pars(fit, pars = c("Intercept_1", "Intercept_2"), type = "hex")`.
  This may give you ideas how to change the model.

- You can set the initial values for the JAGS sampler using, e.g.,
  `mcp(..., inits = list(cp_1 = 30, sigma_1 = 2))`. Names must be
  parameters in `mcp_pars(fit)`. Parameters with constant priors and
  parameters reparameterized for faster sampling cannot be initialized.
  Initial values for change points must be ordered.

  - By default, all chains start at the same change points. Change
    points cannot pass each other, so a chain that starts with change
    points in the wrong places may get stuck there. If you give chains
    different change point inits, disagreement between chains may
    reflect poor starting points rather than a problem with the model.

## Speed

A lot of data and complicated models will slow down fitting.

- Run the chains in parallel by selecting a future backend before
  fitting, e.g., `future::plan(future::multisession, workers = 4)`. The
  worker pool is reused across calls. Shut it down afterward with
  `future::plan(future::sequential)`. Without a parallel future plan,
  chains run sequentially.

- While more data usually slows warmup and sampling, it also typically
  has better identifiability and faster convergence. You can therefore
  try reducing both, e.g., `mcp(..., warmup = 300, iter = 1000)`. For
  JAGS, warmup iterations are used for sampler adaptation and discarded
  before sampling.

## Errors or won’t run

Most of these problems should stem from inappropriate priors and such
problems may be exacerbated by fragile link functions (e.g.,
`binomial(link = "identity")`). The article on [priors in
mcp](https://lindeloev.github.io/mcp/dev/articles/priors.md) may be
helpful, but in particular:

- Errors on “directed cycle” usually stem from using parameters in
  priors. For example, if you set
  `prior = list(Intercept_1 = "dnorm(Intercept_2, 1)", Intercept_2 = "dnorm(Intercept_1, 1)")`
  this is an infinite regress.

- Errors on “incompatible with parent nodes” usually stem from
  impossible values. For example, setting `prior = list(sigma_1 = -1)`
  when the sigma-link is identity or setting
  `prior = list(cp_1 = 10, cp_2 = 5)`, which are not ordered.

If you encounter problems not readily resolved, don’t hesitate to [raise
a Github Issue](https://github.com/lindeloev/mcp/issues), asking for
help or filing a bug report.
