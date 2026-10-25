# mcp: Regression with Multiple Change Points

`mcp` is regression with multiple change points, inferring full Bayesian
uncertainty for all parameters, including the change points themselves.
At its simplest, it feels like [`lm()`](https://rdrr.io/r/stats/lm.html)
with one formula per segment. At its most flexible, it aims to be
**[`brms`](https://paulbuerkner.com/brms/) for change points**:

- *GLMM and distributional regression:* Gaussian, binomial, Bernoulli,
  Poisson, and Negative Binomial families. Group-level (random) effects,
  AR/MA, and regression on distributional parameters (e.g., `sigma`,
  `shape`).

- *Express expert knowledge:* Specify the number of change points, the
  form of each segment, and restrict coefficients via priors. Compare
  hypotheses using
  [`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
  and
  [`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md).

- *Native feel:*
  [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
  [`plot_pars()`](https://lindeloev.github.io/mcp/dev/reference/plot_pars.md),
  [`pp_check()`](https://lindeloev.github.io/mcp/dev/reference/pp_check.md),
  and
  [`as_draws_df()`](https://lindeloev.github.io/mcp/dev/reference/as_draws.mcpfit.md)
  follow the conventions of [`brms`](https://paulbuerkner.com/brms/),
  [`rstanarm`](https://mc-stan.org/rstanarm/), and
  [`posterior`](https://mc-stan.org/posterior/articles/posterior.html),
  adapted to change point regression. Inspect the sampled JAGS model in
  `fit$jags_code`.

You will find many aliases for change point regression in the
literature, including piecewise linear, discontinuity, segmented, and
broken stick regression.

See [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) to
fit models,
[mcp-formula](https://lindeloev.github.io/mcp/dev/reference/mcp-formula.md)
for the model and formula syntax, and
[mcp-priors](https://lindeloev.github.io/mcp/dev/reference/mcp-priors.md)
for priors.

## Extended formulas and model features - and links to worked examples

- *Distributional regression:* Model residual variance and other
  distributional parameters across segments, e.g., `~ sigma(1 + x)` or
  `~ shape(1)`. [Read
  more](https://lindeloev.github.io/mcp/articles/dpar.html).

- *Time-series residuals:* Model serial dependence using `ar(p)` and
  `ma(q)` terms with a generalized link-scale recurrence. [Read
  more](https://lindeloev.github.io/mcp/articles/arma.html).

- *Group-level effects:* Add hierarchical (random) intercepts, slopes,
  and change points, e.g., `~ 1 + (1|id)` or `~ cp((1|id)) + 0 + x`.
  [Read
  more](https://lindeloev.github.io/mcp/articles/group_effects.html).

- *Families:* Binomial, Bernoulli, Poisson, and negative-binomial
  responses with several link functions. [Read
  more](https://lindeloev.github.io/mcp/articles/families.html).

- *Model comparison and hypothesis testing:* Compare models using
  PSIS-LOO
  ([`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md))
  or WAIC
  ([`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)),
  and test hypotheses
  ([`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md))
  using Savage-Dickey density ratios for point-null tests or
  posterior-to-prior odds for directional tests. [Read
  more](https://lindeloev.github.io/mcp/articles/comparison.html).

- *Fits, predictions, and forecasts:*
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  and
  [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  including forecasts beyond the data. [Read
  more](https://lindeloev.github.io/mcp/articles/predict.html).

- *Missing responses:* Missing values in the response are imputed
  jointly with the parameters. [Read
  more](https://lindeloev.github.io/mcp/articles/missing.html).

## Options

- `mcp_quiet`: Logical. The default of `mcp(..., quiet)`.
  `options(mcp_quiet = TRUE)` suppresses routine JAGS output and mcp
  sampling-status messages for the session.

- `mcp_reparam`: Logical. `TRUE` (default) lets
  [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) sample
  some default priors on scales that mix faster. This changes only
  `fit$jags_code` and sampling efficiency, not the model, the priors, or
  the posterior. `options(mcp_reparam = FALSE)` writes the plain model,
  which is mostly useful for seeing how `mcp` works. Ignored when
  `jags_code` is supplied. Read more in [From formulas to sampler
  code](https://lindeloev.github.io/mcp/articles/implementation.html#reparameterizations).

## References

- Lindeløv, J. K. (2020). mcp: An R Package for Regression With Multiple
  Change Points. *OSF Preprints*.
  [doi:10.31219/osf.io/fzqxv](https://doi.org/10.31219/osf.io/fzqxv)
  Introduces the `mcp` package, formula syntax, default priors, and
  workflow for regression with multiple change points across generalized
  linear and time-series models. Newer-than-2020 versions of the paper
  may be available at that link. Please cite the newest version.

- Carlin, B. P., Gelfand, A. E., & Smith, A. F. (1992). Hierarchical
  Bayesian Analysis of Changepoint Problems. *Applied Statistics*,
  41(2), 389–405. [doi:10.2307/2347570](https://doi.org/10.2307/2347570)
  Introduced Gibbs sampling for hierarchical Bayesian change-point
  models, the computational approach that BUGS and JAGS build on.

## See also

[`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md)

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
