# mcp: Regression with Multiple Change Points

`mcp` is regression with multiple change points. At its simplest, it
feels like [`lm()`](https://rdrr.io/r/stats/lm.html), but with one
formula per segment. At its most flexible, it aims to be
**[`brms`](https://paulbuerkner.com/brms/) for change points**: GLM
families (Gaussian, binomial, Bernoulli, Poisson, negative binomial),
group-level effects (random effects), AR/MA, and regression on
distributional parameters (e.g., `sigma`, `shape`), all inferred with
full Bayesian uncertainty, including the change points themselves.

Posterior predictive intervals are supported - also near the change
points. `mcp` supports hypothesis testing via Savage-Dickey density
ratios, posterior contrasts, and PSIS-LOO/WAIC model comparison.

See [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) to
fit models,
[mcp-formula](https://lindeloev.github.io/mcp/dev/reference/mcp-formula.md)
for the model and formula syntax, and
[mcp-priors](https://lindeloev.github.io/mcp/dev/reference/mcp-priors.md)
for priors.

## Extended formulas and model features

- *Distributional regression:* Model residual variance and other
  distributional parameters across segments, e.g., `~ sigma(1 + x)` or
  `~ shape(1)`.

- *Time-series residuals:* Model serial dependence using `ar(p)` and
  `ma(q)` terms with a generalized link-scale recurrence that spans
  continuously across change points.

- *Group-level effects:* Add hierarchical (random) intercepts, slopes,
  and change points, e.g., `~ 1 + (1|id)` or `1 + (1|id) ~ 0 + x`.

- *Model comparison and hypothesis testing:* Compare models using
  PSIS-LOO
  ([`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md))
  or WAIC
  ([`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)),
  and test hypotheses
  ([`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md))
  using Savage-Dickey density ratios for point-null tests or
  posterior-to-prior odds for directional tests.

See [the mcp website](https://lindeloev.github.io/mcp/) for worked
examples.

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
  Introduced the Gibbs sampling approach for continuous change points in
  segmented regression and hierarchical models, providing the
  computational foundation used by BUGS and JAGS.

## See also

[`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md)

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
