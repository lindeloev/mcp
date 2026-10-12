#' mcp: Regression with Multiple Change Points
#'
#' @description
#' `mcp` is regression with multiple change points. At its simplest, it feels like `lm()`, but with one formula per
#' segment. At its most flexible, it aims to be **[`brms`](https://paulbuerkner.com/brms/) for change points**: GLM
#' families (Gaussian, binomial, Bernoulli, Poisson, negative binomial), group-level effects (random effects), AR/MA,
#' and regression on distributional parameters (e.g., `sigma`, `shape`), all inferred with full Bayesian uncertainty,
#' including the change points themselves.
#'
#' Posterior predictive intervals are supported - also near the change points. `mcp` supports hypothesis testing via Savage-Dickey density
#' ratios, posterior contrasts, and PSIS-LOO/WAIC model comparison.
#'
#' See [mcp()] to fit models, [mcp-formula] for the model and formula syntax, and [mcp-priors] for priors.
#'
#' @section Extended formulas and model features:
#' * *Distributional regression:* Model residual variance and other distributional parameters across segments,
#'   e.g., `~ sigma(1 + x)` or `~ shape(1)`.
#' * *Time-series residuals:* Model serial dependence using `ar(p)` and `ma(q)` terms with a generalized link-scale
#'   recurrence.
#' * *Group-level effects:* Add hierarchical (random) intercepts, slopes, and change points, e.g., `~ 1 + (1|id)`
#'   or `~ cp((1|id)) + 0 + x`.
#' * *Model comparison and hypothesis testing:* Compare models using PSIS-LOO (`loo()`) or WAIC (`waic()`), and test hypotheses
#'   (`hypothesis()`) using Savage-Dickey density ratios for point-null tests or posterior-to-prior odds for directional tests.
#'
#' See [the mcp website](https://lindeloev.github.io/mcp/) for worked examples.
#'
#' @section Options:
#' * `mcp_quiet`: Logical. The default of `mcp(..., quiet)`. `options(mcp_quiet = TRUE)` suppresses routine JAGS output
#'   and mcp sampling-status messages for the session.
#' * `mcp_reparam`: Logical. `TRUE` (default) lets [mcp()] sample some default priors on scales that mix
#'   faster. This changes only `fit$jags_code` and sampling efficiency, not the model, the priors, or the posterior.
#'   `options(mcp_reparam = FALSE)` writes the plain model, which is mostly useful for seeing how `mcp` works.
#'   Ignored when `jags_code` is supplied. Read more in
#'   [From formulas to sampler code](https://lindeloev.github.io/mcp/articles/implementation.html#reparameterizations).
#'
#' @inherit mcp references
#'
#' @seealso \code{\link{mcp}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
"_PACKAGE"

