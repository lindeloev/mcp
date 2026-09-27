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
#'   recurrence that spans continuously across change points.
#' * *Group-level effects:* Add hierarchical (random) intercepts, slopes, and change points, e.g., `~ 1 + (1|id)`
#'   or `1 + (1|id) ~ 0 + x`.
#' * *Model comparison and hypothesis testing:* Compare models using PSIS-LOO (`loo()`) or WAIC (`waic()`), and test hypotheses
#'   (`hypothesis()`) using Savage-Dickey density ratios for point-null tests or posterior-to-prior odds for directional tests.
#'
#' See [the mcp website](https://lindeloev.github.io/mcp/) for worked examples.
#'
#' @inherit mcp references
#'
#' @seealso \code{\link{mcp}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
"_PACKAGE"


.onLoad = function(libname, pkgname) {
  if (requireNamespace("posterior", quietly = TRUE)) {
    registerS3method("as_draws", "mcpfit", as_draws.mcpfit, envir = asNamespace("posterior"))
    registerS3method("as_draws_df", "mcpfit", as_draws_df.mcpfit, envir = asNamespace("posterior"))
    registerS3method("as_draws_array", "mcpfit", as_draws_array.mcpfit, envir = asNamespace("posterior"))
    registerS3method("as_draws_matrix", "mcpfit", as_draws_matrix.mcpfit, envir = asNamespace("posterior"))
    registerS3method("as_draws_rvars", "mcpfit", as_draws_rvars.mcpfit, envir = asNamespace("posterior"))
    registerS3method("nchains", "mcpfit", nchains.mcpfit, envir = asNamespace("posterior"))
    registerS3method("ndraws", "mcpfit", ndraws.mcpfit, envir = asNamespace("posterior"))
    registerS3method("niterations", "mcpfit", niterations.mcpfit, envir = asNamespace("posterior"))
  }
  if (requireNamespace("coda", quietly = TRUE)) {
    registerS3method("as.mcmc", "mcpfit", as.mcmc.mcpfit, envir = asNamespace("coda"))
  }
  if (requireNamespace("rstantools", quietly = TRUE)) {
    registerS3method("posterior_epred", "mcpfit", posterior_epred.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("posterior_predict", "mcpfit", posterior_predict.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("posterior_linpred", "mcpfit", posterior_linpred.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("log_lik", "mcpfit", log_lik.mcpfit, envir = asNamespace("rstantools"))
  }
  setHook(packageEvent("rstantools", "onLoad"), function(...) {
    registerS3method("posterior_epred", "mcpfit", posterior_epred.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("posterior_predict", "mcpfit", posterior_predict.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("posterior_linpred", "mcpfit", posterior_linpred.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("log_lik", "mcpfit", log_lik.mcpfit, envir = asNamespace("rstantools"))
  })
}
