#' mcp: Multiple Change Point Regression in R
#'
#' @description
#' Flexible and informed regression with Multiple Change Points. `mcp` can infer change points in models on means,
#' variances (and other distributional parameters), autocorrelation structure, and any combination of these, as well as the parameters of the
#' segments in between. All parameters are estimated with uncertainty, and posterior predictive intervals are
#' supported - also near the change points. `mcp` supports hypothesis testing via Savage-Dickey density
#' ratios, posterior contrasts, and PSIS-LOO/WAIC model comparison.
#'
#' @inheritSection mcp The mcp model
#' @inheritSection mcp Time-series residuals (link-scale observation-driven GARMA)
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
