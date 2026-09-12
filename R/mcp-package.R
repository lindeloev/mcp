#' mcp: Multiple Change Point Regression in R
#'
#' @description
#' Flexible and informed regression with Multiple Change Points. `mcp` can infer change points in models on means,
#' variances (and other distributional parameters), autocorrelation structure, and any combination of these, as well as the parameters of the
#' segments in between. All parameters are estimated with uncertainty, and posterior predictive intervals are
#' supported - also near the change points. `mcp` supports hypothesis testing via Savage-Dickey density
#' ratios, posterior contrasts, and PSIS-LOO/WAIC model comparison.
#'
#' @details
#' **The mcp model**
#'
#' An `mcp` model divides a continuous predictor \eqn{x} into \eqn{K} segments separated by
#' ordered change points \eqn{\Delta_1 < \dots < \Delta_{K-1}}. In each segment \eqn{k \in \{1, \dots, K\}},
#' the linear predictor \eqn{\eta_i} is evaluated directly from the segment-local distance \eqn{(x_i - \Delta_{k-1})}:
#'
#' \deqn{\eta_i = \alpha_k + \beta_{k,1} (x_i - \Delta_{k-1}) \quad (\text{with } \Delta_0 = 0)}
#'
#' where the segment-start level \eqn{\alpha_k} is freely estimated for the first and disjoined segments, as in non-segmented regression, and inherited continuously for joined segments:
#'
#' \deqn{\alpha_k = \begin{cases}
#'   \beta_{k,0}, & \text{Disjoined segments } (\sim \texttt{1 + x}, \text{ including } k = 1) \\
#'   \alpha_{k-1} + \beta_{k-1,1} (\Delta_{k-1} - \Delta_{k-2}), & \text{Joined segments } (k \ge 2, \sim \texttt{0 + x})
#' \end{cases}}
#'
#' In all segments, estimated slope and intercept parameters are absolute values (not changes relative to the preceding segment).
#'
#' If additional continuous covariates or categorical factors are included (e.g., `+ z + group`),
#' they enter additively on their original scale (\eqn{\dots + \sum \gamma_{k,j} z_{j,i}}); only the
#' change-point predictor \eqn{x} is converted to segment-local coordinates.
#'
#' Distributional parameters (\code{sigma()}, \code{shape()}, etc.) and autoregressive terms (\code{ar()}, \code{ma()})
#' follow this exact same segmented structure on their respective link scales.
#'
#' **Extended formulas and model features**
#'
#' * *Distributional regression:* Model residual variance and other distributional parameters across segments,
#'   e.g., `~ sigma(1 + x)` or `~ shape(1)`.
#' * *Time-series residuals:* Model serial dependence using `ar(p)` and `ma(q)` terms with a generalized link-scale
#'   recurrence that spans continuously across change points.
#' * *Group-level effects:* Add hierarchical (random) intercepts, slopes, and change points, e.g., `~ 1 + (1|id)`
#'   or `1 + (1|id) ~ 0 + x`.
#' * *Model comparison and hypothesis testing:* Compare models using PSIS-LOO (`loo()`) or WAIC (`waic()`), and test directional
#'   or point-null hypotheses (`hypothesis()`) using Savage-Dickey density ratios.
#'
#' See [the mcp website](https://lindeloev.github.io/mcp/) for worked examples and vignettes.
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
