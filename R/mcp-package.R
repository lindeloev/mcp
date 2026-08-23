#' mcp: Multiple Change Point Regression in R
#'
#' @description
#' Flexible and informed regression with Multiple Change Points. `mcp` can infer change points in models on means,
#' variances (and other distributional parameters), autocorrelation structure, and any combination of these, as well as the parameters of the
#' segments in between. All parameters are estimated with uncertainty, and prediction intervals are
#' supported - also near the change points. `mcp` supports hypothesis testing via Savage-Dickey density
#' ratios, posterior contrasts, and PSIS-LOO/WAIC model comparison.
#'
#' @details
#' **The mcp model**
#'
#' For a continuous predictor \eqn{x}, segment 1 (\eqn{x \le \Delta_1}) is a standard regression
#' model with the intercept defined at \eqn{x = 0}:
#'
#' \deqn{\eta_{1, i} = f_1(x_i, \mathbf{\beta}_1) = \beta_{1,0} + \beta_{1,1} x_i + \dots}
#'
#' For subsequent segments \eqn{k \in \{2, \dots, K\}}, each change point \eqn{\Delta_{k-1}} marks the
#' onset of a new segment. The model defines a segment-local predictor \eqn{X_{k,i}} measured from the
#' segment onset \eqn{\Delta_{k-1}} and capped at the segment end \eqn{\Delta_k}:
#'
#' \deqn{X_{k,i} = \min(x_i, \Delta_k) - \Delta_{k-1}}
#'
#' with \eqn{\Delta_K = \max(\mathbf{x})}. The local predictor \eqn{X_{k,i}} is zero at segment onset
#' (\eqn{x_i \le \Delta_{k-1}}) and plateaus at \eqn{\Delta_k - \Delta_{k-1}} for \eqn{x_i \ge \Delta_k}.
#'
#' **Joined (continuous) segments (`~ 0 + x`)**
#'
#' In joined models, segment \eqn{k \ge 2} has no intercept and specifies a new segment slope \eqn{\beta_{k,1}}:
#'
#' \deqn{\eta_i = \eta_{1,i}(X_{1,i}) + \sum_{k=2}^K [x_i \ge \Delta_{k-1}] \cdot \beta_{k,1} X_{k,i}}
#'
#' where \eqn{X_{1,i} = \min(x_i, \Delta_1)}. Because earlier segments plateau at their boundary values,
#' the expectation is automatically continuous across change points without requiring parameter constraints.
#' In both joined and disjoined segments, slope and intercept parameters are absolute values (not differences
#' relative to the previous segment).
#'
#' **Disjoined (discontinuous) segments (`~ 1 + x`)**
#'
#' When segment \eqn{k} contains an explicit intercept (\eqn{\beta_{k,0}}), it introduces a step change at
#' \eqn{x_i = \Delta_{k-1}} and previous segments are truncated at \eqn{[x_i < \Delta_{k-1}]}.
#'
#' **Extended formulas and model features**
#'
#' * *Distributional regression:* Model residual variance and other distributional parameters across segments,
#'   e.g., `~ sigma(1 + x)` or `~ shape(1)`.
#' * *Time-series residuals:* Model serial dependence using `ar(p)` and `ma(q)` terms with a generalized link-scale
#'   recurrence that spans continuously across change points.
#' * *Group-level effects:* Add hierarchical (random) intercepts, slopes, and change points, e.g., `~ 1 + (1|id)`
#'   or `1 + (1|id) ~ 0 + x`.
#' * *Model comparison and hypothesis testing:* Compare models using LOO-CV (`loo()`, `waic()`) and test directional
#'   or point-null hypotheses (`hypothesis()`) using Savage-Dickey density ratios.
#'
#' See [the mcp website](https://lindeloev.github.io/mcp/) for worked examples and vignettes.
#'
#' @references
#' * Lindeløv, J. K. (2020). mcp: An R Package for Regression With Multiple Change Points. *OSF Preprints*. [doi:10.31219/osf.io/fzqxv](https://doi.org/10.31219/osf.io/fzqxv)
#' * Carlin, B. P., Gelfand, A. E., & Smith, A. F. (1992). Hierarchical Bayesian Analysis of Changepoint Problems. *Applied Statistics*, 41(2), 389–405. [doi:10.2307/2347570](https://doi.org/10.2307/2347570)
#' * Stephens, D. A. (1994). Bayesian retrospective multiple-changepoint identification. *Applied Statistics*, 43(1), 159–178. [doi:10.2307/2986119](https://doi.org/10.2307/2986119)
#'
#' @seealso \code{\link{mcp}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
"_PACKAGE"
