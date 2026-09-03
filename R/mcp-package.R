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
#' An `mcp` model divides a continuous predictor \eqn{x} into \eqn{K} segments separated by
#' ordered change points \eqn{\Delta_1 < \dots < \Delta_{K-1}}. In each segment \eqn{k \in \{1, \dots, K\}},
#' the linear predictor \eqn{\eta_i} is evaluated directly from the segment-local distance \eqn{(x_i - \Delta_{k-1})}:
#'
#' \deqn{\eta_i = \begin{cases}
#'   \beta_{1,0} + \beta_{1,1} x_i, & \text{First segment } (k = 1) \\
#'   \beta_{k,0} + \beta_{k,1} (x_i - \Delta_{k-1}), & \text{Disjoined segments } (k \ge 2, \sim \texttt{1 + x}) \\
#'   \alpha_k + \beta_{k,1} (x_i - \Delta_{k-1}), & \text{Joined segments } (k \ge 2, \sim \texttt{0 + x})
#' \end{cases}}
#'
#' where \eqn{\alpha_k} is the baseline level inherited continuously from the preceding segment:
#'
#' \deqn{\alpha_k = \alpha_{k-1} + \beta_{k-1,1} (\Delta_{k-1} - \Delta_{k-2}) \quad (\text{with } \Delta_0 = 0)}
#'
#' In all segments, slope and intercept parameters are absolute values (not changes relative to the preceding segment).
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
