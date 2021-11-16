#' Regression with Multiple Change Points
#'
#' @docType package
#' @name mcp-package
#' @aliases mcp-package
#'
#' @description
#' The \pkg{mcp} package provides an interface to fit regression models with
#' multiple change points between generalized linear segments, optionally with
#' per-segment variance and autocorrelation structures.
#'
#' The main function of \pkg{mcp} is the `mcp()` function, which uses a formula
#' syntax to specify a wide range of change point models. Based on the supplied
#' data, formulas, and additional information, it writes JAGS code on the fly
#' and use \pkg{rstan} to fit the model, optionally in parallel to speed up
#' sampling. You will need to install JAGS for `mcp()` to work.
#'
#' A large number of post-processing methods can be applied. These include
#'
#' * Summarise fits using `summary()`, `fixef()`, and `ranef()`.
#' * Visualize fits using `plot()` and individual parameters using `plot_pars()`.
#' * Test hypotheses using `hypothesis()` and `loo()`.
#'
#' Extensive documentation with worked examples is available at the
#' [mcp website](https://lindeloev.github.io/mcp).
NULL
