#' Priors in mcp
#'
#' @description
#' The `prior` argument of [mcp()] is a named list. Names are parameter names (`cp_1`, `Intercept_1`, `x_2`,
#' `sigma_1`, etc.; see [mcp-formula] for how they arise) and the values specify their priors.
#' Use [prior_summary()] to see the resolved priors of a fitted model, including defaults.
#'
#' @details
#' Each value in the list is either
#'
#'  * A distribution in `mcp`'s JAGS-string syntax (e.g.,
#'      `Intercept_1 = "dnorm(0, 1) T(0,)"`) indicating a
#'      conventional prior distribution. Data-calibrated, regularizing defaults
#'      are used where priors are not specified. These are designed for stable
#'      estimation and prediction, but should be justified before hypothesis testing.
#'      `mcp` uses conventional distribution scales rather than JAGS precision:
#'      SD for `dnorm()`, scale for `dt()`, `ddexp()`, and `dlogis()`, and
#'      log-SD for `dlnorm()`.
#'  * A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed value.
#'  * A model parameter name (e.g., `Intercept_2 = "Intercept_1"`), equating parameters via JAGS.
#'      Note that sharing coefficients via `same()` in the segment formulas (e.g., `~ same(1)`)
#'      is generally preferred. If two group-level deviations are shared via the prior,
#'      they will need to have the same grouping variable.
#'
#' The default prior on change points is `dirichlet(1)` (uniform order statistics).
#' For a single change point, this is the Beta(1, 1) / Uniform distribution over `[min(x), max(x)]`.
#' For multiple change points, it corresponds to a flat Dirichlet distribution over segment lengths.
#' You can also explicitly set `cp_i = "dirichlet(alpha)"` with the same positive `alpha` for all
#' change points to regularize spacing (`alpha > 1` penalizes change points from occurring close together,
#' while `alpha < 1` favors clustering). Under the hood, this is parameterized as an exact sequential
#' stick-breaking Beta chain for fast and robust sampling.
#' [Read more](https://lindeloev.github.io/mcp/articles/priors.html).
#'
#' @section Notes on priors:
#'
#'   * *Ordered change point priors:* Default population-level `cp_i` priors are ordered and the ordering is imposed through the priors. For user-defined priors,
#'       `mcp` adds truncation (e.g., `T(cp_1, )`) only when the prior has neither
#'       explicit truncation nor an inherently bounded form such as `dunif()` or
#'       `dirichlet()`.
#'   * *Data-dependent terms:* If `mcp` encounters a data-dependent term like `min(time)`, `max(time)`, `median(response)`, or `mad(response)` in the prior string, they are resolved from the model data so a numerical value is passed to JAGS. The following terms are also allowed: `n_segments()` and `n_cp()`.
#'       The older constants `MINX`, `MAXX`, `MEANX`, `SDX`, `MINY`, `MAXY`, `MEANY`, `SDY`, and `N_CP` remain accepted with a deprecation warning.
#'   * *Group-level change points:* Group-specific locations follow a hierarchical
#'       normal distribution around their population change point, truncated so
#'       that realized locations remain in the observed range and ordered.
#'   * *Parameterization:* Prior strings use conventional scale parameterizations. `mcp` converts
#'       these to the parameterization required by JAGS when generating code:
#'       inverse variance for `dnorm()`, `dt()`, and `dlnorm()`, and inverse
#'       scale for `ddexp()` and `dlogis()`. Use
#'       `prior_summary(fit)` for resolved priors and
#'       `prior_summary(fit, verbose = TRUE)` for their rules and descriptions.
#'
#' @seealso [mcp()], [prior_summary()], [mcp-formula]
#' @name mcp-priors
#' @aliases mcp_priors
#' @examples
#' model = list(
#'   response ~ 1,
#'   ~ 0 + time,
#'   ~ 1 + time
#' )
#' prior = list(
#'   cp_1 = "dunif(20, 50)",         # Distribution
#'   time_2 = "dnorm(0, 2) T(0, )",  # Truncated distribution
#'   Intercept_3 = 15,               # Fixed value
#'   time_3 = "time_2"               # Shared with another parameter
#' )
#' fit = mcp(model, data = mcp_example_data("demo"), prior = prior, sample = "none")
#' prior_summary(fit)
NULL
