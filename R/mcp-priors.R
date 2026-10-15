#' Priors in mcp
#'
#' @description
#' The `prior` argument of [mcp()] is a named list. Names are parameter names (`cp_1`, `Intercept_1`, `x_2`,
#' `sigma_1`, etc.; see [mcp-formula] for how they arise) and the values specify their priors.
#' Use [prior_summary()] to see the resolved priors of a fitted model, including defaults.
#' See [Working with priors](https://lindeloev.github.io/mcp/articles/priors.html) for worked examples
#' and prior predictive checks.
#'
#' @details
#' Each value in the list is either
#'
#'  * A distribution string (e.g., `Intercept_1 = "dnorm(0, 1) T(0, )"`). Distributions use JAGS names
#'      and argument order, including `T()` truncation, but with conventional scales instead of JAGS
#'      precision: SD for `dnorm()`, SD on the log scale for `dlnorm()` (`sdlog` in [stats::dlnorm()]),
#'      and scale for `dt()`, `ddexp()`, and `dlogis()`. Change points also accept `dirichlet(alpha)`,
#'      and all priors accept data-dependent terms such as `sd(y)` (see below).
#'  * A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed value.
#'  * A model parameter name (e.g., `Intercept_2 = "Intercept_1"`), sharing that parameter's value.
#'      Note that sharing coefficients via `same()` in the segment formulas (e.g., `~ same(1)`)
#'      is generally preferred. If two group-level deviations are shared via the prior,
#'      they will need to have the same grouping variable.
#'
#' Weakly informative defaults are used where priors are not specified. These are designed for stable
#' estimation and prediction, but should be justified before hypothesis testing.
#' Default priors follow `rstanarm` for intercepts and coefficients and `brms` for scale parameters, with
#' a few documented deviations. For example, with an identity link:
#'
#'  * *Coefficients* are autoscaled as in `rstanarm`, e.g., `normal(0, 2.5 * sd(y) / sd(x))` for a slope on `x`.
#'      This makes the prior largely invariant to the units of `x` and `y`.
#'  * *Intercepts* are `normal(mean(y), 2.5 * sd(y))`. Like `rstanarm` and `brms`, the default prior on `Intercept_1`
#'      applies at the start of the data, `min(x)`. The reported `Intercept_1` is still at `x = 0`, as in [lm()], so
#'      `fit$prior` shows its prior given the x-terms, e.g., `dnorm(mean(y) - x_1 * min(x), 2.5 * sd(y))`.
#'      Other predictors are not centered. User-specified priors on `Intercept_1` apply directly at `x = 0`.
#'  * *`sigma` and group-level SDs* follow `brms`, e.g., `student_t(3, 0, max(2.5, mad(y))) T(0, )`.
#'
#' Use `prior_summary(fit, verbose = TRUE)` to see the default priors of a model and the rules behind them.
#' See [Default priors](https://lindeloev.github.io/mcp/articles/default_priors.html) for all families, links,
#' and distributional parameters, and for how the defaults deviate from `rstanarm` and `brms`.
#'
#' The default prior on change points is `dirichlet(1)` (uniform order statistics).
#' For a single change point, this is the Beta(1, 1) / Uniform distribution over `[min(x), max(x)]`.
#' For multiple change points, it corresponds to a flat Dirichlet distribution over segment lengths.
#' You can also explicitly set `cp_i = "dirichlet(alpha)"` with the same positive `alpha` for all
#' change points to regularize spacing (`alpha > 1` penalizes change points from occurring close together,
#' while `alpha < 1` favors clustering). Under the hood, `dirichlet(1)` is sampled as sorted uniforms
#' and other `alpha` as an exact sequential stick-breaking Beta chain for fast and robust sampling.
#' [Read more](https://lindeloev.github.io/mcp/articles/priors.html#cp_prior).
#'
#' @section Notes on priors:
#'
#'   * *Ordered change point priors:* Default population-level `cp_i` priors are ordered and the ordering is imposed through the priors. For user-defined priors,
#'       `mcp` adds truncation (e.g., `T(cp_1, )`) only when the prior has neither
#'       explicit truncation nor an inherently bounded form such as `dunif()` or
#'       `dirichlet()`.
#'   * *Data-dependent terms:* If `mcp` encounters a data-dependent term like `min(time)`, `max(time)`, `mean(response)`, `median(response)`, `sd(response)`, or `mad(response)` in the prior string, they are resolved from the model data so a numerical value is passed to JAGS. The following terms are also allowed: `n_segments()` and `n_cp()`.
#'       The older constants `MINX`, `MAXX`, `MEANX`, `SDX`, `MINY`, `MAXY`, `MEANY`, `SDY`, and `N_CP` remain accepted with a deprecation warning.
#'   * *Group-level change points:* Group-specific locations follow a hierarchical
#'       normal distribution around their population change point, truncated so
#'       that realized locations remain in the observed range and ordered.
#'       [Read more](https://lindeloev.github.io/mcp/articles/default_priors.html#group-level-change-points).
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
