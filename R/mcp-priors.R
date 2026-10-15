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
#' @section Regression priors:
#' 
#' Weakly informative defaults are used where priors are not specified. These are the conventional defaults of
#' `rstanarm` for intercepts and coefficients and of `brms` for scale parameters, with a few documented deviations.
#' They are designed for stable estimation and prediction, but should be justified before hypothesis testing.
#' Use `prior_summary(fit, verbose = TRUE)` to see the default priors of a model and the rules behind them.
#' For more details, read [the article on mcp default priors](https://lindeloev.github.io/mcp/articles/default_priors.html) for all families, links,
#' and distributional parameters, and for how the defaults deviate from `rstanarm` and `brms`.
#'
#' @section Change point priors: 
#' 
#' Change points are always ordered: each change point keeps its specified prior, restricted to ordered change points.
#' The default is a flat Dirichlet prior on the segment lengths, written as `dunif(min(x), max(x))` on each change point
#' ([read more](https://lindeloev.github.io/mcp/articles/default_priors.html#change-points)).
#' Set a minimum segment length using the previous change point in the lower limit, e.g., `cp_2 = "dunif(cp_1 + 20, max(x))"`.
#' [Read more](https://lindeloev.github.io/mcp/articles/priors.html#cp_prior) about setting priors on change points.
#'
#' @section Notes on priors:
#'
#'   * *Change point priors* without `dunif()` or `T()` are truncated to the observed range of x,
#'       e.g., `dnorm(30, 10)` becomes `dnorm(30, 10) T(min(x), max(x))`.
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
