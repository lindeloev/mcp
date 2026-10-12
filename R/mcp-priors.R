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
#'  * A distribution (e.g.,
#'      `Intercept_1 = "dnorm(0, 1) T(0,)"`) indicating a
#'      conventional prior distribution. This is not strictly JAGS, since `mcp` uses conventional distribution scales rather than JAGS precision:
#'      SD for `dnorm()`, scale for `dt()`, `ddexp()`, and `dlogis()`, and
#'      log-SD for `dlnorm()`. Weakly informative defaults
#'      are used where priors are not specified. These are designed for stable
#'      estimation and prediction, but should be justified before hypothesis testing.
#'  * A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed value.
#'  * A model parameter name (e.g., `Intercept_2 = "Intercept_1"`), equating parameters via JAGS.
#'      Note that sharing coefficients via `same()` in the segment formulas (e.g., `~ same(1)`)
#'      is generally preferred. If two group-level deviations are shared via the prior,
#'      they will need to have the same grouping variable.
#'
#' Default priors follow `rstanarm` for intercepts and coefficients and `brms` for scale parameters, with
#' a few documented deviations:
#'
#'  * *Coefficients* are autoscaled as in `rstanarm`. For the mean (`mu`),
#'      \deqn{\beta \sim \mathrm{Normal}\left(0, \frac{2.5 \cdot s_y}{\mathrm{sd}(\mathrm{column})}\right)}{beta ~ Normal(0, 2.5 * s_y / sd(column))}
#'      where \eqn{s_y = \mathrm{sd}(y)}{s_y = sd(y)} for identity links and \eqn{s_y = 1}{s_y = 1} for other links.
#'      \eqn{\mathrm{sd}(\mathrm{column})}{sd(column)} is the SD of the coefficient's model-matrix column, including
#'      the 0/1 dummies of factors. This makes the prior largely invariant to the units of `x` and `y`. Logit and
#'      probit links use 1.5 and 1 instead of 2.5. For x-terms, `sd(x)` is taken over all data, and powers such
#'      as `I(x^2)` use `sd((x - min(x))^2)`.
#'  * *Intercepts* are `normal(mean(y), 2.5 * sd(y))` for identity links, `normal(0, 1.5)` for logit, and
#'      `normal(0, 1)` for probit. The latter two are (approximately) uniform on the probability scale.
#'      For log links (Gaussian, Poisson, and negative binomial), intercepts are
#'      `normal(median(log(y)), max(2.5, mad(log(y))))`, rounded to one decimal and with `y <= 0` replaced by 0.1 so that the
#'      log is defined. This robustly centers the prior on the typical log-response, as in `brms` (which uses a
#'      Student-t shape and replaces zeros the same way), whereas `rstanarm` uses `normal(0, 2.5)`.
#'      For a factor without an intercept, e.g., `y ~ 0 + state`, each factor level's coefficient is an intercept and
#'      gets this prior where the segment starts at zero (segment 1 or joined to nothing).
#'  * *`sigma` and group-level SDs* follow `brms`, e.g., `student_t(3, 0, max(2.5, mad(y)))`, since
#'      half-Student-t priors are a standard weakly informative choice for scale parameters (Gelman, 2006).
#'
#' Like `rstanarm` and `brms`, which place intercept priors on centered predictors, the default prior on
#' the segment-1 intercept applies to the level at the start of the segment, `min(x)`. The reported
#' `Intercept_1` is still at `x = 0`, as in [lm()], so its prior is given the x-terms, e.g.,
#' `dnorm(mean(y) - x_1 * min(x), 2.5 * sd(y))` in `fit$prior`. Later segments' intercepts are at their change point.
#' Unlike in `rstanarm` and `brms`, other predictors are not centered, so the prior applies where they are zero.
#' If such a predictor is far from zero, center it before fitting or set a prior on `Intercept_1`.
#' User-specified intercept priors apply to `Intercept_1` as written. The autoscaling does not account for
#' segment length. If you expect steep slopes or strong curvature within short segments, consider setting a
#' wider prior.
#'
#' The default prior on change points is `dirichlet(1)` (uniform order statistics).
#' For a single change point, this is the Beta(1, 1) / Uniform distribution over `[min(x), max(x)]`.
#' For multiple change points, it corresponds to a flat Dirichlet distribution over segment lengths.
#' You can also explicitly set `cp_i = "dirichlet(alpha)"` with the same positive `alpha` for all
#' change points to regularize spacing (`alpha > 1` penalizes change points from occurring close together,
#' while `alpha < 1` favors clustering). Under the hood, `dirichlet(1)` is sampled as sorted uniforms
#' and other `alpha` as an exact sequential stick-breaking Beta chain for fast and robust sampling.
#' [Read more](https://lindeloev.github.io/mcp/articles/priors.html).
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
