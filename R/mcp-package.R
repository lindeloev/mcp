#' mcp: Regression with Multiple Change Points
#'
#' @description
#' `mcp` is regression with multiple change points, inferring full Bayesian uncertainty for all parameters, including
#' the change points themselves. At its simplest, it feels like `lm()` with one formula per segment. At its most
#' flexible, it aims to be **[`brms`](https://paulbuerkner.com/brms/) for change points**:
#'
#' * *GLMM and distributional regression:* Gaussian, binomial, Bernoulli, Poisson, and Negative Binomial families.
#'   Group-level (random) effects, AR/MA, and regression on distributional parameters (e.g., `sigma`, `shape`).
#' * *Express expert knowledge:* Specify the number of change points, the form of each segment, and constrain
#'   coefficients via priors. Compare hypotheses using `loo()` and `hypothesis()`.
#' * *Native feel:* `summary()`, `plot_pars()`, `pp_check()`, and `as_draws_df()` follow the conventions of
#'   [`brms`](https://paulbuerkner.com/brms/), [`rstanarm`](https://mc-stan.org/rstanarm/), and
#'   [`posterior`](https://mc-stan.org/posterior/articles/posterior.html), adapted to change point regression.
#'   Inspect the sampled JAGS model in `fit$jags_code`.
#'
#' You will find many aliases for change point regression in the literature, including piecewise linear,
#' discontinuity, segmented, and broken stick regression.
#'
#' See [mcp()] to fit models, [mcp-formula] for the model and formula syntax, and [mcp-priors] for priors.
#'
#' @section Extended formulas and model features - and links to worked examples:
#' * *Distributional regression:* Model residual variance and other distributional parameters across segments,
#'   e.g., `~ sigma(1 + x)` or `~ shape(1)`. [Read more](https://lindeloev.github.io/mcp/articles/dpar.html).
#' * *Time-series residuals:* Model serial dependence using `ar(p)` and `ma(q)` terms with a generalized link-scale
#'   recurrence. [Read more](https://lindeloev.github.io/mcp/articles/arma.html).
#' * *Group-level effects:* Add hierarchical (random) intercepts, slopes, and change points, e.g., `~ 1 + (1|id)`
#'   or `~ cp((1|id)) + 0 + x`. [Read more](https://lindeloev.github.io/mcp/articles/group_effects.html).
#' * *Families:* Binomial, Bernoulli, Poisson, and negative-binomial responses with several link functions.
#'   [Read more](https://lindeloev.github.io/mcp/articles/families.html).
#' * *Model comparison and hypothesis testing:* Compare models using PSIS-LOO (`loo()`) or WAIC (`waic()`), and test hypotheses
#'   (`hypothesis()`) using Savage-Dickey density ratios for point-null tests or posterior-to-prior odds for directional tests.
#'   [Read more](https://lindeloev.github.io/mcp/articles/comparison.html).
#' * *Fits, predictions, and forecasts:* `fitted()`, `predict()`, and `residuals()`, including forecasts beyond the data.
#'   [Read more](https://lindeloev.github.io/mcp/articles/predict.html).
#' * *Missing responses:* Missing values in the response are imputed jointly with the parameters.
#'   [Read more](https://lindeloev.github.io/mcp/articles/missing.html).
#'
#' @section Options:
#' * `mcp_quiet`: Logical. The default of `mcp(..., quiet)`. `options(mcp_quiet = TRUE)` suppresses routine JAGS output
#'   and mcp sampling-status messages for the session.
#' * `mcp_reparam`: Logical. `TRUE` (default) lets [mcp()] sample some default priors on scales that mix
#'   faster. This changes only `fit$jags_code` and sampling efficiency, not the model, the priors, or the posterior.
#'   `options(mcp_reparam = FALSE)` writes the plain model, which is mostly useful for seeing how `mcp` works.
#'   Ignored when `jags_code` is supplied. Read more in
#'   [From formulas to sampler code](https://lindeloev.github.io/mcp/articles/implementation.html#reparameterizations).
#'
#' @inherit mcp references
#'
#' @seealso \code{\link{mcp}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
"_PACKAGE"

