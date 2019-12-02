# mcp 0.2
The API and internal structure should be stable now.

## New features: 

 * Quadratic and other terms using `I(x^2)`, `I(x^3.24)`, `sin(x)`, `sqrt(x)`, etc.
 * Model variance for `family = gaussian()` using `~ sigma([formula here])`.
 * Do order-N autoregressive models (AR(N)) using e.g.,`~ ar([order here])`. Useful for time series.
 * Plot prediction intervals using `plot(fit, quantiles = TRUE, quantiles_type = "predict")`.
 * Now respects `options(mc.cores = 3)`. All guides have been updated to recommend this as a default.

## Other changes:

 * `fit$func_y` has been renamed to `fit$simulate`.
 * `plot()` now only visualizes the total fit while `plot_pars()` only visualize individual parameters. These functions were mixed in `plot()` before.
 * The argument `update` has been discarded from `mcp()` (it's all on `adapt` now) and `inits` has been added.
 * Many internal changes to prepare for upcoming features. The biggest internal change is that `rjags` and `future` replace the `dclone` package. This gives faster and cleaner installations, is more future-proof, and avoids the need for a temporary file on the disk when sampling in parallel.
 * Updated documentation and website.
 

# mcp 0.1.0
First public release.

 * Varying change points
 * Basic GLM: Gaussian, binomial, Bernoulli, and Poisson, and associated vignettes.
 * summary(fit), fixef(fit), and ranef(fit)
 * plot(fit, "segments") and plot(fit, "bayesplot-name-here") with some options
 * 1000+ basic unit tests to ensure non-breaking code for a wide variety of models.
 * Testing and model comparison using `loo` and `hypothesis`
