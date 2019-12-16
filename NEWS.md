
# mcp 0.2.0
The API and internal structure should be stable now. v0.2.0 will be released on CRAN.

## New features: 

 * Model quadratic and other terms using `I(x^2)`, `I(x^3.24)`, `sin(x)`, `sqrt(x)`, etc.
 * Model variance for `family = gaussian()` using `~ sigma([formula here])`.
 * Model Nth order autoregressive models using `~ ar(order, formula)`, typically like `y ~ 1 + x + ar(2)` for AR(2). `plot()` visualize posteriors for AR(N) models. Simulate AR(N) models from scratch or given known data with `fit$simulate()`. The [article on AR(N)](https://lindeloev.github.io/mcp/articles/arma.html) has more details and examples. AR(N) models are popular to detect changes in time-series.
 * Plot prediction intervals using `plot(fit, quantiles = TRUE, quantiles_type = "predict")`.
 * Use `options(mc.cores = 3)` for considerable speed gains. All vignettes/articles have been updated to recommend this as a default, though serial sampling is still the technical default.
 * `fit$simulate()` adds the simulation parameters as an attribute (`attr(y, "simulate")`) to the predicted variable. `summary()` recognizes this and adds the simulated values to the results table (columns `sim` and `match`) so that one can inspect whether the values were recovered.
 * Use `plot(fit, which_y = "sigma")` to plot the residual standard deviation on the y-axis. It works for AR(N) as well (`which_y = "ar1"`, `which_y = "ar2"`, etc.). This is useful to visualize change points in variance and autocorrelation. The vignettes on variance and autocorrelations have been updated with worked examples.
 * Set a Dirichlet prior on the change points using `prior = list(cp_1 = "dirichlet(1)", cp_2 = ...)`. [Read pros and cons here](https://lindeloev.github.io/mcp/articles/priors.html).

## Other changes:

 * `fit$func_y()` has been renamed to `fit$simulate()`.
 * `plot()` only visualize the total fit while `plot_pars()` only visualize individual parameters. These functions were mixed in `plot()` previously.
 * The argument `update` has been discarded from `mcp()` (it's all on `adapt` now) and `inits` has been added.
 * Many internal changes to make `mcp` more future proof. The biggest internal change is that `rjags` and `future` replace the `dclone` package. Among other things, this gives faster and cleaner installations.
 * Many more informative error messages to help you quickly understand and solve errors.
 * Updated documentation and website.
 

# mcp 0.1.0
First public release.

 * Varying change points
 * Basic GLM: Gaussian, binomial, Bernoulli, and Poisson, and associated vignettes.
 * summary(fit), fixef(fit), and ranef(fit)
 * plot(fit, "segments") and plot(fit, "bayesplot-name-here") with some options
 * 1000+ basic unit tests to ensure non-breaking code for a wide variety of models.
 * Testing and model comparison using `loo` and `hypothesis`
