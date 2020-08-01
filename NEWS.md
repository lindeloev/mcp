# mcp 0.3.0

## New features:

 * Get fits and predictions for in-sample and out-of-sample data. [Read more in the article on these functions](https://lindeloev.github.io/mcp/articles/predict.html).
   - Use `predict(fit)` to get predicted values and quantiles. 
   - Use `fitted(fit)` to get estimated values and quantiles. 
   - Use `residuals(fit)` to get residuals and quantiles.
   
   All of the above functions include many arguments that align with (and extends) the options already in `plot.mcpfit()`, including getting fits/predictions for sigma (`which_y = "sigma"`), for the prior (`prior = TRUE`), and arbitrary quantiles (`probs = c(0.1, 0.5, 0.999)`). Use the `newdata` argument to get out-of-sample fitted/predicted values. Set `summary = FALSE` to get per-draw values.
 
 * Added support for weighted regression for gaussian families: `model = list(y | weights(weight_column) ~ 1 + x)`. Weights are visualized as dot sizes in `plot(fit)`.

 * Support for more link functions across families (e.g., `family = gaussian(link = "log")`):
   - `gaussian`: "identity", "log"
   - `binomial`: "logit", "probit", "identity"
   - `bernoulli`: "logit", "probit", "identity"
   - `poisson`: "log", "identity"

 * New argument `scale` in `fitted()`, `plot()`, and `fit$simulate()`. When `scale = "response"` (default), they return fits on the observed scale. When `scale = "linear"`, they return fits on the parameter scale where the linear trends are. Useful for model understanding and debugging.
 
 * Use `pp_check(fit)` to do prior/posterior predictive checking. See `pp_check(fit, type = "x")` for a list of plot types. `pp_check(fit, facet_by = "varying_column")` facets by a data column.

 * Improvements to `plot()`:
   - Change point densities are now computed on a per-panel basis in `plot(fit, facet_by = "varying_column")`. Previous releases only displayed population-level change points.
   - You can now plot varying effects with `rate = FALSE` for binomial models.
   - Change point densities in `plot(fit)` are not located directly on the x-axis. They were "floating" 5% above the x-axis in the previous releases.
 
 * New argument `nsamples` reduces the number of samples used in most functions to speed up processing. `nsamples = NULL` uses all samples for maximum accuracy.
 
 * New argument `arma` in many functions toggles whether autoregressive effects should be modelled.
 
 * Although the API is still in alpha, feel free to try extracting samples using `mcp:::tidy_samples(fit)`. This is useful for further processing using `tidybayes`, `bayesplot`, etc. and is used extensively internally in `mcp`. One useful feature is computing absolute values for varying change points: `mcp:::tidy_samples(fit, population = FALSE, absolute = TRUE)`. Feedback is appreciated before `tidy_samples` will to become part of the `mcp` API in a future release.


## Other changes

 * Change point densities in `plot(fit)` are now scaled to 20% of the plot for each chain X changepoint combo. This adresses a common problem where a wide posterior was almost invisibly low when a narrow posterior was present. This means that heights should only be compared *within* each chain x changepoint combo - not across.
 * Removed the implicit ceiling of 1000 lines and samples in `plot.mcpfit()`.
 * Rownames are removed from `ranef()` and `fixef()` returns.
 * A major effort has been put into making `mcp` robust and agile to develop. `mcp` now use defensive programming with helpful error messages. The Test suite includes 3600+ tests.
 * `plot()`, `predict()`, etc. are now considerably faster for AR(N) due to vectorization of the underlying code. 


## Bug fixes

 * Sigma is now forced to stay positive via a floor at 0.
 * Fixed: support and require dplyr 1.0.0. Now also requires tidybayes 2.0.3.
 * Fixed: Parallel sampling sometimes produced identical chains.
 * Fixed several small bugs


# mcp 0.2.0
The API and internal structure should be stable now. v0.2.0 will be released on CRAN.

## New features: 

 * Model quadratic and other terms using `I(x^2)`, `I(x^3.24)`, `sin(x)`, `sqrt(x)`, etc.
 * Model variance for `family = gaussian()` using `~ sigma([formula here])`.
 * Model Nth order autoregressive models using `~ ar(order, formula)`, typically like `y ~ 1 + x + ar(2)` for AR(2). Simulate AR(N) models from scratch or given known data with `fit$simulate()`. The [article on AR(N)](https://lindeloev.github.io/mcp/articles/arma.html) has more details and examples. AR(N) models are popular to detect changes in time-series.
 * Many updates to `plot()`.
   - Includes the posterior densities of the change point(s). Disable using `plot(fit, cp_dens = FALSE)`.
   - Supports AR(N) models (see above).
   - Plot posterior parameter intervals using `plot(fit, q_fit = TRUE)`. `plot(fit, q_fit = c(0.025, 0.5, 0.975))` plots 95% HDI and the median.
   - Plot prediction intervals using `plot(fit, q_predict = TRUE)`.
   - Choose data geom. Currently takes "point" (default) and "line" (`plot(fit, geom_data = "line")`). The latter is useful for time series. Disable using `geom_data = FALSE`.
 * Use `options(mc.cores = 3)` for considerable speed gains for the rest of the session. All vignettes/articles have been updated to recommend this as a default, though serial sampling is still the technical default. `mcp(..., cores = 3)` does the same thing on a call-by-ball basis.
 * `fit$simulate()` adds the simulation parameters as an attribute (`attr(y, "simulate")`) to the predicted variable. `summary()` recognizes this and adds the simulated values to the results table (columns `sim` and `match`) so that one can inspect whether the values were recovered.
 * Use `plot(fit, which_y = "sigma")` to plot the residual standard deviation on the y-axis. It works for AR(N) as well, e.g., `which_y = "ar1"`, `which_y = "ar2"`, etc. This is useful to visualize change points in variance and autocorrelation. The vignettes on variance and autocorrelations have been updated with worked examples.
 * Much love for the priors:
   - Set a Dirichlet prior on the change points using `prior = list(cp_1 = "dirichlet(1)", cp_2 = ...)`. [Read pros and cons here](https://lindeloev.github.io/mcp/articles/priors.html).
   - The default prior has been changed from "truncated-uniforms" to a "t-tail" prior to be more uninformative while still sampling effectively. [Read more here](https://lindeloev.github.io/mcp/articles/priors.html).
   - You can now sample the prior using `mcp(..., sample = "prior")` or `mcp(..., sample = "both")` and most methods can now take the prior: `plot(fit, prior = TRUE)`, `plot_pars(fit, prior = TRUE)`, `summary(fit, prior = TRUE)`, `ranef(fit, prior = TRUE)`.
 * `mcp` can now be cited! Call `citation("mcp")` or see the pre-print here: [https://osf.io/fzqxv](https://osf.io/fzqxv).

## Other changes:
 * Some renaming: "segments" --> "model". `fit$func_y()` --> `fit$simulate()`.
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
