# mcp 0.4.0

## Major new features

-   **Multiple regression:** `mcp` now supports several continuous predictors, categorical predictors, interactions, etc. for all terms on RHS. E.g., `~ 1 + x + x:group + sigma(1 + group) + ar(2, 0 + z)`. Basically, it now aims to "feels" like `lm()` or `glm()` for each distributional parameter in each segment. Explore `ex = mcp_example("multiple")` to see it in action. Default priors generally align with brms, with some adjustments to accommodate the change-point model.

-   **Group-level effects on RHS:** Predictor formulas now support group-level effects (sometimes called random effects) using familiar `lme4` and `brms` syntax. `(1 | group)` specifies a group-level intercept, while `(1 + x || group)` and `(factor || group)` support independent coefficients, including slopes and factors. This also works inside distributional formulas such as `sigma(1 + (factor || id))`. As with `ar()`, an effect carries into later segments until it is redefined or disabled with `(0 | group)`. See `mcp_example("group_mu")` for a worked example. Correlated multi-coefficient terms are not yet supported.

-   `mcpfit`s now work natively with R generics, `{posterior}` and `{tidybayes}` posterior draw and prediction API. Changes include:
    - `summary(fit)` now reports rank-normalized split-`rhat`, `ess_bulk`, and `ess_tail` from `{posterior}` and central quantile intervals instead of HDIs.
    - Adding `as_draws(fit)`, `as_draws_df(fit)`, `tidy_draws(fit)`, `ndraws(fit)`, rstantools linpred/epred, etc. with full S3 generic registration for `{posterior}`, `{rstantools}`, and `{tidybayes}` (`spread_draws()`, `gather_draws()`). 
    - Per-draw methods (`summary = FALSE` in `fitted()`, `predict()`, `residuals()`, `log_lik()`) now return dot-prefixed columns (`.epred`, `.prediction`, `.residual`, `.loglik`) for `{tidybayes}` / `{ggdist}` compatibility.
    - `nsamples` soft-deprecated in favor of `ndraws`; `which_y` deprecated in favor of `dpar`.
    - Added R generics for mcpfits: `formula(fit)`, `family(fit)`, `model.frame(fit)`, `nobs(fit)`, `vcov(fit)`, and `confint(fit)`. Fits now also store a proper matched `$call`. Printed summaries now report posterior `sd` and has a new layout.

-   **Negative binomial and GARMA:** You can now do `mcp(..., family = negbinomial())`. Autoregression (`ar()`) has been generalized to GARMA link-scale residuals for Gaussian, binomial, Poisson, and negative-binomial models with their default links, using `ar(..., boundary = 0.1)` by default to keep zero and boundary counts finite. Added moving-average terms with `ma(q)`, which can be used alone or combined with `ar(p)` in each segment.

-   **Missing response imputation:** Missing responses are now retained as posterior JAGS imputations. `predict(fit) |> filter(is.na(y))` can be used to see imputed values - similarly for `fitted()`. It supports AR/MA too. Missing values in AR/MA are not supported in `log_lik()`, `loo()`, and `waic()` who will throw an informative error. See more details in the new vignette/article on missing data.


## Major breaking changes

mcp v0.4 is a major breaking change with the aim of remaining relatively stable going forward towards version 1.0. Although a lot has been updated, the parameter estimates in v0.4.0 remain practically identical to the previous public release (v0.3.4). Deprecation detections were added until we reach 1.0.

-   Renamed parameters to be more consistent with brms: `int_i` --> `Intercept_i`; `x_1_E2` --> `xE2_1`; `x_1_sin` --> `sinx_1`, etc.

-   `fit$pars` was removed. Use `mcp_pars(fit)` for a canonical table of model parameters and `mcp_columns(fit)` for resolved data-column roles, including the automatically chosen `par_x`. 

-   In general, renamed `"varying"` to `"group"`. The example names `"varying_mu"` and `"varying_cp"` are now `"group_mu"` and `"group_cp"`. In `plot_pars()`, use `pars = "group"`; the old `"varying"` selector remains as a deprecated alias. The established `varying =` method argument remains unchanged for now.

-   `summary()`, `fixef()`, `ranef()`, `prior_summary()`, and everything else now return rows in a canonical order instead of the previous incidental (near-alphabetical) order. Use `verbose = TRUE` with `summary()`, `fixef()`, or `ranef()` to include `segment` and `dpar` columns.

-   `Rhat` was renamed to `rhat` in `summary()`, `fixef()`, and `ranef()`, to match `{posterior}` standard.

-   `plot()` is now split into `plot()` for plotting full fits while `plot_dpar()` plots one distributional parameter (`mu`, `sigma`, `shape`, `ar1`, etc.). The argument order was changed too. The new coloring function (`plot(fit, color_by = "column")`) is particularly useful when models include categorical predictors or rhs group-level effects. See `mcp_example("group_mu")` for a worked example.

-   `hypothesis()` now restricts equality tests to simple parameter contrasts and gives clearer guidance on Savage-Dickey Bayes factors. Non-linear transformations in hypotheses like `x_1 / x_2 = 1` could previously be run, even though Savage-Dickey does not support it. Equality tests now return `p = NA` because their Bayes factor is a model comparison, not a posterior parameter probability. Density estimates are more robust and warn about sparse tails.

-   `fixef()` now reports only population-level effects for `mu` (the primary response parameter), excluding change points, other distributional parameters, AR/MA parameters, and group-effect SDs. Use `summary(fit)` to get all parameters.

-   `y | weights(w)` now specifies brms-aligned observation log-likelihood weights rather than Gaussian precision weights (which previously scaled the residual SD as `sigma / sqrt(w)`). Predictive draws and expectations now use `sigma` directly while `log_lik()` multiplies observation log-densities by `w`. Additionally, the weight column is no longer removed from `fitted()` and `predict()` output so that observation weights remain available in returned data.

-   AR and MA intercepts now have zero-centered, regularizing `dnorm(0, 0.5) T(-1, 1)` priors, replacing independent uniform priors. Their categorical contrasts and numeric slopes now use modest normal priors instead of heavy-tailed Student-t priors. Coefficients remain direct and are not jointly constrained to stationary or invertible regions.

-   Gaussian models with an explicit `sigma()` formula now model the residual standard deviation with a log link, as in `brms`. Models without an explicit `sigma()` are unchanged: `sigma_1` remains the residual SD itself with the response-scale half-Student-t prior.

-   Formula transformations now use the original predictor values, as in `lm()`, `glm()`, and `brms`. Previously, transformations of the change-point predictor such as `sin(x)` and `exp(x)` restarted at each segment onset. Bare `par_x` and polynomial bases such as `I(par_x^2)` remain segment-local to support joined segment shapes.

-   Group-level change points `(1|id)` remain exactly zero-centered to identify population-level change points and support efficient sampling. For models with multiple group-level change points, their realized locations are now constrained to remain ordered. Read more in the vignette on group-level effects.

-   After sampling, `mcp()` now verifies that population- and group-level change points are strictly ordered in every draw and that population-level change points lie within the observed predictor range. If ordering is broken, wrong segments apply.


-   Parallel sampling is now controlled exclusively through the active `{future}` plan. The `cores` argument to `mcp()` is deprecated and ignored, but remains available for backwards compatibility. Use `future::plan(future::multisession, workers = 3)` before calling `mcp()` to sample chains in parallel, and `future::plan(future::sequential)` to shut down those workers. Without a parallel future plan, chains are sampled sequentially.

-   `mcp(..., adapt = 1000)` is soft-deprecated in favor of `mcp(..., warmup = 1000)` for compatibility with additional future samplers.

-   Dropped support for `rel()` in formulas. This was ambiguous for interaction terms and made the code hard to maintain. Another way of achieving the same functionality via the priors may be added in future versions.

-   The arguments for `fit$simulate()` have all changed to accommodate multiple (categorical) predictors. `fit$simulate(fit, data, ..., .type = "predict")` is the new argument structure. Note that (1) it now requires `fit` as the first argument, (2) it requires `data.frame` or `tibble` as the second argument instead of just a vector of `par_x`, (3) further arguments are prefixed with a "." to avoid name conflicts internally in mcp. `...` are the model parameters as usual.

-   Dropped support for `mcp(..., data = NULL)`. You now must provide some mock-up data to inform `mcp` about the types and levels of the predictor columns. See, e.g., `mcp_example("intercepts")$example_code` for a simple example or `mcp_example("multiple")$example_code` for a more involved example. All docs have been updated appropriately.

-   mcp no longer exports `phi`, `logit`, `ilogit`, or `probit`.

-   `fit = mcp_example("name")` now returns the fit directly instead of a list with a `$fit` entry. It now defaults to sampling the model (`sample = "post"`) and the `sample` argument is now directly passed to `mcp(..., sample = sample)` so `sample = TRUE` is deprecated.


## Other new features

-   Added `prior_summary(fit)`. Its compact output shows each parameter's resolved prior and bounds; `prior_summary(fit, verbose = TRUE)` also shows the data-dependent rule, a plain-language description, source, and kind (`distribution`, `alias`, `expression`, or `constant`). Default priors are now resolved before JAGS code is generated, so `fit$prior` and generated code no longer depend on opaque data constants.

-   With the introduction of multiple regression including categorical predictors (including for distributional parameters) come default priors. For non-`par_x` terms, numeric coefficients are now auto-scaled to data and model. Briefly, based on its observed range when it has two values, and two standard deviations otherwise, aligning with [Gelman (2008)](https://doi.org/10.1002/sim.3107) and the prior autoscaling in [`rstanarm`](https://mc-stan.org/rstanarm/reference/priors.html). See `prior_summary(fit, verbose = TRUE)` for details.

-   Default `plot(fit)` style has been updated in many ways to accommodate multiple regression and group-effects.

-   Added `mcp(..., series = "data_column")` to identify independent series in models with AR/MA terms.

-   Added option to test hypotheses on the prior using `hypothesis(fit, prior = TRUE)`.

-   `mcp_example(...)` now shows an illustrative plot of the fitted example by default. Disable using `plot = FALSE`.

-   Methods like `fitted()` and `predict()` now accept `fitted(fit, varying = "cp")` and `fitted(fit, varying = "predictor")` as fast selectors for group-level effects in the corresponding formula part. Exact group-level parameter names remain supported too, `varying = TRUE` selects all, and `ranef()` continues to return all group-level effects.

-   Added `log_lik(fit)` which returns a draws-by-observation matrix by default. This is the usual `{brms}` return shape and is accepted directly by `{loo}`. `log_lik()` supports the same arguments as e.g. `fitted()`.

-   Memory improvement: The `mcpfit` is now \< 10% of the size as before because the log-likelihood is not stored. Use `log_lik(fit)` to compute it, or call `loo(fit)` or `waic(fit)` directly.

-   Sampling is now 1-10% faster due to a new formalization of the underlying JAGS code.

-   Use `mcp(..., seed = 42)` for reproducible JAGS sampling. See `mcp_example("demo")$example_code` how to ensure reproducibility across simulation, fit, and plotting.

-   `mcp(..., quiet = TRUE)` suppresses routine JAGS output and mcp sampling-status messages while preserving warnings and errors.

-   Added warnings when various thresholds are exceeded. `mcp(..., diagnostics = list(rhat = 1.01, ess_bulk = 400, ess_tail = 400, ar = 0.10, ma = 0.10))` controls fit diagnostics with configurable thresholds. Partial lists override these defaults, individual diagnostics can be disabled with `NULL`, and `diagnostics = FALSE` disables diagnostic warnings. `summary(fit)` inherits these settings and accepts its own override for the diagnostic footer.

-   Added AR/MA warnings: (1) for AR/MA models to `loo()`, `predict()`, etc. where the serial dependence is currently ignored. Proper handling requires leave-future-out or blocked cross-validation, which are not currently implemented in `mcp`. (2) when posterior violation probabilities exceed the configurable `ar` or `ma` diagnostic thresholds (10% by default).

-   Several new arguments to `loo`. `loo(fit, pointwise = TRUE)` uses `loo::loo.function()` for more memory-efficient (but slower) computation of LOO. Other new arguments include the usual from `fitted()` etc.: `loo(fit, ndraws = 1000, arma = FALSE, varying = FALSE)`.

-   Added `interpolate_newdata(fit, by = NULL)` which generates a data.frame with all combinations of categorical predictors along with interpolated continuous predictors. Use `by` to include grouping factors. The documentation shows how this can be useful for generating custom plots when simple tweaking `plot()` is not enough.

-   Fits made with custom `jags_code` now warn when calling R-side simulation, prediction, or model-evaluation methods, because those methods continue to evaluate the supplied formulas rather than the custom JAGS code.


## Minor breaking changes

-   mcp now requires R >= 4.1.0, matching its use of the native pipe (`|>`).

-   Corrected the JAGS translation of prior scales. `ddexp()` and `dlogis()` now convert their conventional scale to inverse scale rather than inverse variance; The exported `sd_to_prec()` helper is soft-deprecated because prior translation is now an internal, sampler-specific step.

-   Removed the non-functional `exponential()` family. It had no default priors and could not be fitted. Exponential survival models would require dedicated support for censoring and survival likelihoods rather than the previous ordinary-response placeholder.

-   Default coefficient priors are now more consistent with `brms` defaults while retaining proper priors for the `mcp` parameterization. For single-segment models, coefficient priors are directly aligned with `brms` / Gelman (2008) scaling, using `max(x) - min(x)` (or `2 * sd(z)`) as the reference predictor change across all segments. This ensures that slope priors remain identical and comparable across models regardless of the number of change points. For non-small datasets, this should have minimal influence since priors remain minimally informative. See priors using `prior_summary(fit, verbose = TRUE)`. The implicit Gaussian `sigma_1` prior is now the same response-calibrated half-Student-t used by `brms`, `dt(0, max(2.5, round(mad(y), 1)), 3) T(0, )`. Version 0.3.4 instead used a response-SD-calibrated half-normal prior. Both are weakly informative so will have negligible impact on fits.

-   Uppercase data-dependent constants for priors (`MINX`, `MAXX`, etc.) remain temporarily supported with a deprecation warning. They are now replaced with readable expressions such as `min(time)`, `max(time) - min(time)`, `mad(response)`, `n_segments()`, and `n_cp()`.

-   Binomial and Bernoulli models with logit or probit links now use a narrower, weakly regularizing `dt(0, 1.5, 3)` rather than `dt(0, 2.5, 3)` for intercepts and categorical contrasts.

-   Gaussian models with a log link is now more aligned with `brms`: calibrate priors using data, replacing zeros with 0.1 and using finite location and scale fallbacks.

-   The term "ct" (for "central tendency") has been replaced with "mu", e.g., in `plot(fit, which_y = "mu")`. These were defaults, so hopefully no-one will notice the renaming.

-   `fit$data` now only contains the data columns that are used in the model.

-   Removed `which_y` argument from `predict()`.

-   New `mcpfit` objects no longer include empty `$loo` and `$waic` components.

-   `draws_format` replaces `samples_format` in `fitted()`, `predict()`, and `log_lik()`. The old argument remains available with a soft deprecation.

## Bug fixes

-   Formula environments are now preserved, so locally defined transformations work in `mcp()` formulas and when applying fitted predictors to new data.

-   `hypothesis(..., prior = TRUE)` now reports `BF = NA` rather than comparing prior draws with themselves and returning a Bayes factor of 1.

-   Directional Bayes factors now divide posterior odds by prior odds. Previously, they reported posterior odds alone, which is only a Bayes factor when the prior probability of the hypothesis is 0.5.

-   `residuals()` now uses the conventional `observed - fitted` sign. Previously, `mcp` residuals used `fitted - observed`.

-   Fixed several posterior predictive check and information-criterion bugs present in v0.3.4. From most to least serious:
    - For missing data: Missing responses were scored in WAIC/LOO as if their JAGS-imputed values had been observed. Fix: missing responses remain latent in JAGS but are excluded from observed-data PPC, log-likelihood, WAIC, and LOO calculations.
    
    - For models with group-level change points, LOO checks with `facet_by != NULL` could pair group-specific predictions with weights from the wrong observations.
    - For AR/MA models, `pp_check()` conditioned every prediction on the observed response history. They now generate each posterior replication recursively as a fresh time series.
    - LOO checks with `prior = TRUE` incorrectly combined prior predictions with posterior LOO weights. A rare use case.
    - `pp_check(..., ndraws = NULL)` errored. Fixed.

-   Weighted regression: While JAGS correctly modelled weights, R-side simulation/generation ignored it. This means `predict()`, PPCs, `log_lik()`, WAIC, and LOO were incorrect when using weighted regression. Weighted Gaussian posterior predictions and log-likelihoods now use the observation-level standard deviation `sigma / sqrt(weight)`, matching the JAGS precision `weight / sigma^2`. This makes `predict()`, posterior predictive checks, `log_lik()`, WAIC, and LOO consistent with the fitted model.

-   In models going from higher-order to lower-order, (`~ ar(2), ~ ar(1)`), the higher-order components were not "turned off".

-   Bug only noticeable for very small samples: For Gaussian identity-link AR models, R-side calculations could leak residuals between posterior draws and omitted available partial lags for the first observations of AR(2+) models. Each draw is now evaluated as a separate series and uses the same partial-lag recurrence as JAGS.

-   Conditional R-side GARMA histories are now evaluated across posterior draws rather than through a scalar loop over every draw and observation. This substantially speeds up `fitted()`, `predict()`, `residuals()`, and `log_lik()` for AR/MA models.

-   The quantiles from `fitted()` and `predict()` for group-level change-point models ignored the group level and were identical across levels.

-   Now works for 200+ characters formulas too.

-   Fixed `cores = "all"` failed. Thanks for reporting, @m-r-munroe!

## Behind the scenes

In general, every effort has been made to anticipate future developments and build the structure needed for that.

-  The `mcpfamily` contains everything which was previously scattered across multiple functions. This means it is easy to add families and link functions - also enabling custom families in a future release.

-   Major changes in how the model is translated into JAGS code. The JAGS code is quite different but functionally equivalent. This also opens for adding a new sampler in a future release.

-   More thorough defensive coding.

-   Much expanded test suite (now 7.000+ tests when run in full). The test suite now includes external validation of inference and simulation: AR against `arima()`/`arima.sim()`, binomial against `glm()` / `rbinom()`, and changepoints against `{segmented}`.

-   Fewer imports to userspace. This minimizes the risk of name conflicts.


# mcp 0.3.4

This is a bug fix release.

## Bug fixes

-   Now respects the `cores` argument to `mcp()`.
-   Document all function arguments and remove documentation for removed arguments.

# mcp 0.3.3

This is a bug fix release.

## Bug fixes

-   Support `ggplot >= 3.4.0`, `tidyselect >= 1.2.0`, and newer `future` by replacing deprecated functions.
-   Accept `mcp(..., cores = "all")`.
-   Fix documentation of `iter` argument to `mcp()`.
-   Other small fixes to deployment and documentation.

# mcp 0.3.2

This release contains no user-facing changes. The test suite suite is now compatible with dplyr 1.0.8, which caused the test suite to fail. This, in turn, would trigger the removal of mcp from CRAN.

# mcp 0.3.1

This is mostly a bug fix release.

## New features:

-   `ex = mcp_example("demo", sample = TRUE)` is the new interface that replaces the `ex_*` datasets in prior versions. This reduces clutter of the namespace/documentation and the size of the package. It also gives the user richer details on the simulation and analyses. For "demo", the `ex_demo` dataset is now `ex$data` and the `ex_fit` is `ex$fit`.

-   Nicer printing of lists and texts all over. E.g., try `print(demo_fit$jags_code)` and `print(demo_fit$pars)`.

## Bug fixes

-   Support breaking changes in `tidybayes >= 3.0.0` and `dplyr >= 1.0.6`

# mcp 0.3.0

## New features:

-   Get fits and predictions for in-sample and out-of-sample data. [Read more in the article on these functions](https://lindeloev.github.io/mcp/articles/predict.html).

    -   Use `predict(fit)` to get predicted values and quantiles.
    -   Use `fitted(fit)` to get estimated values and quantiles.
    -   Use `residuals(fit)` to get residuals and quantiles.

    All of the above functions include many arguments that align with (and extends) the options already in `plot.mcpfit()`, including getting fits/predictions for sigma (`which_y = "sigma"`), for the prior (`prior = TRUE`), and arbitrary quantiles (`probs = c(0.1, 0.5, 0.999)`). Use the `newdata` argument to get out-of-sample fitted/predicted values. Set `summary = FALSE` to get per-draw values.

-   Added support for weighted regression for gaussian families: `model = list(y | weights(weight_column) ~ 1 + x)`. Weights are visualized as dot sizes in `plot(fit)`.

-   Support for more link functions across families (e.g., `family = gaussian(link = "log")`):

    -   `gaussian`: "identity", "log"
    -   `binomial`: "logit", "probit", "identity"
    -   `bernoulli`: "logit", "probit", "identity"
    -   `poisson`: "log", "identity"

-   New argument `scale` in `fitted()`, `plot()`, and `fit$simulate()`. When `scale = "response"` (default), they return fits on the response scale. When `scale = "linear"`, they return fits on the linear-predictor (link) scale. Useful for model understanding and debugging.

-   Use `pp_check(fit)` to do prior/posterior predictive checking. See `pp_check(fit, type = "x")` for a list of plot types. `pp_check(fit, facet_by = "group_column")` facets by a grouping column.

-   Improvements to `plot()`:

    -   Change-point densities are now computed on a per-panel basis in `plot(fit, facet_by = "group_column")`. Previous releases only displayed population-level change points.
    -   You can now plot group-level effects with `rate = FALSE` for binomial models.
    -   Change point densities in `plot(fit)` are not located directly on the x-axis. They were "floating" 5% above the x-axis in the previous releases.

-   New argument `nsamples` reduces the number of samples used in most functions to speed up processing. `nsamples = NULL` uses all samples for maximum accuracy.

-   New argument `arma` in many functions toggles whether autoregressive effects should be modelled.

-   Although the API is still in alpha, feel free to try extracting samples using `mcp:::tidy_samples(fit)`. This is useful for further processing using `tidybayes`, `bayesplot`, etc. and is used extensively internally in `mcp`. One useful feature is computing absolute values for group-specific change points: `mcp:::tidy_samples(fit, population = FALSE, absolute = TRUE)`. Feedback is appreciated before `tidy_samples` becomes part of the `mcp` API in a future release.

## Other changes

-   Change point densities in `plot(fit)` are now scaled to 20% of the plot for each chain X changepoint combo. This addresses a common problem where a wide posterior was almost invisibly low when a narrow posterior was present. This means that heights should only be compared *within* each chain x changepoint combo - not across.
-   Removed the implicit ceiling of 1000 lines and samples in `plot.mcpfit()`.
-   Rownames are removed from `ranef()` and `fixef()` returns.
-   A major effort has been put into making `mcp` robust and agile to develop. `mcp` now use defensive programming with helpful error messages. The Test suite includes 3600+ tests.
-   `plot()`, `predict()`, etc. are now considerably faster for AR(N) due to vectorization of the underlying code.

## Bug fixes

-   Sigma is now forced to stay positive via a floor at 0.
-   Fixed: support and require dplyr 1.0.0. Now also requires tidybayes 2.0.3.
-   Fixed: Parallel sampling sometimes produced identical chains.
-   Fixed several small bugs

# mcp 0.2.0

The API and internal structure should be stable now. v0.2.0 will be released on CRAN.

## New features:

-   Model quadratic and other terms using `I(x^2)`, `I(x^3.24)`, `sin(x)`, `sqrt(x)`, etc.

-   Model residual standard deviation for `family = gaussian()` using `~ sigma([formula here])`.

-   Model Nth order autoregressive models using `~ ar(order, formula)`, typically like `y ~ 1 + x + ar(2)` for AR(2). Simulate AR(N) models from scratch or given known data with `fit$simulate()`. The [article on AR(N)](https://lindeloev.github.io/mcp/articles/arma.html) has more details and examples. AR(N) models are popular to detect changes in time-series.

-   Many updates to `plot()`.

    -   Includes the posterior densities of the change point(s). Disable using `plot(fit, cp_dens = FALSE)`.
    -   Supports AR(N) models (see above).
    -   Plot posterior parameter intervals using `plot(fit, q_fit = TRUE)`. `plot(fit, q_fit = c(0.025, 0.5, 0.975))` plots 95% HDI and the median.
    -   Plot posterior predictive intervals using `plot(fit, q_predict = TRUE)`.
    -   Choose data geom. Currently takes "point" (default) and "line" (`plot(fit, geom_data = "line")`). The latter is useful for time series. Disable using `geom_data = FALSE`.

-   Use `options(mc.cores = 3)` for considerable speed gains for the rest of the session. All vignettes/articles have been updated to recommend this as a default, though serial sampling is still the technical default. `mcp(..., cores = 3)` does the same thing on a call-by-ball basis.

-   `fit$simulate()` adds the simulation parameters as an attribute (`attr(y, "simulate")`) to the predicted variable. `summary()` recognizes this and adds the simulated values to the results table (columns `sim` and `match`) so that one can inspect whether the values were recovered.

-   Use `plot(fit, which_y = "sigma")` to plot the residual standard deviation on the y-axis. It works for AR(N) as well, e.g., `which_y = "ar1"`, `which_y = "ar2"`, etc. This is useful to visualize changes in residual standard deviation and autocorrelation. The relevant vignettes have been updated with worked examples.

-   Much love for the priors:

    -   Set a Dirichlet prior on the change points using `prior = list(cp_1 = "dirichlet(1)", cp_2 = ...)`. [Read pros and cons here](https://lindeloev.github.io/mcp/articles/priors.html).
    -   The default prior has been changed from truncated uniforms to a regularizing t-tail prior that samples effectively. [Read more here](https://lindeloev.github.io/mcp/articles/priors.html).
    -   You can now sample the prior using `mcp(..., sample = "prior")` or `mcp(..., sample = "both")` and most methods can now take the prior: `plot(fit, prior = TRUE)`, `plot_pars(fit, prior = TRUE)`, `summary(fit, prior = TRUE)`, `ranef(fit, prior = TRUE)`.

-   `mcp` can now be cited! Call `citation("mcp")` or see the pre-print here: <https://osf.io/fzqxv>.

## Other changes:

-   Some renaming: "segments" --> "model". `fit$func_y()` --> `fit$simulate()`.
-   `plot()` only visualize the total fit while `plot_pars()` only visualize individual parameters. These functions were mixed in `plot()` previously.
-   The argument `update` has been discarded from `mcp()` (it's all on `adapt` now) and `inits` has been added.
-   Many internal changes to make `mcp` more future proof. The biggest internal change is that `rjags` and `future` replace the `dclone` package. Among other things, this gives faster and cleaner installations.
-   Many more informative error messages to help you quickly understand and solve errors.
-   Updated documentation and website.

# mcp 0.1.0

First public release.

-   Varying change points
-   Basic GLM: Gaussian, binomial, Bernoulli, and Poisson, and associated vignettes.
-   summary(fit), fixef(fit), and ranef(fit)
-   plot(fit, "segments") and plot(fit, "bayesplot-name-here") with some options
-   1000+ basic unit tests to ensure non-breaking code for a wide variety of models.
-   Testing and model comparison using `loo` and `hypothesis`
