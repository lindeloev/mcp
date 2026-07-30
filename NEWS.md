# mcp 0.4.0

## Major new features

-   Supports several continuous predictors, categorical predictors, interactions, etc. for all terms on RHS. E.g., `~ 1 + x + x:group + sigma(1 + group) + ar(2, 0 + z)`. Basically, it now "feels" like `lm()` for each distributional parameter in each segment. All mcp functions support this now, including `plot()`, `fit$simulate()`, `predict(fit, newdata = ...)`, `hypothesis()`, `pp_check()`, etc. Explore `ex = mcp_example("multiple")` to see it in action. Default priors generally align with brms, with some adjustments to the change-point model.

-   Predictor formulas now support group-level effects (sometimes called random or varying effects) using familiar `lme4` and `brms` syntax. `(1 | group)` specifies a group-level intercept, while `(1 + x || group)` and `(factor || group)` support independent coefficients, including slopes and factors. This also works inside distributional formulas such as `sigma(1 + (factor || id))`. As with `ar()`, an effect carries into later segments until it is redefined or disabled with `(0 | group)`. See `mcp_example("group")` for a worked example. Correlated multi-coefficient `|` terms are not yet supported.

-   Prediction methods now accept `varying = "cp"` and
    `varying = "predictor"` as fast selectors for group-level effects in the
    corresponding formula part. Exact varying-parameter names remain
    supported too, `varying = TRUE` selects all, and `ranef()` continues to return
    all group-level effects.

## Major breaking changes

-   AR and MA intercepts now have zero-centered, regularizing
    `dnorm(0, 0.5) T(-1, 1)` priors, replacing independent uniform priors.
    Their categorical contrasts and numeric slopes now use modest normal priors
    instead of heavy-tailed Student-t priors. Coefficients remain direct and
    are not jointly constrained to stationary or invertible regions.

-   Gaussian models with an explicit `sigma()` formula now model the residual
    standard deviation with a log link, as in `brms`. Models without an explicit
    `sigma()` are unchanged: `sigma_1` remains the residual SD itself with the
    response-scale half-Student-t prior.

-   Formula transformations now use the original predictor values, as in
    `lm()`, `glm()`, and `brms`. Previously, transformations of the change-point
    predictor such as `sin(x)` and `exp(x)` restarted at each segment onset.
    Bare `par_x` and polynomial bases such as `I(par_x^2)` remain segment-local
    to support joined segment shapes.

-   Posterior summaries now report central quantile intervals instead of
    highest-density intervals. Sampling diagnostics now use the rank-normalized
    split-Rhat, bulk ESS, and tail ESS from `posterior`; the `n.eff` column has
    been replaced by `ess_bulk` and `ess_tail`. Intervals will be highly similar
    for near-symmetrical posteriors.

-   Parallel sampling is now controlled exclusively through the active `{future}` plan. The `cores` argument to `mcp()` is deprecated and ignored, but remains available for backwards compatibility. Use `future::plan(future::multisession, workers = 3)` before calling `mcp()` to sample chains in parallel, and `future::plan(future::sequential)` to shut down those workers. Without a parallel future plan, chains are sampled sequentially.

-   Dropped support for `rel()` in formulas. This was ambiguous for interaction terms and made the code hard to maintain. Another way of achieving the same functionality via the priors may be added in future versions.

-   Renamed parameters to be more consistent with brms: `int_i` --> `Intercept_i`; `x_1_E2` --> `xE2_1`; `x_1_sin` --> `sinx_1`, etc.

-   The arguments for `fit$simulate()` have all changed to accommodate multiple (categorical) predictors. `fit$simulate(fit, data, ..., .type = "predict")` is the new argument structure. Note that (1) it now requires `fit` as the first argument, (2) it requires `data.frame` or `tibble` as the second argument instead of just a vector of `par_x`, (3) further arguments are prefixed with a "." to avoid name conflicts internally in mcp. `...` are the model parameters as usual.

-   Dropped support for `mcp(..., data = NULL)`. You now must provide some mock-up data to inform `mcp` about the types and levels of the predictor columns. See, e.g., `mcp_example("intercepts")$call` for a simple example or `mcp_example("multiple")$call` for a more involved example. All docs have been updated appropriately.

-   Changed arguments and their order in `plot()`.

    -   The most used arguments are first and the new `color_by` has a nice place in the sequence of arguments. From `plot(fit, lines, geom_data, cp_dens, q_fit, q_predict, ...)` to `plot(fit, q_fit, q_predict, facet_by, color_by, lines, geom_data, ...)`

    -   Dropped arguments `which_y`, `scale` which only made sense for distributional parameters (dpars) and where several of the other arguments (e.g., `q_predict`) were insensible. Use the new `plot_dpar()` for this.

-   `which_y` is deprecated across `fitted()`, `predict()`, `plot()`, `fit$simulate()`, etc. Use `dpar` instead, following the naming used by `brms`/`posterior`. Its default also changed from `dpar = "mu"` to `dpar = "epred"`, i.e., `fitted()` now returns the expected value of the response by default rather than the central-tendency parameter on the response scale. For the currently supported families these coincide, but `epred` is the forward-compatible choice for families where the mean depends on more than one distributional parameter (e.g., lognormal).

-   `fit = mcp_example("name")` now returns the fit directly instead of a list with a `$fit` entry. It now defaults to sampling the model (`sample = "post"`) and the `sample` argument is now directly passed to `mcp(..., sample = sample)` so `sample = TRUE` is deprecated.

## Other new features


-   The user-facing `nsamples` argument is soft-deprecated in favor of `ndraws`, matching the terminology used by `posterior`, `rvar`, and related packages. Existing calls continue to work for now.

-   Extended autoregression (`ar()`) to GARMA link-scale residuals for Gaussian, binomial, Poisson, and negative-binomial models with their default links, using `ar(..., boundary = 0.1)` by default to keep zero and boundary counts finite. Added moving-average terms with `ma(q)`, which can be used alone or combined with `ar(p)` in each segment. Bernoulli models and non-default links remain unsupported.

-   Added AR/MA warnings: (1) for AR/MA models to `loo()`, `predict()`, etc. where the serial dependence is currently ignored. Proper handling requires leave-future-out or blocked cross-validation, which are not currently implemented in `mcp`. (2) when posterior probability is >10% of violation of AR-stationarity or MA-invertibility.

-   In addition to (segment-wide) intercepts and slopes, there are now default priors for categorical predictors.

-   Memory improvement: The `mcpfit` is now \< 10% of the size as before because the log-likelihood is not computed by default anymore (no `fit$mcmc_loglik` anymore). You can add it using `fit = add_loglik(fit)` (adds `fit$loglik`) but if absent, it is automatically computed when calling relevant functions, e.g., `loo(fit)`.

-   Several new arguments to `loo`. `loo(fit, pointwise = TRUE)` uses `loo::loo.function()` for more memory-efficient (but slower) computation of LOO. Other new arguments include the usual from `fitted()` etc.: `loo(fit, ndraws = 1000, arma = FALSE, varying = FALSE)`.

-   Sampling is now 1-10% faster due to a new formalization of the underlying JAGS code.

-   `plot()` now uses a filled area for the change-point densities to visually distinguish it from color-coded fitted lines.

-   Added `interpolate_newdata(fit, by = NULL)` which generates a data.frame with all combinations of categorical predictors along with interpolated continuous predictors. Use `by` to include varying-effect groups. The documentation shows how this can be useful for generating custom plots when simple tweaking `plot()` is not enough.

-   Added `niterations(fit)` and `nchains(fit)` for convenience.

-   Added `log_lik(fit)` which is analogous to e.g., `fitted(fit)`.

-   Added `prior_summary(fit)`. Its compact output shows each parameter's resolved prior and bounds; `prior_summary(fit, verbose = TRUE)` also shows the data-dependent rule, a plain-language description, source, and kind (`distribution`, `alias`, `expression`, or `constant`). Default priors are now resolved before JAGS code is generated, so `fit$prior` and generated code no longer depend on opaque data constants.

-   Added option to test hypotheses on the prior using `hypothesis(fit, prior = TRUE)`.

-   `mcp()` now warns when posterior samples show poor mixing.

-   With the introduction of multiple regression including categorical predictors (including for distributional parameters) come default priors. For term involving local `par_x`, this scale is multiplied by the expected segment width "as usual". For non-`par_x` terms, numeric coefficients are scaled to data and model. Briefly, based on its observed range when it has two values, and two standard deviations otherwise, aligning with [Gelman (2008)](https://doi.org/10.1002/sim.3107) and the prior autoscaling in [`rstanarm`](https://mc-stan.org/rstanarm/reference/priors.html).


## Minor breaking changes

-   `residuals()` now uses the conventional `observed - fitted` sign. Previously, `mcp` residuals used `fitted - observed`.

-   Removed the non-functional `exponential()` family. It had no default priors and could not be fitted. Exponential survival models would require dedicated support for censoring and survival likelihoods rather than the previous ordinary-response placeholder.

-   Default coefficient priors are now more consistent with `brms` defaults while retaining proper priors and change-point-aware scaling for the `mcp` parameterization. For non-small datasets, this should have minimal influence since priors remain minimally informative. See priors using `prior_summary(fit, verbose = TRUE)`.

-   Binomial and Bernoulli models with logit or probit links now use a narrower (but still very uninformative)
    `dt(0, 1.5, 3)` rather than `dt(0, 2.5, 3)` for intercepts and categorical
    contrasts.

-   Gaussian models with a log link is now more aligned with `brms`: calibrate priors using data, replacing zeros with 0.1 and using finite location and scale fallbacks.

-   The implicit Gaussian `sigma_1` prior is now the same response-calibrated half-Student-t used by `brms`, `dt(0, max(2.5, round(mad(y), 1)), 3) T(0, )`. Version 0.3.4 instead used a response-SD-calibrated half-normal prior. Both are weakly informative so will have negligible impact on fits.

-   Uppercase data-dependent constants for priors (`MINX`, `MAXX`, etc.) remain temporarily supported with a deprecation warning. They are now replaced with readable expressions such as `min(time)`, `max(time) - min(time)`, `mad(response)`, `segment_width(time)`, `n_segments()`, and `n_cp()`.

-   Default priors for coefficients involving local `par_x` (except AR/MA) now broaden with
    more change points because narrower intervals allow steeper slopes within
    the same overall spread.

-   `summary()`, `fixef()`, `ranef()`, `prior_summary()`, and everything else now return rows in
    a canonical order instead of the previous incidental (near-alphabetical)
    order. They also gained columns `segment` and `dpar` columns. This may affect scripts that index results positionally.

-   The term "ct" (for "central tendency") has been replaced with "mu", e.g., in `plot(fit, which_y = "mu")`. These were defaults, so I hope no one will notice the renaming.

-   `fit$data` now only contains the data columns that are used in the model.

-   Removed `which_y` argument from `predict()`.

## Bug fixes

-   Directional Bayes factors now divide posterior odds by prior odds. Previously,
    they reported posterior odds alone, which is only a Bayes factor when the
    prior probability of the hypothesis is 0.5.

-   Fixed several posterior predictive check and information-criterion bugs present in v0.3.4. From most to least serious:
    - For missing data: Missing responses were scored in WAIC/LOO as if their JAGS-imputed values had been observed. Fix: missing responses remain latent in JAGS but are excluded from observed-data PPC, log-likelihood, WAIC, and LOO calculations.
    - For models with varying change points: LOO checks with `facet_by != NULL` could pair group-specific predictions with weights from the wrong observations.
    - LOO checks with `prior = TRUE` incorrectly combined prior predictions with posterior LOO weights. A rare use case.
    - `pp_check(..., ndraws = NULL)` errored. Fixed.

-   Weighted regression: While JAGS correctly modelled weights, R-side simulation/generation ignored it. This means `predict()`, PPCs, `log_lik()`, WAIC, and LOO were incorrect when using weighted regression. Weighted Gaussian posterior predictions and log-likelihoods now use the observation-level standard deviation `sigma / sqrt(weight)`, matching the JAGS precision `weight / sigma^2`. This makes `predict()`, posterior predictive checks, `log_lik()`, WAIC, and LOO consistent with the fitted model.

-   Bug only noticeable for very small samples: For Gaussian identity-link AR models, R-side calculations could leak residuals between posterior draws and omitted available partial lags for the first observations of AR(2+) models. Each draw is now evaluated as a separate series and uses the same partial-lag recurrence as JAGS.

-   The quantiles for `fitted()` and `predict()` for varying-changepoint models ignored the varying level - they were identical across levels.

-   Now works for 200+ characters formulas too.

-   Fixed #131 (`cores = "all"` failed). Thanks for reporting, @m-r-munroe!

## Behind the scenes

-   Built-in response families now carry their distributional parameters, default priors, response validation, R-side expected-value/log-likelihood/random-generation functions, JAGS likelihood, and optional GARMA behavior in one internal `mcpfamily` specification. The model machinery delegates to that contract instead of branching on family names.

-   Major changes in how the model is translated into JAGS code. The JAGS code is quite different but functionally equivalent.

-   More thorough defensive coding.

-   Much expanded test suite (now 4.000+ tests when run in full).

-   The test suite now includes external validation of inference and simulation: AR against `arima()`/`arima.sim()`, binomial against `glm()` / `rbinom()`.

-   Fewer imports to userspace. This minimizes the risk of name conflicts.

-   Many small improvements in efficiency and code simplicity.

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

-   New argument `scale` in `fitted()`, `plot()`, and `fit$simulate()`. When `scale = "response"` (default), they return fits on the observed scale. When `scale = "linear"`, they return fits on the parameter scale where the linear trends are. Useful for model understanding and debugging.

-   Use `pp_check(fit)` to do prior/posterior predictive checking. See `pp_check(fit, type = "x")` for a list of plot types. `pp_check(fit, facet_by = "varying_column")` facets by a data column.

-   Improvements to `plot()`:

    -   Change point densities are now computed on a per-panel basis in `plot(fit, facet_by = "varying_column")`. Previous releases only displayed population-level change points.
    -   You can now plot varying effects with `rate = FALSE` for binomial models.
    -   Change point densities in `plot(fit)` are not located directly on the x-axis. They were "floating" 5% above the x-axis in the previous releases.

-   New argument `nsamples` reduces the number of samples used in most functions to speed up processing. `nsamples = NULL` uses all samples for maximum accuracy.

-   New argument `arma` in many functions toggles whether autoregressive effects should be modelled.

-   Although the API is still in alpha, feel free to try extracting samples using `mcp:::tidy_samples(fit)`. This is useful for further processing using `tidybayes`, `bayesplot`, etc. and is used extensively internally in `mcp`. One useful feature is computing absolute values for varying change points: `mcp:::tidy_samples(fit, population = FALSE, absolute = TRUE)`. Feedback is appreciated before `tidy_samples` will to become part of the `mcp` API in a future release.

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

-   Model variance for `family = gaussian()` using `~ sigma([formula here])`.

-   Model Nth order autoregressive models using `~ ar(order, formula)`, typically like `y ~ 1 + x + ar(2)` for AR(2). Simulate AR(N) models from scratch or given known data with `fit$simulate()`. The [article on AR(N)](https://lindeloev.github.io/mcp/articles/arma.html) has more details and examples. AR(N) models are popular to detect changes in time-series.

-   Many updates to `plot()`.

    -   Includes the posterior densities of the change point(s). Disable using `plot(fit, cp_dens = FALSE)`.
    -   Supports AR(N) models (see above).
    -   Plot posterior parameter intervals using `plot(fit, q_fit = TRUE)`. `plot(fit, q_fit = c(0.025, 0.5, 0.975))` plots 95% HDI and the median.
    -   Plot prediction intervals using `plot(fit, q_predict = TRUE)`.
    -   Choose data geom. Currently takes "point" (default) and "line" (`plot(fit, geom_data = "line")`). The latter is useful for time series. Disable using `geom_data = FALSE`.

-   Use `options(mc.cores = 3)` for considerable speed gains for the rest of the session. All vignettes/articles have been updated to recommend this as a default, though serial sampling is still the technical default. `mcp(..., cores = 3)` does the same thing on a call-by-ball basis.

-   `fit$simulate()` adds the simulation parameters as an attribute (`attr(y, "simulate")`) to the predicted variable. `summary()` recognizes this and adds the simulated values to the results table (columns `sim` and `match`) so that one can inspect whether the values were recovered.

-   Use `plot(fit, which_y = "sigma")` to plot the residual standard deviation on the y-axis. It works for AR(N) as well, e.g., `which_y = "ar1"`, `which_y = "ar2"`, etc. This is useful to visualize change points in variance and autocorrelation. The vignettes on variance and autocorrelations have been updated with worked examples.

-   Much love for the priors:

    -   Set a Dirichlet prior on the change points using `prior = list(cp_1 = "dirichlet(1)", cp_2 = ...)`. [Read pros and cons here](https://lindeloev.github.io/mcp/articles/priors.html).
    -   The default prior has been changed from "truncated-uniforms" to a "t-tail" prior to be more uninformative while still sampling effectively. [Read more here](https://lindeloev.github.io/mcp/articles/priors.html).
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
