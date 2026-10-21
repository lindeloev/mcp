# mcp decisions
This document contains decisions made about mcp, which are knowingly not perfect, but where perfection seem impossible or not worth the added complexity.

See LEARNINGS.md for similar considerations.


## 2026-08-01: Identity links purposefully kept
Identity links can make impossible predictions for Bernoulli, binomial, Poisson, and NB. I am aware and kept them anyway. It is up to the user to make sensible models.

# 2026-08-01: Intercept meaning in each segment
Intercept_1 is at x = 0. Later intercepts are at the onset of the segment, which is random. The Intercept_1 means the parameters of the first segment is compatible with regular regression as we know them from lm(), etc. But brms mean-centers it's intercept, which is nice for interpretability. However, given the random location of the change points, mean-centering within a segment is not trivial - at least I did not find a good solution. Especially given the "handover" nature of mcp-model between segments, I find it defensible to have Intercept_2+ represent the beginning of those intervals. See LEARNINGS.md for another note on mean-centering.

## 2026-08-05: Bridge sampling
Requires a lot of work on tracking priors, including their truncation etc. Since bridge sampling is very sensitive to priors anyway, and I doubt users will put a lot of thought into priors, I have opted not to do it.

## 2026-08-08: Check of simulation recovery
If the simulated changepoint location is the same as the location of an actual data point, the changepoint posterior can be confined to just *before* that data point. The current check of recovery is a bit hacky, but everything else explodes in complexity, as far as I can see.

## 2026-08-14: LOO for time series / GARMA
Making `loo` or something loo-like support leave-future-out or blocked sampling would be a large undertaking and it would be slow. One complication is that mcp defaults to data-dependent priors, so changing data changes the priors and hence fit/posterior. It would involve iterations of fitting to before-cutpoint data and doing ELPD on future data. One could ignore and just keep the default priors constant, creating `newdata` that represents leave-future-out and calling `log_lik` on history+newdata. Feasible but slow.

## 2026-08-14: R generics purposefully not implemented
`terms()`, `model.matrix()`, `simulate()`, `update()`. I am not yet ready to commit to a format for these. `fit$simulate()` will probably always replace `simulate()`

## 2026-08-14: fit$jags_code and fit$simulate()
 mcp currently simply throw a warning if using fit$simulate() with custom jags_code. Supplying custom `mcp(..., jags_code)` will make fit$simulate() out of sync. This is a known error but requires too much extra tooling to address for now. Options include user-supplied "r_code" to match or a very capable JAGS -> R translator (mcp already does this for it's own models internally).

## 2026-08-17: What is population-level?
`summary.mcpfit()` calls both distributional parameters (`dpar`, e.g., sigma, shape, ...) and AR/MA parameters "population-level". Only change points and group-level effects are not called population-level. This may not be classical wording, but having sections for each `dpar` seems like bloat.

## 2026-08-07: GARMA order > data length
Running ar(4) on length-3 data fails in JAGS with uninformative error. This is such a rare edge case that I will not add ~15 lines to protect against it for each series x group combo.

# 2026-09-03: no novel group predictions
At least for now, fit$simulate() will not support novel categorical levels.

# 2026-09-03: no inits for change points
`mcp` internally assigns good inits for default change points. Adding user-override for these would add more complexity than I like for such an edge case.

# 2026-09-03: keep mcp priors for sigma
The prior for intercept-sigma and for `sigma(1)` (still Intercept, but on log-scale) are quite different. But since these are brms-defaults, this behavior is kept.

# 2026-09-07: There will be no support for 0-parameter models
Fitting a model like `y ~ 0` (only one segment) is not what mcp is for, and it would add a lot of code complexity to support it.

# 2026-09-07: Will not check user priors for change points
Users can specify non-ordered change points (`prior = list(cp_1 = 20, cp_2 = 10)`). mcp should not spend code and complexity checking this and raising informative errors.

# 2026-09-12: residuals() are relative to fitted values - not predicted values
`brms` defaults to `resid = y - predicted` whereas mcp use `resid = y - fitted`, aligning with base R and classical statistics. This is more R-native.

# 2026-09-23: Do not warn on negative support for group-SD
Users overriding default priors can specify priors that result in negative sigma. E.g. `prior = list(Intercept_1_id_sd = "dnorm(0, 1)")`. Checking this complicates too much compared to the relatively infrequent use case and will fail during sampling. For now, we trust users to be informed. Similarly, users specifying negative SD (`dnorm(0, -1)`) will not get an error because mcp squares SD in conversion to JAGS precision.

# 2026-10-10: Residual SD floor of 1e-5 and median inits
The `sigma` prior is truncated at zero (`T(0, )`), but the residual SD in the likelihood is floored at `max(1e-5, sigma)` in JAGS and mirrored R-side in `fit$simulate()`. Lower floors (1e-6 and below) break the JAGS glm module (CHOLMOD "not positive definite") when the residual SD approaches zero.

# 2026-10-10: No init dispersion
JAGS 5.0 inits bounded parameters near the prior mode (zero for half-t SDs). mcp initializes all scalar nodes with numeric priors at their prior medians, keeping equispaced change points. Dispersed inits (prior quantiles per chain) were tried, but started on impossible values in too many cases (extreme scales; identity-link; etc.).

# 2026-10-12: iter excludes warmup
This behavior breaks with brms and some others. But I think it is the correct one so it is purpusefully kept like this.

# 2026-10-12: Formula environments will not be included
Base R allows using variables from the calling environment in formulas. mcp requires formula terms to be present in data. For now, this behavior will be kept. Including the environment in formulas leads to strange situations, e.g., shipping development environments to CRAN; shipping potential sensitive data, etc.

# 2026-10-16: fit$simulate(fit, ...) is preferred over simulate(fit)
Even though it looks odd to write `fit` twice, `fit$simulate()` gives code completion for model parameters. I will likely not support `simulate(fit, ...)` for this reason.

# 2026-10-22: vcov() follows fixef(); confint() covers all population-level parameters
This mirrors `lme4`, where `vcov()` covers the fixed effects and `confint()` also covers variance components and `sigma`. That is a better precedent than base R `lm()` and `glm()`, from which this behavior differ. 

mcp defaults to return summaries for `dpar = "mu"` because the constant `sigma_1` is also a coefficient in mcp. `brms` `vcov()` and `fixef()` include all modeled dpars (e.g., `sigma_Intercept`).