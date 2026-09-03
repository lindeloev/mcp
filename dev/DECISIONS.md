# mcp decisions
This document contains decisions made about mcp, which are knowingly not perfect, but where perfection seem impossible or not worth the added complexity.

See LEARNINGS.md for similar considerations.


## 2026-08-01: Identity links purposefully kept
Identity links can make impossible predictions for Bernoulli, binomial, Poisson, and NB. I am aware and kept them anyway. It is up to the user to make sensible models.

# 2026-08-01: Intercept meaning in each segment
Intercept_1 is at x = 0. Later intercepts are at the onset of the segment, which is random. The Intercept_1 means the parameters of the first segment is compatible with regular regression as we know them from lm(), etc. But brms mean-centers it's intercept, which is nice for interpretability. However, given the random location of the change points, mean-centering within a segment is not trivial - at least I did not find a good solution. Especially given the "handover" nature of mcp-model between segments, I find it defensible to have Intercept_2+ represent the beginning of those intervals. See LEARNINGS.md for another note on mean-centering.

## 2026-08-05: Prediction quantile precision
The `q_predict = TRUE` currently simply runs `quantile(simulated_outcomes, c(0.025, 0.975))`. For small number of draws, it is highly uncertaint. It added ~100 lines of code to improve precision 4-fold per second using rao-blackwell. I opted not to do it. Other solutions are welcome.

## 2026-08-05: Bridge sampling
Requires a lot of work on tracking priors, including their truncation etc. Since bridge sampling is very sensitive to priors anyway, and I doubt users will put a lot of thought into priors, I have opted not to do it.

## 2026-08-08: Check of simulation recovery
If the simulated changepoint location is the same as the location of an actual data point, the changepoint posterior can be confined to just *before* that data point. The current check of recovery is a bit hacky, but everything else explodes in complexity, as far as I can see.

## 2026-08-14: LOO for time series / GARMA
Making `loo` or something loo-like support leave-future-out or blocked sampling would be a large undertaking and it would be slow. One complication is that mcp defaults to data-dependent priors, so changing data changes the priors and hence fit/posterior. It would involve iterations of fitting to before-cutpoint data and doing ELPD on future data. One could ignore and just keep the default priors constant, creating `newdata` that represents leave-future-out and calling `log_lik` on history+newdata. Feasible but slow.

## 2026-08-14: R generics purposefully not implemented
`terms()`, `model.matrix()`, `simulate()`, `update()`. I am not yet ready to commit to a format for these.

coef() will be implemented later. I have not settled on what set of parameters to include.

## 2026-08-14: fit$jags_code and fit$simulate()
Supplying custom `mcp(..., jags_code)` will make fit$simulate() out of sync. This is a known error but requires too much extra tooling to address for now. Options include user-supplied "r_code" to match or a very capable JAGS -> R translator (mcp already does this for it's own models internally).

## 2026-08-17: What is population-level?
`summary.mcpfit()` calls both distributional parameters (`dpar`, e.g., sigma, shape, ...) and AR/MA parameters "population-level". Only change points and group-level effects are not called population-level. This may not be classical wording, but having sections for each `dpar` seems like bloat.

## 2026-08-07: GARMA order > data length
Running ar(4) on length-3 data fails in JAGS with uninformative error. This is such a rare edge case that I will not add ~15 lines to protect against it for each series x group combo.

# 2026-09-03: no novel group predictions
At least for now, fit$simulate() will not support novel categorical levels.