# Changelog

## mcp (development version)

## mcp 0.4.0

CRAN release: 2026-10-26

### Major new features

- **Multiple regression:** `mcp` now supports several continuous
  predictors, categorical predictors, interactions, etc. for all terms
  on RHS. E.g., `~ 1 + x + x:state + sigma(1 + state) + ar(2, 0 + z)`.
  Basically, it now aims to “feel” like
  [`lm()`](https://rdrr.io/r/stats/lm.html) or
  [`glm()`](https://rdrr.io/r/stats/glm.html) for each distributional
  parameter in each segment. Explore `ex = mcp_example("multiple")` to
  see it in action and read more in [Understanding mcp
  formulas](https://lindeloev.github.io/mcp/articles/formulas.html).
  Default priors follow `rstanarm` and `brms` where possible. With this
  possibility for more parameters also comes the convenient `same()`,
  which shares coefficients across segments - including categorical
  predictors and group-level effects ([read
  more](https://lindeloev.github.io/mcp/articles/formulas.html#reusing-coefficients-between-segments)).

- **Group-level effects on RHS:** Predictor formulas now support
  group-level effects (random effects) using familiar `lme4` and `brms`
  syntax. `(1 | group)` specifies a group-level intercept, while
  `(1 + x || group)` and `(factor || group)` support independent
  coefficients, including slopes and factors. This also works inside
  distributional formulas such as `sigma(1 + (factor || id))`. See
  `mcp_example("group_mu")` and the [article on group-level
  effects](https://lindeloev.github.io/mcp/articles/group_effects.html)
  for worked examples. Correlated multi-coefficient terms are not yet
  supported.

- **Native feel:** `mcpfit`s now work natively with R generics,
  [posterior](https://mc-stan.org/posterior/) and
  [tidybayes](https://mjskay.github.io/tidybayes/) posterior draw and
  prediction API. Changes include:

  - `summary(fit)` now reports rank-normalized split-`rhat`, `ess_bulk`,
    and `ess_tail` from [posterior](https://mc-stan.org/posterior/) and
    central quantile intervals instead of HDIs.

  - Adding `as_draws(fit)`, `as_draws_df(fit)`, `tidy_draws(fit)`,
    `ndraws(fit)`, `variables(fit)`, `posterior_epred(fit)`,
    `posterior_predict(fit)`, `posterior_linpred(fit)`,
    [`loo_compare()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
    etc. with full S3 generic registration for
    [posterior](https://mc-stan.org/posterior/),
    [rstantools](https://mc-stan.org/rstantools/), and
    [tidybayes](https://mjskay.github.io/tidybayes/)
    ([`spread_draws()`](https://mjskay.github.io/tidybayes/reference/spread_draws.html),
    [`gather_draws()`](https://mjskay.github.io/tidybayes/reference/spread_draws.html)).
    `bayesplot::mcmc_trace(fit)` etc. work too now.

  - Per-draw methods (`summary = FALSE` in
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md))
    now return dot-prefixed columns (`.epred`, `.prediction`,
    `.residual`, `.loglik`) for
    [tidybayes](https://mjskay.github.io/tidybayes/) /
    [ggdist](https://mjskay.github.io/ggdist/) compatibility.
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    returns `.linpred` with `scale = "linear"` and `[.dpar]` for other
    dpars (e.g., `.sigma`, `.ar1`, or `mu`). See [Fits and
    predictions](https://lindeloev.github.io/mcp/articles/predict.html).

  - `nsamples` soft-deprecated in favor of `ndraws`; `which_y`
    deprecated in favor of `dpar`.

  - Added R generics for mcpfits: `formula(fit)`, `family(fit)`,
    `model.frame(fit)`, `nobs(fit)`, `vcov(fit)`, `confint(fit)`, and
    `coef(fit)`. As in `lme4` and `brms`,
    [`coef()`](https://rdrr.io/r/stats/coef.html) returns absolute
    group-specific coefficients and change points, while
    [`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
    returns relative deviations. Fits now also store a proper matched
    `$call`. Printed summaries now report posterior `sd` and have a new
    layout.

- **Dirichlet is now the default prior on change points:** For
  one-change-point models, this prior is identical to v0.3.4. The
  default change point prior is now the mathematically beautiful and
  principled Dirichlet prior, which was previously only supported as a
  secondary option. It is uniform order statistics, i.e., uniform priors
  on all change points, restricted to observed data range,
  `cp_i = "dunif(min(x), max(x))"`, and ordered: `cp_1 < cp_2 < ...`.
  After ordering, this results in a scaled beta prior for each change
  point. [See a more formal description
  here](https://lindeloev.github.io/mcp/articles/default_priors.html#dirichlet).
  While the dunif()-specification is the new default, the old
  `cp_i = "dirichlet(alpha)"` syntax remains supported (see breaking
  changes below). New: set a minimum segment length, e.g.,
  `cp_2 = "dunif(cp_1 + 20, max(x))"`. For models with 2+ change points,
  the previous sequentially truncated t-distributions centered at
  `min(x)` have been removed. Truncated-t only served the purpose of
  higher ESS/s at the expense of being mathematically unprincipled. A
  new implementation in JAGS makes it fast. Although this is a
  mathematical change, parameter posteriors remain essentially identical
  to previous fits for all
  [`mcp_example()`](https://lindeloev.github.io/mcp/dev/reference/mcp_example.md).
  Any `hypothesis(fit, "cp_1 = 20")` (Savage-Dickey Bayes factor) on
  default priors will likely differ, but BF is ill advised on default
  priors anyway - that is why
  [`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md)
  now throws a warning if this is attempted.

- **Negative binomial and offset():** You can now do
  `mcp(..., family = negbinomial())`. NB (and Poisson) can now model
  rates and exposure because
  [`offset()`](https://rdrr.io/r/stats/offset.html) is now supported in
  formulas - also for other families. If one segment includes an
  [`offset()`](https://rdrr.io/r/stats/offset.html), all segments must
  include one. See [the Poisson and negative-binomial
  article](https://lindeloev.github.io/mcp/articles/poisson.html#modeling-rates-using-offset)
  for an applied example.

- **Missing response imputation:** Missing responses are now retained as
  posterior JAGS imputations and used to complete response histories in
  AR/MA models. Predictions generate fresh replicated outcomes at every
  row; `predict(fit) |> filter(is.na(y))` can be used to see predictive
  intervals and draws at missing points, and similarly
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  for expected responses. Missing AR/MA histories remain unsupported in
  [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  [`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
  and
  [`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md).
  See [the article on missing
  data](https://lindeloev.github.io/mcp/articles/missing.html).

- **Improved speed and robustness:** ESS/s has improved more than
  10-fold and convergence is more likely. This is due to
  [reparameterizations](https://lindeloev.github.io/mcp/articles/implementation.html)
  that keep intercepts and slopes relatively independent of earlier
  segments during sampling and due to (selectively) using the JAGS glm
  module. The reparameterizations do not change the model or its priors.

### Major breaking changes

mcp v0.4 is a major breaking change with the aim of remaining relatively
stable going forward towards version 1.0. Although a lot has been
updated, the parameter estimates in v0.4.0 remain practically identical
to the previous public release (v0.3.4). Deprecation detections have
been added until we reach 1.0.

- **Renamed parameters and column names** to align with
  [brms](https://github.com/paul-buerkner/brms): `int_i` –\>
  `Intercept_i`; `x_1_E2` –\> `xE2_1`; `x_1_sin` –\> `sinx_1`, `cp_1_sd`
  –\> `cp_1_id_sd` (for `(1 | id)`), etc. Interactions join their parts
  with a dot, e.g., `x:stateB` –\> `x.stateB_1`.

- **Renamed return column names** to align with
  [posterior](https://mc-stan.org/posterior/) and others:

  - [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`fixef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    and
    [`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
    now name the parameter identifier column `parameter` instead of
    `name`, matching
    [`prior_summary()`](https://lindeloev.github.io/mcp/dev/reference/prior_summary.md)
    and
    [`mcp_pars()`](https://lindeloev.github.io/mcp/dev/reference/mcp_pars.md).

  - [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    with `summary = TRUE` now take `width = 0.95` instead of `probs`
    (deprecated with warning) and return the columns `.mean`, `.sd`,
    `.lower`, and `.upper` instead of
    `fitted`/`predict`/`residuals`/`loglik`, `error`, `Q2.5`, and
    `Q97.5`. This matches
    [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
    (`width`, `mean`, `sd`, `lower`, `upper`), while the dot prefix
    prevents collisions with data columns, as in `broom::augment()`.
    Consequently, data columns named `sd`, `fitted`, `predict`, etc. are
    no longer rejected. For other quantiles, summarise the draws from
    `summary = FALSE`. The evaluation-row key `data_row` is now `.row`,
    matching [tidybayes](https://mjskay.github.io/tidybayes/).

  - `plot(fit, q_fit, q_predict)` is deprecated in favor of
    `plot(fit, interval = c("fitted", "predict"), width = 0.95)`. This
    matches the new `width` in
    [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).
    Add other quantiles as a ggplot2 layer from
    `fitted(..., summary = FALSE)` instead.

  - `Rhat` was renamed to `rhat` in
    [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`fixef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    and
    [`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    to match [posterior](https://mc-stan.org/posterior/) standard.

  - [`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md)
    returns column `prob` instead of `p` for the posterior (or prior)
    probability of directional hypotheses, preventing confusion with
    frequentist p-values. It also returns the posterior `sd` of the
    hypothesis, like
    [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md).

- **Change points are modeled using `cp()`:** Segments 2+ are now
  one-sided formulas, and change points are modeled using `cp()`, e.g.,
  `list(y ~ 1, ~ cp((1 | id)) + 0 + x)` instead of
  `list(y ~ 1, 1 + (1 | id) ~ 0 + x)`. This avoids confusing the
  left-hand side with the response variable. The old syntax is
  deprecated and works with a warning that shows the new formula.
  `formula(fit)` and `summary(fit)` show the new syntax.

- **Parameters and outputs:**

  - `fit$pars` and `fit$mcmc_loglik` were removed and now raise helpful
    errors. `fit$pars` is now replaced with `mcp_pars(fit)` and
    `mcp_columns(fit)`. `fit$mcmc_loglik` is no longer stored - generate
    using `log_lik(fit)`.

  - `mcp:::tidy_samples()` is removed. Use `as_draws_df(fit)` or
    [`tidybayes::spread_draws()`](https://mjskay.github.io/tidybayes/reference/spread_draws.html).

  - [`fixef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
    now reports only population-level effects for `mu` (the primary
    response parameter), excluding change points, other distributional
    parameters, AR/MA parameters, and group-effect SDs. Use
    `summary(fit)` to get all parameters.

  - In general, renamed `"varying"` to `"group"`, e.g., it is now
    `plot_pars(fit, pars = "group")` and
    `fitted(fit, group_effects = "x_2_id")`. `varying` is a deprecated
    alias with a warning.

  - [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`fixef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`prior_summary()`](https://lindeloev.github.io/mcp/dev/reference/prior_summary.md),
    [`as_draws()`](https://lindeloev.github.io/mcp/dev/reference/as_draws.mcpfit.md),
    `fit$prior`, `fit$simulate()`, and everything else now order
    parameters canonically as in
    [`mcp_pars()`](https://lindeloev.github.io/mcp/dev/reference/mcp_pars.md)
    instead of the previous incidental (near-alphabetical) order; group
    levels follow the grouping factor’s levels in the fitting data.
    Draws from
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    etc. with `summary = FALSE` are ordered by data row, then draw, as
    in [tidybayes](https://mjskay.github.io/tidybayes/). Use
    `verbose = TRUE` with
    [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`fixef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    or
    [`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
    to include `segment` and `dpar` columns.

  - Tibbles are generally returned instead of data.frame() for
    readability. Tibbles are data frames, but `x[, "col"]` now returns a
    one-column tibble instead of a vector. Use `x$col` or `x[["col"]]`
    to get a vector. This change involves
    [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`fixef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md),
    and
    [`interpolate_newdata()`](https://lindeloev.github.io/mcp/dev/reference/interpolate_newdata.md).

  - mcp no longer exports `logit`, `ilogit`, or `probit`.
    [`phi()`](https://lindeloev.github.io/mcp/dev/reference/mcp-deprecated.md)
    is deprecated in favor of
    [`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html).

- **Functions, plotting, and workflow:**

  - [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
    is now split into
    [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
    for plotting full fits while
    [`plot_dpar()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
    plots one distributional parameter (`mu`, `sigma`, `shape`, `ar1`,
    etc.). The argument order was changed too. The new coloring function
    (`plot(fit, color_by = "column")`) is particularly useful when
    models include categorical predictors or rhs group-level effects.
    The default `color_by = "auto"` colors by a sole categorical
    predictor. See `mcp_example("group_mu")` and the [plotting
    article](https://lindeloev.github.io/mcp/articles/plotting.html) for
    worked examples.

  - Parallel sampling is now controlled exclusively through the active
    [future](https://future.futureverse.org) plan. The `cores` argument
    to [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md)
    is deprecated and ignored, but remains available for backwards
    compatibility. Use `future::plan(future::multisession, workers = 3)`
    before calling
    [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) to
    sample chains in parallel, and `future::plan(future::sequential)` to
    shut down those workers. Without a parallel future plan, chains are
    sampled sequentially. Chain-specific initial values
    (`inits = list(list(...), list(...))`) are now supported under
    parallel execution and when supplying an explicit `seed`.

  - **Binomial predictions now default to counts:**
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    returns expected successes and
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    returns predicted successes (`rate = FALSE`), matching
    [`posterior_epred()`](https://lindeloev.github.io/mcp/dev/reference/posterior_epred.mcpfit.md),
    [`posterior_predict()`](https://lindeloev.github.io/mcp/dev/reference/posterior_epred.mcpfit.md),
    and brms. In v0.3.4 these methods defaulted to proportions. Update
    old scripts with explicit `rate = TRUE` to retain proportions, or
    `rate = FALSE` to select counts. An omitted `rate` warns once per
    session per function when the requested counts differ from
    proportions.

  - **ar(p) no longer continues into later segments:** Only terms
    included in a segment get coefficients there. While important for
    the new multiple-regression functionality, the only backward
    incompatibility is that in a model like
    `model = list(y ~ 1 + ar(1), ~ 0)`, the AR would continue into
    segment 2. The equivalent model in v0.4+ is `~ 0 + ar(1, 0)` in
    segment 2.

  - Explicit prior draw selection: Methods requiring draws now require
    posterior draws by default and error if they are unavailable. Use
    `prior = TRUE` to select prior draws; prior-only fits no longer
    trigger an automatic fallback. This also applies to draw-count
    methods, and `print(fit)` shows only the model unless
    `prior = TRUE`.

  - The arguments for `fit$simulate()` have all changed to accommodate
    multiple (categorical) predictors.
    `fit$simulate(fit, newdata, ..., .type = "predict")` is the new
    argument structure. Note that (1) it now requires `fit` as the first
    argument, (2) it requires `data.frame` or `tibble` as the second
    argument instead of just a vector of `par_x`, (3) further arguments
    are prefixed with a “.” to avoid name conflicts internally in mcp.
    `...` are the model parameters as usual. `fit$simulate()` now also
    preserve shared group effects specified in priors (e.g.,
    `prior = list(Intercept_2_id = "Intercept_1_id")`).

  - [`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md)
    takes the fit as its first argument `x` (was `fit`), matching
    [`brms::hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html).
    It now restricts equality tests to simple parameter contrasts and
    gives clearer guidance on Savage-Dickey Bayes factors. Non-linear
    transformations in hypotheses like `x_1 / x_2 = 1` could previously
    be run, even though Savage-Dickey does not support it. Equality
    tests now return `prob = NA` because their Bayes factor is a model
    comparison, not a posterior parameter probability. Density estimates
    are more robust, warn when evaluating Savage-Dickey Bayes factors on
    default priors, and warn about sparse tails. A Savage-Dickey `BF` is
    `0` when no posterior draws lie on one side of the tested value, as
    for directional hypotheses. See [Savage-Dickey Bayes
    factors](https://lindeloev.github.io/mcp/articles/comparison.html#bayes-factors-using-savage-dickey-density-ratios).

  - `mcp(..., adapt = 1000)` is soft-deprecated in favor of
    `mcp(..., warmup = 1000)` for compatibility with additional future
    samplers.

- **Formula and modeling semantics:**

  - Dropped support for
    [`rel()`](https://ggplot2.tidyverse.org/reference/element.html) in
    formulas. This was ambiguous for interaction terms and made the code
    hard to maintain. Compute relative changes from the posterior draws
    as shown in [Relative changes between
    segments](https://lindeloev.github.io/mcp/articles/formulas.html#relative-changes-between-segments).

  - Dropped support for `mcp(..., data = NULL)`. You now must provide
    some mock-up data to inform `mcp` about the types and levels of the
    predictor columns. See, e.g.,
    `mcp_example("intercepts")$example_code` for a simple example or
    `mcp_example("multiple")$example_code` for a more involved example.
    All docs have been updated appropriately.

  - Gaussian models with an explicit
    [`sigma()`](https://rdrr.io/r/stats/sigma.html) formula now model
    the residual standard deviation with a log link, as in `brms`. To
    make the scale visible, its intercepts are named `sigma_Intercept_i`
    (previously `sigma_i`), matching `Intercept_i`, while slopes keep
    their names (e.g., `sigma_x_2`) and are now on the log scale. Models
    without an explicit [`sigma()`](https://rdrr.io/r/stats/sigma.html)
    are unchanged: `sigma_1` remains the residual SD itself with the
    response-scale half-Student-t prior. See [distributional
    parameters](https://lindeloev.github.io/mcp/articles/dpar.html#priors-for-distributional-parameters).

  - `y | weights(w)` now specifies observation log-likelihood weights
    rather than Gaussian-only precision weights (which previously scaled
    the residual SD as `sigma / sqrt(w)`). Unlike `brms`, weights must
    be strictly positive due to JAGS requirements (cannot be zero).
    Predictive draws and expectations now use the distribution
    parameters directly while
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    multiplies observation log-densities by `w`.

  - Formula transformations now use the original predictor values, as in
    [`lm()`](https://rdrr.io/r/stats/lm.html),
    [`glm()`](https://rdrr.io/r/stats/glm.html), and `brms`. Previously,
    transformations of the change-point predictor such as `sin(x)` and
    `exp(x)` restarted at each change point. Bare `par_x` and polynomial
    bases such as `I(par_x^2)` are still measured from the change point
    to support joined segments.

- **Fitting and priors:**

  - The prior for user-specified priors on 2+ change points will differ
    from mcp v0.3.4, though for informative data, posteriors will be
    highly similar. User-specified priors on change points are now used
    unchanged, though restricted to ordered change points. v0.3.4
    imposed ordering by truncating priors using
    [`T()`](https://rdrr.io/r/base/logical.html) at the previous change
    point, e.g., `dnorm(40, 10) T(cp_1, max(x))`, which pushed wide
    priors towards `max(x)`.
    [`dunif()`](https://rdrr.io/r/stats/Uniform.html) priors and priors
    with [`T()`](https://rdrr.io/r/base/logical.html) were used
    unchanged, which allowed label switching. Fixed change points, e.g.,
    `cp_2 = 80`, are now ordered relative to the other change points
    too. Priors without
    [`dunif()`](https://rdrr.io/r/stats/Uniform.html) or
    [`T()`](https://rdrr.io/r/base/logical.html) are truncated to the
    observed range of `x`. See [priors on change
    points](https://lindeloev.github.io/mcp/articles/priors.html#cp_prior).

  - The lower and upper bounds of change point priors, i.e., the
    arguments of [`dunif()`](https://rdrr.io/r/stats/Uniform.html) and
    [`T()`](https://rdrr.io/r/base/logical.html), can no longer refer to
    other change points, except to set a minimum segment length, e.g.,
    `cp_2 = "dunif(cp_1 + 20, max(x))"`. Change points in the
    distribution’s parameters remain allowed, e.g.,
    `cp_3 = "dnorm(cp_2 + 30, 10)"`.

  - `dirichlet(alpha)` now takes one shared `alpha >= 1` (non-integer
    allowed) for all change points and applies it to all segments,
    including the last. v0.3.4 allowed non-identical integer `alpha` per
    change point and silently used `alpha = 1` for the last segment, so
    `dirichlet(2)` now gives a different prior. Non-identical alphas now
    error.

  - Group-level change points `(1|id)` now use a standard hierarchical
    normal model sampled efficiently through their absolute locations.
    This makes `cp_i_sd` the SD of the normal distribution before
    truncation rather than the realized sample SD.

  - The default prior on group-level change-point SDs (`cp_i_id_sd`) is
    now `dnorm(0, (max(x) - min(x)) / (2 * (n_cp() + 1))) T(0, )`, i.e.,
    a quarter of the expected interval in which each group’s change
    point can vary. Under the old prior, group-level change points were
    barely pooled, and `cp_i` and `cp_i_id_sd` were weakly identified.
    The new prior is approximately uniform on the *realized* SD of the
    group-level change points, which differs from `cp_i_id_sd`, the SD
    *before truncation*. See [default
    priors](https://lindeloev.github.io/mcp/articles/default_priors.html#group-level-change-points)
    for the justification.

  - AR and MA intercepts now have zero-centered, regularizing
    `dnorm(0, 0.5) T(-1, 1)` priors, replacing independent uniform
    priors. Their categorical contrasts and numeric slopes now use
    modest normal priors instead of heavy-tailed Student-t priors.
    Coefficients remain direct and are not jointly constrained to
    stationary or invertible regions ([read
    more](https://lindeloev.github.io/mcp/articles/arma.html#ar_prior)).

  - Corrected the JAGS translation of prior scales. `ddexp()` and
    [`dlogis()`](https://rdrr.io/r/stats/Logistic.html) now convert
    their conventional scale to inverse scale rather than inverse
    variance; The exported
    [`sd_to_prec()`](https://lindeloev.github.io/mcp/dev/reference/sd_to_prec.md)
    helper is soft-deprecated because prior translation is now an
    internal, sampler-specific step.

### Other new features

- **Generalized AR/MA:**

  - Autoregression ([`ar()`](https://rdrr.io/r/stats/ar.html)) has been
    generalized to link-scale observation-driven GARMA residuals for
    Gaussian, binomial, Bernoulli, Poisson, and negative-binomial models
    with their default links, using `ar(..., threshold = 0.1)` by
    default to keep observations finite on the link scale. Added
    moving-average terms with `ma(q)`, which can be used alone or
    combined with `ar(p)` in each segment. Independent time series can
    be separated using `ar(..., series = <column>)`. Set `series` and
    `threshold` once in any [`ar()`](https://rdrr.io/r/stats/ar.html) or
    `ma()` term; they apply to all segments. See the [article on time
    series](https://lindeloev.github.io/mcp/articles/arma.html).

  - GARMA observation thresholds are now validated against the support
    of the response: strictly below 0.5 for `family = bernoulli()`, and
    strictly below half the smallest trial count for
    `family = binomial()`, preventing thresholds from erasing the binary
    response history.

  - Explicit AR/MA prediction history control:
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    and the new
    [`posterior_predict()`](https://lindeloev.github.io/mcp/dev/reference/posterior_epred.mcpfit.md)
    support `arma = "replicate"` to generate fresh response series
    recursively, in addition to `arma = TRUE` (default) to condition on
    observed response histories and `arma = FALSE` to omit AR/MA
    effects. This applies to both prior and posterior prediction.
    Previously, AR history handling depended on whether response values
    were supplied.

  - Added AR/MA warnings: (1) for AR/MA models to
    [`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    etc. where the serial dependence is currently ignored. Proper
    handling requires leave-future-out or blocked cross-validation,
    which are not currently implemented in `mcp`. (2) when posterior
    violation probabilities exceed the configurable `ar` or `ma`
    diagnostic thresholds (10% by default).

- **New features in
  [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md):**

  - Added warnings when various thresholds are exceeded.
    `mcp(..., diagnostics = list(rhat = 1.01, ess_bulk = 400, ess_tail = 400, ar = 0.10, ma = 0.10))`
    controls fit diagnostics with configurable thresholds. Partial lists
    override these defaults, individual diagnostics can be disabled with
    `NULL`, and `diagnostics = FALSE` disables diagnostic warnings.
    `summary(fit)` inherits these settings and accepts its own override
    for the diagnostic footer.

  - Use `mcp(..., seed = 42)` or `set.seed(42)` before
    [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) for
    reproducible JAGS sampling, whether chains run sequentially or in
    parallel. See `mcp_example("demo")$example_code` how to ensure
    reproducibility across simulation, fit, and plotting.

  - `mcp(..., quiet = TRUE)` suppresses routine JAGS output and mcp
    sampling-status messages while preserving warnings and errors. Use
    `options(mcp_quiet = TRUE)` to make this the default for the
    session. The JAGS model-compilation text is no longer printed by
    default.

  - `fit$version` records the versions of mcp, R, rjags, and JAGS used
    to make the fit, e.g., `fit$version$mcp >= "0.4.0"`.

  - Memory improvement: The `mcpfit` is now \< 10% of the size as before
    because the log-likelihood is not stored. Use `log_lik(fit)` to
    compute it, or call `loo(fit)` or `waic(fit)` directly.

- **Model evaluation and prediction:**

  - Added `log_lik(fit)` which returns a draws-by-observation matrix by
    default. This is the usual
    [brms](https://github.com/paul-buerkner/brms) return shape and is
    accepted directly by [loo](https://mc-stan.org/loo/).
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    supports the same arguments as
    e.g. [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).

  - Methods like
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    and
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    now select group-level effects using `group_effects` (previously
    `varying`). `group_effects = FALSE` or `NA` gives population-level
    predictions, `"cp"` and `"predictor"` select all group-level effects
    in that formula part, and parameter names accept `*` wildcards,
    e.g., `fitted(fit, group_effects = c("cp_*", "state*_1_id"))`.
    [`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
    continues to return all group-level effects.

  - Fits made with custom `jags_code` now warn when calling R-side
    simulation, prediction, or model-evaluation methods, because those
    methods continue to evaluate the supplied formulas rather than the
    custom JAGS code. `jags_code` is planned for deprecation once a
    future release brings custom families and tooling to run `jags_code`
    directly in JAGS. If you need this, please [open an
    issue](https://github.com/lindeloev/mcp/issues) describing your use
    case.

  - [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
    etc. are now more memory efficient (e.g., ~2.5 GB instead of ~12 GB
    for 2000 observations). Use `loo(fit, ndraws = 1000)` to compute LOO
    on a subset of draws.

- **Priors and hypotheses:**

  - Added `prior_summary(fit)`. Its compact output shows each
    parameter’s resolved prior and its lower and upper bounds, marking
    change points as ordered; `prior_summary(fit, verbose = TRUE)` also
    shows the data-dependent rule, a plain-language description, source,
    and kind (`distribution`, `alias`, `expression`, or `constant`). See
    [Working with
    priors](https://lindeloev.github.io/mcp/articles/priors.html#inspect-resolved-priors-and-their-bounds).

  - With the introduction of multiple regression including categorical
    predictors (including for distributional parameters) come default
    priors. Coefficients are autoscaled as in
    [`rstanarm`](https://mc-stan.org/rstanarm/reference/priors.html),
    i.e., divided by the SD of their model-matrix column (see “Priors
    and distributions” below). See `prior_summary(fit, verbose = TRUE)`
    for details.

  - Added option to test hypotheses on the prior using
    `hypothesis(fit, prior = TRUE)`.

  - `y | weights(w_col)` now supported for all families
    ([`gaussian()`](https://rdrr.io/r/stats/family.html),
    [`binomial()`](https://rdrr.io/r/stats/family.html),
    [`bernoulli()`](https://lindeloev.github.io/mcp/dev/reference/bernoulli.md),
    [`poisson()`](https://rdrr.io/r/stats/family.html),
    [`negbinomial()`](https://lindeloev.github.io/mcp/dev/reference/negbinomial.md)).
    Applies log-likelihood weights as brms.

- **Plotting and examples:**

  - Default `plot(fit)` style has been updated in many ways to
    accommodate multiple regression and group-effects.

  - Predictive summaries are now exact and faster. The intervals
    (`plot(fit, interval = "predict")` and
    `predict(fit)$.lower`/`$.upper`) and moments
    (`predict(fit)$.mean`/`$.sd`) are computed from the posterior
    predictive distribution instead of from random predictions.

  - Added `interpolate_newdata(fit, by = NULL)` which generates a
    data.frame with all combinations of categorical predictors along
    with interpolated continuous predictors. Use `by` to include
    grouping factors. The documentation shows how this can be useful for
    generating [custom
    plots](https://lindeloev.github.io/mcp/articles/plotting.html#making-your-own-plots)
    when simple tweaking
    [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
    is not enough.

  - `mcp_example(...)` now shows an illustrative plot of the fitted
    example by default. Disable using `plot = FALSE`.

  - New help topics `?mcp-formula` (model formulas and the underlying
    model) and `?mcp-priors` (prior syntax), split out from
    [`?mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md).

### Minor breaking changes

- **Priors and distributions:**

  - Default coefficient priors now follow the autoscaled defaults of
    `rstanarm`: `normal(0, 2.5 * sd(y) / sd(x))` for identity links and
    `normal(0, 2.5 / sd(x))` for log links, where `sd(x)` is the
    standard deviation of the predictor’s model-matrix column, including
    dummy variables for factors. For the change-point variable, `sd(x)`
    is taken over all data (so slope priors do not depend on the number
    of change points), powers such as `I(x^2)` use `sd((x - min(x))^2)`,
    and interactions such as `x:z` use `sd(x) * sd(z)` (or just `sd(x)`
    for 0/1 dummies as in `x:state`). Version 0.3.4 used
    `dt(0, sd(y) / (max(x) - min(x)), 3)` for Gaussian slopes, which is
    about 9 times narrower for evenly spread `x` and noticeably shrank
    steep slopes in short segments. Intercepts also follow `rstanarm`
    (`normal(mean(y), 2.5 * sd(y))` for Gaussian models), and the
    default prior on the segment-1 intercept applies to the level at
    `min(x)`, just as `rstanarm` and `brms` place intercept priors on
    centered predictors. `Intercept_1` is still reported at `x = 0`, as
    in [`lm()`](https://rdrr.io/r/stats/lm.html), so `fit$prior` and
    [`prior_summary()`](https://lindeloev.github.io/mcp/dev/reference/prior_summary.md)
    show its prior given the x-terms, e.g.,
    `dnorm(mean(y) - x_1 * min(x), 2.5 * sd(y))`. This greatly improves
    sampling when `x` is far from zero, e.g., calendar years. Likewise
    for factor levels and contrasts with their own x-terms, e.g.,
    `stateb_1` in `y ~ 1 + state * x`. `sigma` and group-level SDs
    follow `brms`. See priors using `prior_summary(fit, verbose = TRUE)`
    and read more about [default
    priors](https://lindeloev.github.io/mcp/articles/default_priors.html).
    The implicit Gaussian `sigma_1` prior is now the same
    response-calibrated half-Student-t used by `brms`,
    `dt(0, max(2.5, round(mad(y), 1)), 3) T(0, )`. Version 0.3.4 instead
    used a response-SD-calibrated half-normal prior. Both are weakly
    informative so will have negligible impact on fits.

  - Uppercase data-dependent constants for priors (`MINX`, `MAXX`, etc.)
    remain temporarily supported with a deprecation warning. They are
    now replaced with readable expressions such as `min(time)`,
    `max(time) - min(time)`, `mad(response)`, `n_segments()`, and
    `n_cp()`.

  - Priors can no longer refer to the internal JAGS nodes `cp_0` and
    `cp_K` (e.g., `cp_3` with two change points). Use `min(x)` and
    `max(x)` instead, as also shown in `fit$prior`.

  - Binomial and Bernoulli models with logit or probit links now use
    narrower, weakly regularizing scales than `rstanarm`’s 2.5:
    `normal(0, 1.5)` (logit) and `normal(0, 1)` (probit) for intercepts,
    which are approximately and exactly uniform on the probability
    scale, respectively. Coefficients and categorical contrasts use
    `normal(0, 1.5 / sd(x))` and `normal(0, 1 / sd(x))`. See [the
    justification](https://lindeloev.github.io/mcp/articles/default_priors.html#deviations).

  - **Gaussian models with log link (`gaussian(link = "log")`):** Now
    aligned with `mcp`’s general log-link model standards (`dnorm`
    rather than Student-t to reduce extreme multiplicative tails under
    exponentiation). Population intercepts are robustly calibrated to
    empirical log-data (with non-positive values replaced by 0.1, as
    `brms` does for zeros) and shift with
    [`offset()`](https://rdrr.io/r/stats/offset.html), and coefficients
    and categorical contrasts use `rstanarm`’s `dnorm(0, 2.5 / sd(x))`.

  - **Count model priors
    ([`poisson()`](https://rdrr.io/r/stats/family.html),
    [`negbinomial()`](https://lindeloev.github.io/mcp/dev/reference/negbinomial.md)):**
    In v0.3.4, log-link Poisson models used uncalibrated, unscaled
    `dnorm(0, 10)` for both intercepts and slopes. In v0.4.0, population
    log-mean intercepts are now robustly data-calibrated to
    `dnorm(median(log(y)), ...)` (zeros replaced by 0.1, as in `brms`)
    with a minimum scale of 2.5, and automatically shift to the observed
    log-rate (\log(y) - \text{offset}) when an
    [`offset()`](https://rdrr.io/r/stats/offset.html) is present
    (calibrating location and scale directly to rates, unlike `brms`,
    which subtracts `mean(offset)` from the location but computes the
    scale from the log-counts). Coefficients and categorical contrasts
    now use `rstanarm`’s `dnorm(0, 2.5 / sd(x))`. For the new RHS
    group-level features in v0.4.0, group SDs default to
    `dnorm(0, 2.5) T(0, )` (or scaled half-normals). Broad normal priors
    follow the precedent of `rstanarm`, substantially reducing extreme
    multiplicative effects and explosive prior simulations compared with
    Student-t priors under exponentiation. Similarly, modeled log-shape
    in negative-binomial (`shape(...)`) uses `dnorm(0, 2.5)` and
    `dnorm(0, 2.5 / sd(x))` rather than Student-t priors ([read
    more](https://lindeloev.github.io/mcp/articles/poisson.html#priors-for-poisson-models)).

- **Interface and data structures:**

  - `fit = mcp_example("name")` now returns the fit directly instead of
    a list with a `$fit` entry. It now defaults to sampling the model
    (`sample = "post"`) and the `sample` argument is now directly passed
    to `mcp(..., sample = sample)` so `sample = TRUE` is deprecated.

  - The term “ct” (for “central tendency”) has been replaced with
    `dpar = "epred"` (the default), e.g., in
    `plot_dpar(fit, dpar = "epred")` and `fitted(fit, dpar = "epred")`.
    Deprecated `which_y = "ct"` is translated to `dpar = "epred"` with a
    warning.

  - `mcp(..., inits)` can only set parameters in `mcp_pars(fit)` and
    errors on other names. Parameters reparameterized for faster
    sampling can no longer be initialized. Initial values for change
    points must be ordered ([read
    more](https://lindeloev.github.io/mcp/articles/implementation.html#initial-values)).

  - `fit$data` now only contains the data columns that are used in the
    model.

  - `draws_format` replaces `samples_format` in
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).
    The old argument remains available with a soft deprecation.
    [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    now also takes `draws_format`, e.g.,
    `residuals(fit, summary = FALSE, draws_format = "matrix")` for a
    draws-by-observation matrix like `brms`.

  - Removed `which_y` argument from
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).

  - `rate = TRUE` in
    [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    and
    [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
    now errors for families other than binomial and Bernoulli. In
    v0.3.4, it was silently ignored. For rates per unit of exposure in
    models with `offset(log(exposure))`, set exposure to 1, e.g.,
    `fitted(fit, newdata = transform(df, exposure = 1))`.

  - [`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
    and
    [`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
    always evaluate the fitted model, including all group-level effects
    and AR/MA terms. Use `log_lik(fit, group_effects = FALSE)` or
    `log_lik(fit, arma = FALSE)` for conditional log-likelihoods.
    Additionally, `loo(ndraws = ...)` and `waic(ndraws = ...)` now thin
    draws evenly across chains to preserve MCMC chain identities and
    order.

  - New `mcpfit` objects no longer include empty `$loo` and `$waic`
    components.

  - [`is.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/is_mcpfit.md)
    is deprecated in favor of
    [`is_mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/is_mcpfit.md),
    matching [posterior](https://mc-stan.org/posterior/)’s `is_draws()`.

- **Requirements and removals:**

  - mcp now requires R \>= 4.1.0, matching its use of the native pipe
    (`|>`).

  - Removed
    [`criterion()`](https://lindeloev.github.io/mcp/dev/reference/mcp-defunct.md)
    and
    [`get_segment_table()`](https://lindeloev.github.io/mcp/dev/reference/mcp-defunct.md).
    They now throw an error pointing to their replacements: `loo(fit)`
    or `waic(fit)` instead of `criterion(fit)`, and `mcp_pars(fit)` or
    `mcp_columns(fit)` on `fit = mcp(..., sample = FALSE)` instead of
    [`get_segment_table()`](https://lindeloev.github.io/mcp/dev/reference/mcp-defunct.md).

  - Removed the non-functional `exponential()` family. It had no default
    priors and could not be fitted. Exponential survival models would
    require dedicated support for censoring and survival likelihoods
    rather than the previous ordinary-response placeholder.

### Bug fixes

- Equality hypotheses now warn when all prior draws lie on one side of
  the tested value: kernel smoothing can otherwise give misleading Bayes
  factors for values outside the prior bounds. Boundary tests receive
  the same warning; Bayes factor calculations are unchanged.

- Named arguments in prior distributions or truncations
  (e.g. `dnorm(sd = 1, mean = 10)`) are now rejected explicitly with an
  informative error rather than being silently re-interpreted
  positionally.

- `hypothesis(..., prior = TRUE)` now reports `BF = NA` rather than
  comparing prior draws with themselves and returning a Bayes factor of
  1.

- Directional Bayes factors now divide posterior odds by prior odds.
  Previously, they reported posterior odds alone, which is only a Bayes
  factor when the prior probability of the hypothesis is 0.5.

- [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  now uses the conventional `observed - fitted` sign. Previously, `mcp`
  residuals used `fitted - observed`.

- Poisson and negative binomial predictions no longer fail when a single
  draw implies a count mean above 2,146,275,819. This was common for
  prior predictions with default log-link priors.

- Fixed several posterior predictive check and information-criterion
  bugs present in v0.3.4. From most to least serious:

  - For missing data: Missing responses were scored in WAIC/LOO as if
    their JAGS-imputed values had been observed. Fix: missing responses
    remain latent in JAGS but are excluded from observed-data PPC,
    log-likelihood, WAIC, and LOO calculations.

  - For models with group-level change points, LOO checks with
    `facet_by != NULL` could pair group-specific predictions with
    weights from the wrong observations.

  - For AR/MA models,
    [`pp_check()`](https://lindeloev.github.io/mcp/dev/reference/pp_check.md)
    conditioned every prediction on the observed response history. It
    now defaults to `arma = "replicate"`, generating each posterior
    replication recursively as a fresh time series. Use `arma = TRUE`
    for the previous one-step-ahead check. LOO predictive checks keep
    `arma = TRUE` to match the one-step-ahead log-likelihood behind
    their weights.

  - LOO checks with `prior = TRUE` incorrectly combined prior
    predictions with posterior LOO weights. A rare use case.

  - `pp_check(..., ndraws = NULL)` errored. Fixed.

- Weighted regression now uses brms-style likelihood weights throughout:
  JAGS targets the likelihood raised to `weight`,
  [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  multiplies each log density by `weight`, and predictive draws use the
  unweighted response distribution. In v0.3.4, weights instead acted as
  Gaussian precision weights and R-side predictions and information
  criteria were not fully consistent with the fitted model.

- In models going from higher-order to lower-order (`~ ar(2), ~ ar(1)`),
  the higher-order components were not turned off.

- Change points without a population-level intercept (e.g.,
  `0 + (1 | id) ~ 1`) now raise an error. Previously, `0 +` was silently
  ignored. Group-level change points now also accept `(1 || id)`, which
  is identical to `(1 | id)` for intercepts.

- Bug only noticeable for very small samples: For Gaussian identity-link
  AR models, R-side calculations could leak residuals between posterior
  draws and omitted available partial lags for the first observations of
  AR(2+) models. Each draw is now evaluated as a separate series and
  uses the same partial-lag recurrence as JAGS.

- Conditional R-side GARMA histories are now evaluated across posterior
  draws rather than through a scalar loop over every draw and
  observation. This substantially speeds up
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  [`residuals()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  and
  [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  for AR/MA models.

- The quantiles from
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  and
  [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  for group-level change-point models ignored the group level and were
  identical across levels.

- Explicitly setting `dpar = "mu"` now returns the success probability
  parameter \mu regardless of `rate`. `rate` only affects `"epred"` and
  predictions, where the default `rate = FALSE` returns counts (e.g., in
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  and
  [`posterior_epred()`](https://lindeloev.github.io/mcp/dev/reference/posterior_epred.mcpfit.md)).

- Now works for 200+ characters formulas too.

- Fixed `cores = "all"` failed. Thanks for reporting,
  [@m-r-munroe](https://github.com/m-r-munroe)!

### Behind the scenes

In general, every effort has been made to anticipate future developments
and build the structure needed for that.

- The `mcpfamily` contains everything which was previously scattered
  across multiple functions. This means it is easy to add families and
  link functions - also enabling custom families in a future release.

- Major changes in how the model is translated into JAGS code, which is
  now module-like. The JAGS code is quite different but functionally
  equivalent. This opens for adding a new sampler in a future release.

- More thorough defensive coding.

- Much expanded test suite (now 10.000+ tests when run in full). The
  test suite now includes external validation of inference and
  simulation: AR against
  [`arima()`](https://rdrr.io/r/stats/arima.html)/[`arima.sim()`](https://rdrr.io/r/stats/arima.sim.html),
  binomial against [`glm()`](https://rdrr.io/r/stats/glm.html) /
  [`rbinom()`](https://rdrr.io/r/stats/Binomial.html), and changepoints
  against `{segmented}`.

- Fewer imports to userspace. This minimizes the risk of name conflicts.

## mcp 0.3.4

CRAN release: 2024-03-17

This is a bug fix release.

### Bug fixes

- Now respects the `cores` argument to
  [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md).
- Document all function arguments and remove documentation for removed
  arguments.

## mcp 0.3.3

CRAN release: 2023-03-22

This is a bug fix release.

### Bug fixes

- Support `ggplot >= 3.4.0`, `tidyselect >= 1.2.0`, and newer `future`
  by replacing deprecated functions.
- Accept `mcp(..., cores = "all")`.
- Fix documentation of `iter` argument to
  [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md).
- Other small fixes to deployment and documentation.

## mcp 0.3.2

CRAN release: 2022-02-18

This release contains no user-facing changes. The test suite suite is
now compatible with dplyr 1.0.8, which caused the test suite to fail.
This, in turn, would trigger the removal of mcp from CRAN.

## mcp 0.3.1

CRAN release: 2021-11-17

This is mostly a bug fix release.

### New features:

- `ex = mcp_example("demo", sample = TRUE)` is the new interface that
  replaces the `ex_*` datasets in prior versions. This reduces clutter
  of the namespace/documentation and the size of the package. It also
  gives the user richer details on the simulation and analyses. For
  “demo”, the `ex_demo` dataset is now `ex$data` and the `ex_fit` is
  `ex$fit`.

- Nicer printing of lists and texts all over. E.g., try
  `print(demo_fit$jags_code)` and `print(demo_fit$pars)`.

### Bug fixes

- Support breaking changes in `tidybayes >= 3.0.0` and `dplyr >= 1.0.6`

## mcp 0.3.0

CRAN release: 2020-08-03

### New features:

- Get fits and predictions for in-sample and out-of-sample data. [Read
  more in the article on these
  functions](https://lindeloev.github.io/mcp/articles/predict.html).

  - Use `predict(fit)` to get predicted values and quantiles.
  - Use `fitted(fit)` to get estimated values and quantiles.
  - Use `residuals(fit)` to get residuals and quantiles.

  All of the above functions include many arguments that align with (and
  extends) the options already in
  [`plot.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md),
  including getting fits/predictions for sigma (`which_y = "sigma"`),
  for the prior (`prior = TRUE`), and arbitrary quantiles
  (`probs = c(0.1, 0.5, 0.999)`). Use the `newdata` argument to get
  out-of-sample fitted/predicted values. Set `summary = FALSE` to get
  per-draw values.

- Added support for weighted regression for gaussian families:
  `model = list(y | weights(weight_column) ~ 1 + x)`. Weights are
  visualized as dot sizes in `plot(fit)`.

- Support for more link functions across families (e.g.,
  `family = gaussian(link = "log")`):

  - `gaussian`: “identity”, “log”
  - `binomial`: “logit”, “probit”, “identity”
  - `bernoulli`: “logit”, “probit”, “identity”
  - `poisson`: “log”, “identity”

- New argument `scale` in
  [`fitted()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md),
  and `fit$simulate()`. When `scale = "response"` (default), they return
  fits on the response scale. When `scale = "linear"`, they return fits
  on the linear-predictor (link) scale. Useful for model understanding
  and debugging.

- Use `pp_check(fit)` to do prior/posterior predictive checking. See
  `pp_check(fit, type = "x")` for a list of plot types.
  `pp_check(fit, facet_by = "group_column")` facets by a grouping
  column.

- Improvements to
  [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md):

  - Change-point densities are now computed on a per-panel basis in
    `plot(fit, facet_by = "group_column")`. Previous releases only
    displayed population-level change points.
  - You can now plot group-level effects with `rate = FALSE` for
    binomial models.
  - Change point densities in `plot(fit)` are not located directly on
    the x-axis. They were “floating” 5% above the x-axis in the previous
    releases.

- New argument `nsamples` reduces the number of samples used in most
  functions to speed up processing. `nsamples = NULL` uses all samples
  for maximum accuracy.

- New argument `arma` in many functions toggles whether autoregressive
  effects should be modelled.

- Although the API is still in alpha, feel free to try extracting
  samples using `mcp:::tidy_samples(fit)`. This is useful for further
  processing using `tidybayes`, `bayesplot`, etc. and is used
  extensively internally in `mcp`. One useful feature is computing
  absolute values for group-specific change points:
  `mcp:::tidy_samples(fit, population = FALSE, absolute = TRUE)`.
  Feedback is appreciated before `tidy_samples` becomes part of the
  `mcp` API in a future release.

### Other changes

- Change point densities in `plot(fit)` are now scaled to 20% of the
  plot for each chain X changepoint combo. This addresses a common
  problem where a wide posterior was almost invisibly low when a narrow
  posterior was present. This means that heights should only be compared
  *within* each chain x changepoint combo - not across.
- Removed the implicit ceiling of 1000 lines and samples in
  [`plot.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md).
- Rownames are removed from
  [`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
  and
  [`fixef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
  returns.
- A major effort has been put into making `mcp` robust and agile to
  develop. `mcp` now use defensive programming with helpful error
  messages. The Test suite includes 3600+ tests.
- [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md),
  [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
  etc. are now considerably faster for AR(N) due to vectorization of the
  underlying code.

### Bug fixes

- Sigma is now forced to stay positive via a floor at 0.
- Fixed: support and require dplyr 1.0.0. Now also requires tidybayes
  2.0.3.
- Fixed: Parallel sampling sometimes produced identical chains.
- Fixed several small bugs

## mcp 0.2.0

CRAN release: 2020-01-09

The API and internal structure should be stable now. v0.2.0 will be
released on CRAN.

### New features:

- Model quadratic and other terms using `I(x^2)`, `I(x^3.24)`, `sin(x)`,
  `sqrt(x)`, etc.

- Model residual standard deviation for `family = gaussian()` using
  `~ sigma([formula here])`.

- Model Nth order autoregressive models using `~ ar(order, formula)`,
  typically like `y ~ 1 + x + ar(2)` for AR(2). Simulate AR(N) models
  from scratch or given known data with `fit$simulate()`. The [article
  on AR(N)](https://lindeloev.github.io/mcp/articles/arma.html) has more
  details and examples. AR(N) models are popular to detect changes in
  time-series.

- Many updates to
  [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md).

  - Includes the posterior densities of the change point(s). Disable
    using `plot(fit, cp_dens = FALSE)`.
  - Supports AR(N) models (see above).
  - Plot posterior parameter intervals using `plot(fit, q_fit = TRUE)`.
    `plot(fit, q_fit = c(0.025, 0.5, 0.975))` plots 95% HDI and the
    median.
  - Plot posterior predictive intervals using
    `plot(fit, q_predict = TRUE)`.
  - Choose data geom. Currently takes “point” (default) and “line”
    (`plot(fit, geom_data = "line")`). The latter is useful for time
    series. Disable using `geom_data = FALSE`.

- Use `options(mc.cores = 3)` for considerable speed gains for the rest
  of the session. All vignettes/articles have been updated to recommend
  this as a default, though serial sampling is still the technical
  default. `mcp(..., cores = 3)` does the same thing on a call-by-ball
  basis.

- `fit$simulate()` adds the simulation parameters as an attribute
  (`attr(y, "simulate")`) to the predicted variable.
  [`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
  recognizes this and adds the simulated values to the results table
  (columns `sim` and `match`) so that one can inspect whether the values
  were recovered.

- Use `plot(fit, which_y = "sigma")` to plot the residual standard
  deviation on the y-axis. It works for AR(N) as well, e.g.,
  `which_y = "ar1"`, `which_y = "ar2"`, etc. This is useful to visualize
  changes in residual standard deviation and autocorrelation. The
  relevant vignettes have been updated with worked examples.

- Much love for the priors:

  - Set a Dirichlet prior on the change points using
    `prior = list(cp_1 = "dirichlet(1)", cp_2 = ...)`. [Read pros and
    cons here](https://lindeloev.github.io/mcp/articles/priors.html).
  - The default prior has been changed from truncated uniforms to a
    regularizing t-tail prior that samples effectively. [Read more
    here](https://lindeloev.github.io/mcp/articles/priors.html).
  - You can now sample the prior using `mcp(..., sample = "prior")` or
    `mcp(..., sample = "both")` and most methods can now take the prior:
    `plot(fit, prior = TRUE)`, `plot_pars(fit, prior = TRUE)`,
    `summary(fit, prior = TRUE)`, `ranef(fit, prior = TRUE)`.

- `mcp` can now be cited! Call `citation("mcp")` or see the pre-print
  here: <https://osf.io/fzqxv>.

### Other changes:

- Some renaming: “segments” –\> “model”. `fit$func_y()` –\>
  `fit$simulate()`.
- [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
  only visualize the total fit while
  [`plot_pars()`](https://lindeloev.github.io/mcp/dev/reference/plot_pars.md)
  only visualize individual parameters. These functions were mixed in
  [`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
  previously.
- The argument `update` has been discarded from
  [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) (it’s
  all on `adapt` now) and `inits` has been added.
- Many internal changes to make `mcp` more future proof. The biggest
  internal change is that `rjags` and `future` replace the `dclone`
  package. Among other things, this gives faster and cleaner
  installations.
- Many more informative error messages to help you quickly understand
  and solve errors.
- Updated documentation and website.

## mcp 0.1.0

First public release.

- Varying change points
- Basic GLM: Gaussian, binomial, Bernoulli, and Poisson, and associated
  vignettes.
- summary(fit), fixef(fit), and ranef(fit)
- plot(fit, “segments”) and plot(fit, “bayesplot-name-here”) with some
  options
- 1000+ basic unit tests to ensure non-breaking code for a wide variety
  of models.
- Testing and model comparison using `loo` and `hypothesis`
