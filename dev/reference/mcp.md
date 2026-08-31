# Fit Multiple Linear Segments And Their Change Points

Given a model (a list of segment formulas), `mcp` infers the posterior
distributions of the parameters of each segment as well as the change
points between segments. See details or [the mcp
website](https://lindeloev.github.io/mcp/).

## Usage

``` r
mcp(
  model,
  data,
  prior = list(),
  family = gaussian(),
  par_x = NULL,
  sample = "post",
  cores = NULL,
  chains = 3,
  iter = 3000,
  warmup = 1500,
  adapt = lifecycle::deprecated(),
  inits = NULL,
  jags_code = NULL,
  seed = NULL,
  diagnostics = list(),
  quiet = FALSE
)
```

## Arguments

- model:

  A list of formulas - one for each segment. The many examples on the
  [mcp website](https://lindeloev.github.io/mcp/). But briefly:

  The first formula has the format `response ~ predictors` while the
  following formulas have the format `response ~ cp ~ predictors`. Here,
  `cp` names the change-point part of the formula rather than a literal
  variable. The response and change-point parts can be omitted
  (`cp ~ predictor` assumes the same response; `~ predictor` assumes an
  intercept-only change point). Terms normally carry into later segments
  until redefined (see details).

  The following terms can be modeled:

  - *Regular formulas:* e.g., `~ 1 + x`. [Read
    more](https://lindeloev.github.io/mcp/articles/formulas.html).

  - *Extended formulas:* e.g., `~ x:group + I(x^2) + exp(z)`. [Read
    more](https://lindeloev.github.io/mcp/articles/formulas.html).
    R-side bases such as [`scale()`](https://rdrr.io/r/base/scale.html),
    [`poly()`](https://rdrr.io/r/stats/poly.html), and
    [`splines::ns()`](https://rdrr.io/r/splines/ns.html) are evaluated
    before sampling, and their fitted scaling or basis is reused for
    `newdata`.

  - *Group-level effects (random effects):* e.g., `~ 1 + (1 | id)` for a
    group-level intercept, or `~ 1 + (factor || id)` for independent
    intercept and factor-contrast deviations. [Read
    more](https://lindeloev.github.io/mcp/articles/group_effects.html).

  - *Gaussian residual standard deviation:* e.g., `~sigma(1)` for a
    simple standard-deviation change or `~sigma(1 + x + group)` for more
    advanced structures. Explicit
    [`sigma()`](https://rdrr.io/r/stats/sigma.html) formulas model
    log-SD, while the implicit constant `sigma_1` in a model without
    [`sigma()`](https://rdrr.io/r/stats/sigma.html) remains on the
    response scale. [Read
    more](https://lindeloev.github.io/mcp/articles/dpar.html)

  - *Time-series residuals:* use `ar(p)` and `ma(q)` separately or
    together, e.g., `~ 1 + ar(1) + ma(1)`. Both accept an optional
    regression formula and observation `boundary`. GARMA terms support
    Gaussian, binomial, Poisson, and negative-binomial families with
    their default links. They define a finite conditional recurrence on
    the link scale. For `~ ar(p) + ma(q)`, the model is:

    \$\$\eta_t = b_t + \sum\_{j=1}^{p} \phi\_{j,t}
    \left\[g(y^\*\_{t-j}) - b\_{t-j}\right\] + \sum\_{k=1}^{q}
    \theta\_{k,t} \left\[g(y^\*\_{t-k}) - \eta\_{t-k}\right\]\$\$

    where \\b_t\\ is the linear predictor from the segment formulas,
    \\\phi\_{j,t}\\ is the lag-\\j\\ autoregressive (AR) coefficient at
    time \\t\\, \\\theta\_{k,t}\\ is the lag-\\k\\ moving-average (MA)
    coefficient at time \\t\\, \\g(\cdot)\\ is the link function,
    \\y^\*\_t\\ is the observation, and \\\eta_t\\ is the resulting full
    linear predictor including serial dependence. Note some
    implications:

    - For an N-order component, the last N values *before* the segment
      onset are input to the first \\\eta_t\\ in the segment.

    - AR coefficients are not jointly constrained to stationarity; nor
      MA coefficients to invertibility. [Read
      more](https://lindeloev.github.io/mcp/articles/arma.html)

  - *Likelihood weights:* `y | weights(w) ~ ...` multiplies each
    observation's log-likelihood contribution by `w`, as in `brms`.
    Weights affect posterior inference and
    [`log_lik()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md),
    but not the response distribution used by
    [`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
    or prior/posterior predictive checks. Combine with other auxiliaries
    using `+`, e.g., `y | trials(total) + weights(w) ~ ...`.

  - *Binomial:* use `successes | trials(total) ~ ...` with
    `family = binomial()`.

- data:

  Table-like data in long format (data.frame, tibble, data.table, etc.)
  with syntactic column names. Missing values in the response variable
  are imputed using the posterior predictive.
  [`fitted.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  or
  [`predict.mcpfit`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  details how to see the imputed values.

- prior:

  Named list. Names are parameter names (`cp_i`, `Intercept_i`,
  `xvar_i`, `sigma_1`, etc.) and the values are either

  - A distribution in mcp's JAGS-string syntax (e.g.,
    `Intercept_1 = "dnorm(0, 1) T(0,)"`) indicating a conventional prior
    distribution. Data-calibrated, regularizing defaults are used where
    priors are not specified. These are designed for stable estimation
    and prediction, but should be justified before hypothesis testing.
    `mcp` uses conventional distribution scales rather than JAGS
    precision: SD for [`dnorm()`](https://rdrr.io/r/stats/Normal.html),
    scale for [`dt()`](https://rdrr.io/r/stats/TDist.html), `ddexp()`,
    and [`dlogis()`](https://rdrr.io/r/stats/Logistic.html), and log-SD
    for [`dlnorm()`](https://rdrr.io/r/stats/Lognormal.html). See
    details. With multiple change points, the default is a regularizing
    Student-t prior centered at `min(x)` and sequentially truncated
    between the preceding change point and `max(x)`. User-specified
    [`dunif()`](https://rdrr.io/r/stats/Uniform.html) priors and
    explicit truncations are used as written, so users must supply the
    ordering bounds, e.g., `dunif(cp_1, max(time))` for `cp_2`.

  - A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed
    value.

  - A model parameter name (e.g., `Intercept_2 = "Intercept_1"`),
    indicating that this parameter is shared - typically between
    segments. If two group-level deviations are shared this way, they
    will need to have the same grouping variable.

  - A scaled Dirichlet prior is supported for change points if they are
    all set to `cp_i = "dirichlet(alpha)"` with the same positive
    `alpha` for every change point. `alpha = 1` is most often used. This
    prior is symmetric over segment spacings, unlike the regularizing
    t-tail default for multiple change points, but it mixes less
    efficiently, so you will often need to set `iter` higher. It is
    recommended for hypothesis testing and for the estimation of more
    than 5 change points. [Read
    more](https://lindeloev.github.io/mcp/articles/priors.html).

- family:

  A supported family:
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`bernoulli()`](https://lindeloev.github.io/mcp/dev/reference/bernoulli.md),
  [`poisson()`](https://rdrr.io/r/stats/family.html), or
  [`negbinomial()`](https://lindeloev.github.io/mcp/dev/reference/negbinomial.md),
  with a supported link function; e.g., `gaussian(link = "log")`.

- par_x:

  String (default: `NULL` which is auto-detect).

- sample:

  One of

  - `"post"`: Sample the posterior.

  - `"prior"`: Sample only the prior. Plots, summaries, etc. will use
    the prior. This is useful for prior predictive checks.

  - `"both"`: Sample both prior and posterior. Plots, summaries, etc.
    will default to using the posterior. The prior only has effect when
    doing Savage-Dickey density ratios in
    [`hypothesis`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md).

  - `"none"` or `FALSE`: Do not sample. Returns an mcpfit object without
    sample. This is useful if you only want to check prior strings
    (`fit$prior`), the JAGS model (`fit$jags_code`), etc.

- cores:

  Deprecated and ignored. Configure parallel processing with a
  [future](https://future.futureverse.org/reference/plan.html) plan
  instead, for example
  `future::plan(future::multisession, workers = 3)`. With the default
  future plan, chains are sampled sequentially. The argument remains
  available for backwards compatibility.

- chains:

  Positive integer. Number of chains to run.

- iter:

  Positive integer. Number of post-warmup draws from each chain. The
  total number of draws is `iter * chains`.

- warmup:

  Positive integer. Number of initial iterations per chain which are
  discarded before sampling. Set higher if needed for sampler
  adaptation; use diagnostics to assess convergence.

- adapt:

  Deprecated; use `warmup` instead.

- inits:

  A list if initial values for the parameters. This can be useful if a
  model fails to converge. Read more in
  [`jags.model`](https://rdrr.io/pkg/rjags/man/jags.model.html).
  Defaults to `NULL`, i.e., no inits.

- jags_code:

  String. Pass JAGS code to `mcp` to use directly. This is useful if you
  want to tweak the code in `fit$jags_code` and run it within the `mcp`
  framework. R-side simulation and prediction methods continue to use
  the mcp-default formulas (with warning), so they may no longer match
  custom JAGS code.

- seed:

  `NULL` or a positive integer. Seed for the JAGS random-number
  generators. When supplied, `inits` must be a single named list shared
  by all chains.

- diagnostics:

  Named list of diagnostic warning thresholds. Available elements are
  `rhat = 1.01`, `ess_bulk = 400`, `ess_tail = 400`, `ar = 0.10`, and
  `ma = 0.10`. An empty list uses these defaults; a partial list
  overrides only the supplied values. Set an element to `NULL` to
  disable that diagnostic, or use `FALSE` to disable all configurable
  diagnostic warnings. In
  [`summary.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md),
  `NULL` inherits the settings used to fit the model, while a list or
  `FALSE` overrides the diagnostic footer.

- quiet:

  Logical. Suppress routine JAGS output and mcp sampling-status
  messages? Defaults to `FALSE`.

## Value

An
[`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
object.

## Details

**The mcp model**

For a continuous predictor \\x\\, segment 1 (\\x \le \Delta_1\\) is a
standard regression model with the intercept at \\x = 0\\:

\$\$\eta\_{1,i} = f_1(x_i, \mathbf{\beta}\_1) = \beta\_{1,0} +
\beta\_{1,1} x_i + \dots\$\$

Subsequent segments \\k \in \\2, \dots, K\\\\ are similar, but replace
\\x\\ with a segment-local predictor \\X\_{k,i}\\ measured from the
segment onset \\\Delta\_{k-1}\\ and capped at the segment end
\\\Delta_k\\:

\$\$X\_{k,i} = \min(x_i, \Delta_k) - \Delta\_{k-1}\$\$

with \\\Delta_K = \max(\mathbf{x})\\.

Joined segments (`~ 0 + x`) continue from the plateaued level of earlier
segments with a new segment slope \\\beta\_{k,1} X\_{k,i}\\, enforcing
continuity without parameter constraints. Disjoined segments (`~ 1 + x`)
introduce a new intercept \\\beta\_{k,0}\\ at \\\Delta\_{k-1}\\ and
truncate preceding segments. In both cases, slope and intercept
parameters are absolute values (not differences relative to the previous
segment).

Here, the model was presented for the mean (on the link scale). The
exact same model apply to distributional parameters
([`sigma()`](https://rdrr.io/r/stats/sigma.html), `shape()`, etc.) and
[`ar()`](https://rdrr.io/r/stats/ar.html)/`ma()` too. See more details
on the `mcp` model in
[mcp-package](https://lindeloev.github.io/mcp/dev/reference/mcp-package.md)
and on the [mcp
website](https://lindeloev.github.io/mcp/articles/formulas.html).

**Notes on priors**

- *Ordered change point priors:* Default population-level `cp_i` priors
  are ordered and the ordering is imposed through the priors. For
  user-defined priors, `mcp` adds truncation (e.g., `T(cp_1, )`) only
  when the prior has neither explicit truncation nor an inherently
  bounded form such as [`dunif()`](https://rdrr.io/r/stats/Uniform.html)
  or `dirichlet()`.

- *Data-dependent terms:* If `mcp` encounters a data-dependent term like
  `min(time)`, `max(time)`, `median(response)`, or `mad(response)` in
  the prior string, they are resolved from the model data so a numerical
  value is passed to JAGS. The following terms are also allowed:
  `n_segments()`, and `n_cp()`. The older constants `MINX`, `MAXX`,
  `MEANX`, `SDX`, `MINY`, `MAXY`, `MEANY`, `SDY`, and `N_CP` remain
  accepted with a deprecation warning.

- *Parameterization:* Prior strings use conventional scale
  parameterizations. `mcp` converts these to the parameterization
  required by JAGS when generating code: inverse variance for
  [`dnorm()`](https://rdrr.io/r/stats/Normal.html),
  [`dt()`](https://rdrr.io/r/stats/TDist.html), and
  [`dlnorm()`](https://rdrr.io/r/stats/Lognormal.html), and inverse
  scale for `ddexp()` and
  [`dlogis()`](https://rdrr.io/r/stats/Logistic.html). Use
  `prior_summary(fit)` for resolved priors and
  `prior_summary(fit, verbose = TRUE)` for their rules and descriptions.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# \donttest{
# Define the segments using formulas. A change point is estimated between each formula.
model = list(
  response ~ 1,  # Plateau in the first segment (Intercept_1)
  ~ 0 + time,    # Joined slope (time_2) at cp_1
  ~ 1 + time     # Disjoined slope (Intercept_3, time_3) at cp_2
)

# Fit it and sample the prior too.
# future::plan(future::multisession, workers = 3)  # Uncomment for parallel sampling
data = mcp_example_data("demo")  # Simulated data example
demo_fit = mcp(model, data = data, sample = "both", seed = 42)

# See parameter estimates
summary(demo_fit)
#> Family: gaussian
#> Links: mu = identity; sigma = identity
#> Iterations: 3000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  cp_1        30.85 3.586 24.10 38.315 1.01      407      541 30.0    OK
#>  cp_2        69.78 0.291 69.30 70.261 1.00     5745     7530 70.0    OK
#> 
#> Population-level parameters:
#>         name  mean    sd lower  upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 10.32 0.710  8.91 11.711 1.00     1406     2753 10.0    OK
#>  time_2       0.54 0.065  0.43  0.682 1.01      438      604  0.5    OK
#>  Intercept_3  0.79 1.531 -2.10  3.856 1.00      749     1389  0.0    OK
#>  time_3      -0.24 0.091 -0.42 -0.067 1.00      736     1300 -0.2    OK
#>  sigma_1      4.00 0.299  3.49  4.677 1.00     3844     3900  4.0    OK

# Visual inspection of the results
plot(demo_fit)  # Visualization of model fit/predictions

plot_pars(demo_fit)  # Parameter distributions


pp_check(demo_fit)  # Prior/Posterior predictive checks


# Test a hypothesis
hypothesis(demo_fit, "cp_1 > 10")
#>      hypothesis     mean    lower    upper p  BF
#> 1 cp_1 - 10 > 0 20.84596 14.09683 28.31539 1 Inf

# Make predictions
head(fitted(demo_fit))
#>    response     time    fitted     error      Q2.5     Q97.5
#> 1 32.842651 68.35820 30.374492 1.0363972 28.362893 32.453365
#> 2 -1.160003 87.29038 -3.333719 0.7695936 -4.836521 -1.822012
#> 3 27.564248 69.01173 30.727352 1.0675322 28.660966 32.861266
#> 4 10.062971 11.59361 10.316267 0.7101668  8.907415 11.711335
#> 5 14.056859 19.50091 10.316619 0.7096133  8.907914 11.711335
#> 6 18.292640 46.12009 18.367438 0.9842136 16.134892 20.014110
head(predict(demo_fit))
#>    response     time   predict    error       Q2.5     Q97.5
#> 1 32.842651 68.35820 30.353477 4.210965  22.226246 38.515625
#> 2 -1.160003 87.29038 -3.299188 4.064416 -11.355322  4.704244
#> 3 27.564248 69.01173 30.669097 4.170328  22.564006 38.884471
#> 4 10.062971 11.59361 10.352696 4.047950   2.308171 18.325861
#> 5 14.056859 19.50091 10.346765 4.086618   2.308800 18.326099
#> 6 18.292640 46.12009 18.332877 4.100923  10.222899 26.465992
head(predict(demo_fit, newdata = data.frame(time = c(55.545, 80, 132))))
#>      time    predict    error       Q2.5     Q97.5
#> 1  55.545  23.465561 4.119947  15.440778 31.460977
#> 2  80.000  -1.670004 4.092160  -9.679104  6.435948
#> 3 132.000 -13.894364 5.927965 -25.551752 -2.337486

# Compare to a one-intercept-only model (no change points) with default prior
model_null = list(response ~ 1)
fit_null = mcp(model_null, data = data, par_x = "time")  # fit another model here
demo_loo = loo(demo_fit)
null_loo = loo(fit_null)
loo::loo_compare(demo_loo, null_loo)
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model1       0.0     0.0      NA                    
#>  model2    -108.1     8.8    1.00                    

# Inspect the prior. Useful for prior predictive checks.
summary(demo_fit, prior = TRUE)
#> Family: gaussian
#> Links: mu = identity; sigma = identity
#> Iterations: 3000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name     mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        37.45617 24.88   4.66 91.79 1.00     8902     9027 30.0    OK
#>  cp_2        60.99546 24.31  14.53 96.92 1.00     9183     8085 70.0    OK
#> 
#> Population-level parameters:
#>         name     mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 10.53694 26.30 -29.43 52.10 1.00     9086     8947 10.0    OK
#>  time_2       0.00334  0.24  -0.41  0.43 1.00     9363     8823  0.5      
#>  Intercept_3 10.54960 22.17 -30.77 51.53 1.00     8828     8629  0.0    OK
#>  time_3      -0.00047  0.22  -0.41  0.43 1.00     9183     8660 -0.2    OK
#>  sigma_1     14.04936 16.94   0.45 53.68 1.00     8637     8735  4.0    OK
plot(demo_fit, prior = TRUE)


# Show all priors. Default priors are added where you don't provide any
prior_summary(demo_fit)
#> # A tibble: 7 × 5
#>   parameter   segment dpar  prior                                         bounds
#>   <chr>         <int> <chr> <chr>                                         <chr> 
#> 1 cp_1              2 cp    student_t(df = 1, location = 3.189323, scale… [min(…
#> 2 cp_2              3 cp    student_t(df = 1, location = 3.189323, scale… [cp_1…
#> 3 Intercept_1       1 mu    student_t(df = 3, location = 10.4, scale = 1… none  
#> 4 time_2            2 mu    student_t(df = 3, location = 0, scale = 0.13… none  
#> 5 Intercept_3       3 mu    student_t(df = 3, location = 10.4, scale = 1… none  
#> 6 time_3            3 mu    student_t(df = 3, location = 0, scale = 0.13… none  
#> 7 sigma_1           1 sigma student_t(df = 3, location = 0, scale = 12.8) [0, I…

# Set priors and re-run
prior = list(
  Intercept_1 = 15,
  time_2 = "dt(0, 2, 1) T(0, )",  # t-dist slope. Truncated to positive.
  cp_2 = "dunif(cp_1, 80)",    # change point to segment 2 > cp_1 and < 80.
  Intercept_3 = "Intercept_1"           # Shared intercept between segment 1 and 3
)

fit3 = mcp(model, data = data, prior = prior)

# Show the JAGS model
demo_fit$jags_code
#> model {
#>   # mcp helper values
#>   cp_0 = CONST1_
#>   cp_3 = CONST2_
#> 
#>   # Priors for population-level effects
#>   cp_1 ~ dt(CONST1_, 1/(CONST3_)^2, CONST4_) T(CONST1_,CONST2_)  # Regularizing t-tail within the observed change-point span
#>   cp_2 ~ dt(CONST1_, 1/(CONST3_)^2, CONST4_) T(cp_1,CONST2_)  # Regularizing t-tail ordered after cp_1 within the observed change-point span
#>   Intercept_1 ~ dt(10.4, 1/(12.8)^2, 3)   # Robustly centered mean intercept with a minimum scale of 2.5
#>   time_2 ~ dt(0, 1/(0.1350637)^2, 3)   # Regularizing mean coefficient scaled to a reference predictor change
#>   Intercept_3 ~ dt(10.4, 1/(12.8)^2, 3)   # Robustly centered mean intercept with a minimum scale of 2.5
#>   time_3 ~ dt(0, 1/(0.1350637)^2, 3)   # Regularizing mean coefficient scaled to a reference predictor change
#>   sigma_1 ~ dt(0, 1/(12.8)^2, 3) T(0,)  # Positive residual SD calibrated on the response scale
#> 
#>   # Model and likelihood
#>   for (i_ in 1:length(time)) {
#>     # par_x local to each segment
#>     x_local_1_[i_] = min(time[i_], cp_1)
#>     x_local_2_[i_] = min(time[i_], cp_2) - cp_1
#>     x_local_3_[i_] = min(time[i_], cp_3) - cp_2
#>     
#>     # Formula for mu
#>     link_mu_[i_] =
#>       (time[i_] >= cp_0) * (time[i_] < cp_2) * inprod(rhs_matrix_[i_, c(1)], c(Intercept_1)) * 1 + 
#>       (time[i_] >= cp_1) * (time[i_] < cp_2) * inprod(rhs_matrix_[i_, c(2)], c(time_2)) * x_local_2_[i_] + 
#>       (time[i_] >= cp_2) * inprod(rhs_matrix_[i_, c(3)], c(Intercept_3)) * 1 + 
#>       (time[i_] >= cp_2) * inprod(rhs_matrix_[i_, c(4)], c(time_3)) * x_local_3_[i_]
#>     
#>     # Formula for sigma
#>     link_sigma_[i_] =
#>       (time[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(5)], c(sigma_1)) * 1
#> 
#>     # Likelihood and log-density for family = gaussian()
#>     mu_[i_] = link_mu_[i_]
#>     sigma_[i_] = max(1e-03, link_sigma_[i_])
#>     response[i_] ~ dnorm(mu_[i_], 1 / sigma_[i_]^2)
#>   }
#> }
# }
```
