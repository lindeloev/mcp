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

  - *Weights:* `y | weights(w) ~ ...` specifies observation
    log-likelihood weights.

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
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 100
#>    Unobserved stochastic nodes: 7
#>    Total graph size: 2139
#> 
#> Initializing model
#> 
#> Finished sampling in 2.4 seconds
#> Warning: Some parameters may not have converged well:
#>   * rhat > 1.01 or ess_bulk < 400 or ess_tail < 400: cp_1 and time_2
#> Inspect `summary(fit)` and `plot_pars(fit)`, and consider increasing `iter`/`warmup` or simplifying the model before trusting these results.
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 0
#>    Unobserved stochastic nodes: 107
#>    Total graph size: 2139
#> 
#> Initializing model
#> 
#> Finished sampling in 0 seconds

# See parameter estimates
summary(demo_fit)
#> Family: gaussian(link = 'identity')
#> Iterations: 3000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name  mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        23.10 5.203 12.75 32.82    1      303      414 30.0    OK
#>  cp_2        69.91 0.342 69.35 70.48    1     5513     7615 70.0    OK
#> 
#> Population-level parameters:
#>         name  mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.00 0.928  7.06 10.69    1      454      471 10.0    OK
#>  time_2       0.39 0.051  0.30  0.51    1      387      593  0.5    OK
#>  Intercept_3  1.71 1.212 -0.68  4.09    1      824     1811  0.0    OK
#>  time_3      -0.24 0.067 -0.36 -0.11    1      832     1754 -0.2    OK
#>  sigma_1      3.68 0.273  3.20  4.25    1     4224     4344  4.0    OK
#> 
#> Warning: 2 parameters show poor convergence (rhat > 1.01 or ess_bulk < 400 or ess_tail < 400).

# Visual inspection of the results
plot(demo_fit)  # Visualization of model fit/predictions

plot_pars(demo_fit)  # Parameter distributions


pp_check(demo_fit)  # Prior/Posterior predictive checks


# Test a hypothesis
hypothesis(demo_fit, "cp_1 > 10")
#>      hypothesis     mean    lower   upper         p       BF
#> 1 cp_1 - 10 > 0 13.10434 2.750022 22.8215 0.9877778 17.68719

# Make predictions
head(fitted(demo_fit))
#>     response     time    fitted     error      Q2.5      Q97.5
#> 1 -3.0084198 91.48060 -3.386119 0.7361934 -4.828339 -1.9431978
#> 2 -7.8768640 93.70754 -3.911807 0.8265426 -5.540613 -2.2927652
#> 3 16.3029101 28.61395 11.123294 1.0336692  9.164108 13.0789150
#> 4 -0.0373553 83.04476 -1.394763 0.6392671 -2.640529 -0.1318978
#> 5 27.4463185 64.17455 24.940974 0.8865290 23.196406 26.6955478
#> 6 22.0610004 51.90959 20.115305 0.6000163 18.922358 21.2736971
head(predict(demo_fit))
#>     response     time   predict    error       Q2.5     Q97.5
#> 1 -3.0084198 91.48060 -3.356702 3.719059 -10.713765  3.903302
#> 2 -7.8768640 93.70754 -3.915844 3.781332 -11.078725  3.573730
#> 3 16.3029101 28.61395 11.120751 3.848448   3.588405 18.597786
#> 4 -0.0373553 83.04476 -1.408135 3.708950  -8.747380  5.835234
#> 5 27.4463185 64.17455 24.925625 3.767051  17.515536 32.226158
#> 6 22.0610004 51.90959 20.080061 3.770960  12.562271 27.548763
head(predict(demo_fit, newdata = data.frame(time = c(55.545, 80, 132))))
#>      time    predict    error       Q2.5     Q97.5
#> 1  55.545  21.553020 3.763265  14.115814 28.925727
#> 2  80.000  -0.643233 3.790046  -8.129513  6.735050
#> 3 132.000 -12.898969 4.891453 -22.486915 -3.313545

# Compare to a one-intercept-only model (no change points) with default prior
model_null = list(response ~ 1)
fit_null = mcp(model_null, data = data, par_x = "time")  # fit another model here
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 100
#>    Unobserved stochastic nodes: 2
#>    Total graph size: 622
#> 
#> Initializing model
#> 
#> Finished sampling in 0.4 seconds
demo_loo = loo(demo_fit)
null_loo = loo(fit_null)
loo::loo_compare(demo_loo, null_loo)
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model1       0.0     0.0      NA                    
#>  model2    -104.3     8.5    1.00                    

# Inspect the prior. Useful for prior predictive checks.
summary(demo_fit, prior = TRUE)
#> Family: gaussian(link = 'identity')
#> Iterations: 3000 from 3 chains.
#> Segments:
#>   1: response ~ 1
#>   2: response ~ 1 ~ 0 + time
#>   3: response ~ 1 ~ 1 + time
#> 
#> Change point parameters:
#>         name     mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        35.77146 25.96   1.56 92.46    1     8902     9027 30.0    OK
#>  cp_2        60.32791 25.36  11.85 97.80    1     9183     8085 70.0    OK
#> 
#> Population-level parameters:
#>         name     mean    sd  lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1  9.61982 23.02 -25.35 45.98    1     9086     8947 10.0    OK
#>  time_2       0.00280  0.20  -0.35  0.36    1     9363     8823  0.5      
#>  Intercept_3  9.63090 19.40 -26.52 45.49    1     8828     8629  0.0    OK
#>  time_3      -0.00039  0.19  -0.35  0.36    1     9183     8660 -0.2    OK
#>  sigma_1     12.29319 14.82   0.39 46.97    1     8637     8735  4.0    OK
plot(demo_fit, prior = TRUE)


# Show all priors. Default priors are added where you don't provide any
print(demo_fit$prior)
#> List of 7
#>  $ cp_1       :"dt(0.02388966, 49.43264, 1) T(0.02388966, 98.88917)"
#>  $ cp_2       :"dt(0.02388966, 49.43264, 1) T(cp_1, 98.88917)"
#>  $ Intercept_1:"dt(9.5, 11.2, 3)"
#>  $ time_2     :"dt(0, 0.1132855, 3)"
#>  $ Intercept_3:"dt(9.5, 11.2, 3)"
#>  $ time_3     :"dt(0, 0.1132855, 3)"
#>  $ sigma_1    :"dt(0, 11.2, 3) T(0, )"

# Set priors and re-run
prior = list(
  Intercept_1 = 15,
  time_2 = "dt(0, 2, 1) T(0, )",  # t-dist slope. Truncated to positive.
  cp_2 = "dunif(cp_1, 80)",    # change point to segment 2 > cp_1 and < 80.
  Intercept_3 = "Intercept_1"           # Shared intercept between segment 1 and 3
)

fit3 = mcp(model, data = data, prior = prior)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 100
#>    Unobserved stochastic nodes: 5
#>    Total graph size: 2138
#> 
#> Initializing model
#> 
#> Finished sampling in 1.8 seconds

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
#>   Intercept_1 ~ dt(9.5, 1/(11.2)^2, 3)   # Robustly centered mean intercept with a minimum scale of 2.5
#>   time_2 ~ dt(0, 1/(0.1132855)^2, 3)   # Regularizing mean coefficient scaled to a reference predictor change
#>   Intercept_3 ~ dt(9.5, 1/(11.2)^2, 3)   # Robustly centered mean intercept with a minimum scale of 2.5
#>   time_3 ~ dt(0, 1/(0.1132855)^2, 3)   # Regularizing mean coefficient scaled to a reference predictor change
#>   sigma_1 ~ dt(0, 1/(11.2)^2, 3) T(0,)  # Positive residual SD calibrated on the response scale
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
#>     
#>       # Segment 1: response ~ 1
#>       (time[i_] >= cp_0) * (time[i_] < cp_2) * inprod(rhs_matrix_[i_, c(1)], c(Intercept_1)) * 1 + 
#>     
#>       # Segment 2: response ~ 1 ~ 0 + time
#>       (time[i_] >= cp_1) * (time[i_] < cp_2) * inprod(rhs_matrix_[i_, c(2)], c(time_2)) * x_local_2_[i_] + 
#>     
#>       # Segment 3: response ~ 1 ~ 1 + time
#>       (time[i_] >= cp_2) * inprod(rhs_matrix_[i_, c(3)], c(Intercept_3)) * 1 + 
#>       (time[i_] >= cp_2) * inprod(rhs_matrix_[i_, c(4)], c(time_3)) * x_local_3_[i_]
#>     
#>     # Formula for sigma
#>     link_sigma_[i_] =
#>     
#>       # Segment 1: response ~ 1
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
