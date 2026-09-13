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

  A list of formulas - one for each segment. The general format is
  `response ~ cp ~ predictors` (e.g., `y ~ 1 ~ 1 + x`), except the first
  segment has no change point and uses `response ~ predictors`. The
  response and change-point parts can be omitted (`cp ~ predictor`
  assumes the same response; `~ predictor` assumes an intercept-only
  change point). Non-\$x\$ terms persist into later segments until
  replaced or removed by an intercept reset (see details). See examples
  on the [mcp website](https://lindeloev.github.io/mcp/).

  **1. Response (segment 1 only):**

  - `y ~ ...`: Standard continuous or count response (Gaussian, Poisson,
    Bernoulli).

  - `successes | trials(total) ~ ...`: Binomial response
    (`family = binomial()`).

  - `y | weights(w) ~ ...`: Observation log-likelihood weights
    (multiplies each observation's log-likelihood contribution by
    `w > 0`; affects posterior inference and
    [`log_lik()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
    but not predictions).

  - `y | trials(total) + weights(w) ~ ...`: Combine response auxiliaries
    using `+`.

  **2. Change-point modeling (`cp`, segments 2+):**

  - `~ 1 ~ ...` (or omitted, e.g., `~ x`): Population-level change point
    (default).

  - `1 + (1 | id) ~ ...`: Group-level change-point deviations around the
    population change point. [Read
    more](https://lindeloev.github.io/mcp/articles/group_effects.html).

  **3. Regression formula (all segments):** [Read
  more](https://lindeloev.github.io/mcp/articles/formulas.html)

  - `~ 1 + x`: Disjoined slope with a new segment intercept.

  - `~ 0 + x`: Joined slope (no intercept; continuous from previous
    segment).

  - `~ 1`: Plateau (intercept only, no slope).

  - `~ x:group + I(x^2) + exp(z)`: Extended terms, interactions, and
    R-side bases ([`scale()`](https://rdrr.io/r/base/scale.html),
    [`poly()`](https://rdrr.io/r/stats/poly.html),
    [`splines::ns()`](https://rdrr.io/r/splines/ns.html)). Bases are
    evaluated before sampling and reused for `newdata`.

  - `~ 1 + (1 | id)`: Group-level intercepts (or `(1 + x || id)` for
    independent slopes and intercepts). [Read
    more](https://lindeloev.github.io/mcp/articles/group_effects.html).

  - `~ sigma(1 + x)`: Distributional parameters on the link scale (e.g.,
    log residual SD). [Read
    more](https://lindeloev.github.io/mcp/articles/dpar.html).

  - `~ ar(1) + ma(1)`: Autoregressive and moving-average time-series
    residuals on the link scale (accepts regression formulas,
    `series = id`, and `boundary`; use `ar(0)` or `ma(0)` to turn off in
    later segments). [Read
    more](https://lindeloev.github.io/mcp/articles/arma.html).

- data:

  Table-like data in long format (data.frame, tibble, data.table, etc.)
  with syntactic column names. Missing values in the response variable
  are imputed using the posterior predictive.
  [`fitted.mcpfit`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  or
  [`predict.mcpfit`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  details how to see the imputed values.

- prior:

  Named list. Names are parameter names (`cp_i`, `Intercept_i`,
  `xvar_i`, `sigma_1`, etc.) and the values are either

  - A distribution in `mcp`'s JAGS-string syntax (e.g.,
    `Intercept_1 = "dnorm(0, 1) T(0,)"`) indicating a conventional prior
    distribution. Data-calibrated, regularizing defaults are used where
    priors are not specified. These are designed for stable estimation
    and prediction, but should be justified before hypothesis testing.
    `mcp` uses conventional distribution scales rather than JAGS
    precision: SD for [`dnorm()`](https://rdrr.io/r/stats/Normal.html),
    scale for [`dt()`](https://rdrr.io/r/stats/TDist.html), `ddexp()`,
    and [`dlogis()`](https://rdrr.io/r/stats/Logistic.html), and log-SD
    for [`dlnorm()`](https://rdrr.io/r/stats/Lognormal.html). See
    details.

  - A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed
    value.

  - A model parameter name (e.g., `Intercept_2 = "Intercept_1"`),
    indicating that this parameter is shared - typically between
    segments. If two group-level deviations are shared this way, they
    will need to have the same grouping variable.

  - The default prior on change points is `dirichlet(1)` (uniform order
    statistics). For a single change point, this is the Beta(1, 1) /
    Uniform distribution over `[min(x), max(x)]`. For multiple change
    points, it corresponds to a flat Dirichlet distribution over segment
    lengths. You can also explicitly set `cp_i = "dirichlet(alpha)"`
    with the same positive `alpha` for all change points to regularize
    spacing (`alpha > 1` penalizes change points from occurring close
    together, while `alpha < 1` favors clustering). Under the hood, this
    is parameterized as an exact sequential stick-breaking Beta chain
    for fast and robust sampling. [Read
    more](https://lindeloev.github.io/mcp/articles/priors.html).

- family:

  A supported family:
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`bernoulli()`](https://lindeloev.github.io/mcp/reference/bernoulli.md),
  [`poisson()`](https://rdrr.io/r/stats/family.html), or
  [`negbinomial()`](https://lindeloev.github.io/mcp/reference/negbinomial.md),
  with a supported link function; e.g., `gaussian(link = "log")`.

- par_x:

  String (default: `NULL` which is auto-detect).

- sample:

  One of

  - `"post"`: Sample the posterior.

  - `"prior"`: Sample only the prior. Use `prior = TRUE` explicitly in
    plots, summaries, and other methods to select prior draws.

  - `"both"`: Sample both prior and posterior. Plots, summaries, etc.
    will default to using the posterior. Use `prior = TRUE` to select
    prior draws.

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

  A list of initial values for the parameters. This can be useful if a
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
  generators.

- diagnostics:

  Named list of diagnostic warning thresholds. Available elements are
  `rhat = 1.01`, `ess_bulk = 400`, `ess_tail = 400`, `ar = 0.10`, and
  `ma = 0.10`. An empty list uses these defaults; a partial list
  overrides only the supplied values. Set an element to `NULL` to
  disable that diagnostic, or use `FALSE` to disable all configurable
  diagnostic warnings. In
  [`summary.mcpfit()`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md),
  `NULL` inherits the settings used to fit the model, while a list or
  `FALSE` overrides the diagnostic footer.

- quiet:

  Logical. Suppress routine JAGS output and mcp sampling-status
  messages? Defaults to `FALSE`.

## Value

An [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
object.

## The mcp model

Consider the following model which you can find in `demo_fit` and
`mcp_example("demo")`:

    model = list(
      response ~ 1,  # Plateau in the first segment (Intercept_1)
      ~ 0 + time,    # Joined slope (time_2) in segment 2 which starts at cp_1
      ~ 1 + time     # Disjoined slope (Intercept_3, time_3) at cp_2
    )

![Fitted 3-segment mcp model with a plateau, joined slope, and disjoined
slope](figures/mcp_demo.png)

This model has \\K=3\\ segments separated by \\K-1=2\\ change points:
\\\tau_1\\ and \\\tau_2\\.

More generally, an `mcp` model divides a continuous predictor \\x\\ into
\\K\\ segments separated by ordered change points \\\tau_1 \< \dots \<
\tau\_{K-1}\\. In each segment \\k \in \\1, \dots, K\\\\, the linear
predictor \\\eta_i\\ is evaluated directly from the segment-local
distance \\(x_i - \tau\_{k-1})\\:

\$\$\eta_i = \alpha_k + \beta\_{k,1} (x_i - \tau\_{k-1}) \quad
(\text{with } \tau_0 = 0)\$\$

where the segment-start level \\\alpha_k\\ is freely estimated for the
first and disjoined segments, as in non-segmented regression, and
determined by continuity for joined segments:

\$\$\alpha_k = \begin{cases} \beta\_{k,0}, & \text{Disjoined segments }
(\sim \texttt{1 + x}, \text{ including } k = 1) \\ \alpha\_{k-1} +
\beta\_{k-1,1} (\tau\_{k-1} - \tau\_{k-2}), & \text{Joined segments } (k
\ge 2, \sim \texttt{0 + x}) \end{cases}\$\$

Here, \\\beta\_{k,0}\\ is the segment-start intercept, and
\\\beta\_{k,1}\\ is the slope on \\x\\. In all segments, estimated slope
and intercept parameters are absolute values (not changes relative to
the preceding segment).

If additional continuous covariates or categorical factors are included
(e.g., `+ z + group`), they enter additively on their original scale
(\\\dots + \sum \gamma\_{k,j} z\_{j,i}\\ for covariate \\j\\); only bare
change-point predictor terms \\x\\ and `I(x^k)` are converted to
segment-local coordinates.

The inverse-linked parameter is \\\mu_i = g^{-1}(\eta_i)\\ via link
function \\g(\mu_i) = \eta_i\\, representing the expected response for
most families (or the success probability for binomial models, where the
expected count is \\n_i \mu_i\\). Distributional parameters
([`sigma()`](https://rdrr.io/r/stats/sigma.html), `shape()`, etc.) and
autoregressive terms ([`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`)
follow this exact same segmented structure on their respective link
scales. See more details on the `mcp` model in
[mcp-package](https://lindeloev.github.io/mcp/reference/mcp-package.md)
and on the [mcp
website](https://lindeloev.github.io/mcp/articles/formulas.html).

## Time-series residuals (link-scale observation-driven GARMA)

Autoregressive (`ar(p)`) and moving-average (`ma(q)`) terms define a
finite conditional recurrence on the link scale (generalized
autoregressive moving-average, GARMA). They support Gaussian
(`identity`), binomial (`logit`), Bernoulli (`logit`), Poisson (`log`),
and negative-binomial (`log`) families. If \\\eta^{\text{reg}}\_t\\ is
the ordinary regression predictor from the segment formulas and
\\\eta_t\\ is the predictor including serial dependence, the recurrence
decomposes into components:

\$\$\begin{aligned} \text{AR}\_t &= \sum\_{j=1}^{p} \phi\_{j,t}
\left\[g(y^\*\_{t-j}) - \eta^{\text{reg}}\_{t-j}\right\] \\ \text{MA}\_t
&= \sum\_{k=1}^{q} \theta\_{k,t} \left\[g(y^\*\_{t-k}) -
\eta\_{t-k}\right\] \\ \eta_t &= \eta^{\text{reg}}\_t + \text{AR}\_t +
\text{MA}\_t \end{aligned}\$\$

where \\\phi\_{j,t}\\ is the lag-\\j\\ autoregressive (AR) coefficient
at time \\t\\, \\\theta\_{k,t}\\ is the lag-\\k\\ moving-average (MA)
coefficient at time \\t\\, \\g(\cdot)\\ is the link function, and
\\y^\*\_t\\ is the boundary-constrained observation with pseudo-count
\\b\\ (set via argument `boundary = 0.1` in
[`ar()`](https://rdrr.io/r/stats/ar.html) / `ma()`) to keep residuals
finite on the link scale:

- **Gaussian:** \\y^\*\_t = y_t\\.

- **Poisson / Negative Binomial:** \\y^\*\_t = \max(y_t, b)\\ to prevent
  \\\log(0)\\. Here \\b\\ replaces zero counts with a small positive
  count.

- **Binomial / Bernoulli:** \\y^\*\_t = \min(\max(y_t, b), n_t - b) /
  n_t\\, where \\y_t\\ is observed successes, \\n_t\\ is the number of
  trials (\\n_t = 1\\ for Bernoulli), and \\b\\ clamps counts to the
  interval \\\[b, n_t - b\]\\ before converting to a rate, preventing
  \\\text{logit}(0)\\ and \\\text{logit}(1)\\.

Implications:

- For an \\N\\-order component, the last \\N\\ values *before* the
  segment onset are input to the first \\\eta_t\\ in the segment.

- AR and MA components persist into later segments until replaced or
  turned off via `ar(0)` or `ma(0)`.

- AR coefficients are not jointly constrained to stationarity; nor MA
  coefficients to invertibility.

- See [the arma
  article](https://lindeloev.github.io/mcp/articles/arma.html) for more
  details.

## Notes on priors

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
  `n_segments()` and `n_cp()`. The older constants `MINX`, `MAXX`,
  `MEANX`, `SDX`, `MINY`, `MAXY`, `MEANY`, `SDY`, and `N_CP` remain
  accepted with a deprecation warning.

- *Group-level change points:* Group-specific locations follow a
  hierarchical normal distribution around their population change point,
  truncated so that realized locations remain in the observed range and
  ordered.

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

## References

- Lindeløv, J. K. (2020). mcp: An R Package for Regression With Multiple
  Change Points. *OSF Preprints*.
  [doi:10.31219/osf.io/fzqxv](https://doi.org/10.31219/osf.io/fzqxv)
  Introduces the `mcp` package, formula syntax, default priors, and
  workflow for regression with multiple change points across generalized
  linear and time-series models. Newer-than-2020 versions of the paper
  may be available at that link. Please cite the newest version.

- Carlin, B. P., Gelfand, A. E., & Smith, A. F. (1992). Hierarchical
  Bayesian Analysis of Changepoint Problems. *Applied Statistics*,
  41(2), 389–405. [doi:10.2307/2347570](https://doi.org/10.2307/2347570)
  Introduced the Gibbs sampling approach for continuous change points in
  segmented regression and hierarchical models, providing the
  computational foundation used by BUGS and JAGS.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# \donttest{
# Define the segments using formulas. A change point is estimated between each formula.
model = list(
  response ~ 1,  # Plateau in the first segment (Intercept_1)
  ~ 0 + time,    # Joined slope (time_2) in segment 2 which starts at cp_1
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
#>     variable  mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        31.48 1.744 28.11 34.86 1.01      832     1686 30.0    OK
#>  cp_2        71.13 1.014 69.46 72.78 1.00     5341     7285 70.0    OK
#> 
#> Population-level parameters:
#>     variable  mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 10.03 0.664  8.73 11.32 1.00     1704     2937 10.0    OK
#>  time_2       0.53 0.047  0.44  0.63 1.00     1481     2737  0.5    OK
#>  Intercept_3 17.48 1.224 15.25 20.00 1.00     1196     1403 20.0    OK
#>  time_3      -0.10 0.072 -0.26  0.02 1.00     1255     1492 -0.3      
#>  sigma_1      3.91 0.290  3.39  4.53 1.00     4427     3954  3.5    OK

# Visual inspection of the results
plot(demo_fit)  # Visualization of model fit/predictions

plot_pars(demo_fit)  # Parameter distributions


pp_check(demo_fit)  # Prior/Posterior predictive checks


# Test a hypothesis
hypothesis(demo_fit, "cp_1 > 10")
#>      hypothesis     mean    lower    upper prob  BF
#> 1 cp_1 - 10 > 0 21.47938 18.10775 24.86071    1 Inf

# Make predictions
head(fitted(demo_fit))
#>   response     time   fitted        sd     Q2.5    Q97.5
#> 1 17.23552 76.33986 16.95304 0.9425042 15.17680 18.84139
#> 2 11.35171 83.51711 16.22386 0.7231217 14.82569 17.64630
#> 3 28.04995 60.18529 25.26025 0.8946280 23.51420 26.98299
#> 4 20.68198 74.72964 17.11663 1.0216233 15.20838 19.19173
#> 5 21.21364 85.88256 15.98354 0.7208332 14.60845 17.42739
#> 6 22.32282 40.05069 14.54825 0.6297756 13.25108 15.71204
head(predict(demo_fit))
#>   response     time  predict       sd      Q2.5    Q97.5
#> 1 17.23552 76.33986 16.95592 4.001623  9.031754 24.85364
#> 2 11.35171 83.51711 16.24130 3.989762  8.401454 24.04522
#> 3 28.04995 60.18529 25.22341 4.010209 17.360815 33.14115
#> 4 20.68198 74.72964 17.16957 3.993998  9.157086 25.05408
#> 5 21.21364 85.88256 15.97045 4.000445  8.166262 23.80838
#> 6 22.32282 40.05069 14.57599 4.029706  6.757755 22.33908
head(predict(demo_fit, newdata = data.frame(time = c(55.545, 80, 132))))
#>      time  predict       sd       Q2.5    Q97.5
#> 1  55.545 22.85129 3.965658 14.9562088 30.61321
#> 2  80.000 16.61806 3.999897  8.7254549 24.42420
#> 3 132.000 11.24520 5.215439  0.7408375 21.30547

# Compare to a one-intercept-only model (no change points) with default prior
model_null = list(response ~ 1)
fit_null = mcp(model_null, data = data, par_x = "time")  # fit another model here
demo_loo = loo(demo_fit)
null_loo = loo(fit_null)
loo::loo_compare(demo_loo, null_loo)
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model1       0.0     0.0      NA                    
#>  model2     -54.3     7.8    1.00                    

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
#>     variable     mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  cp_1        33.26717 22.92  2.01 83.05 1.00     8739     8782 30.0    OK
#>  cp_2        65.58075 23.61 15.35 98.41 1.00     8892     8149 70.0    OK
#> 
#> Population-level parameters:
#>     variable     mean    sd lower upper rhat ess_bulk ess_tail  sim match
#>  Intercept_1 14.39840 11.82 -3.47 33.57 1.00     9092     8866 10.0    OK
#>  time_2      -0.00142  0.12 -0.19  0.19 1.00     8721     8667  0.5      
#>  Intercept_3 14.16118 10.50 -4.81 32.93 1.00     9348     8556 20.0    OK
#>  time_3      -0.00018  0.11 -0.20  0.19 1.00     8175     8743 -0.3      
#>  sigma_1      6.63362  7.23  0.19 24.88 1.00     8799     8859  3.5    OK
plot(demo_fit, prior = TRUE)


# Show all priors. Default priors are added where you don't provide any
prior_summary(demo_fit)
#> # A tibble: 7 × 5
#>   parameter   segment dpar  prior                                         bounds
#>   <chr>         <int> <chr> <chr>                                         <chr> 
#> 1 cp_1              2 cp    dirichlet(alpha = 1)                          [min(…
#> 2 cp_2              3 cp    dirichlet(alpha = 1)                          [cp_1…
#> 3 Intercept_1       1 mu    student_t(df = 3, location = 14.2, scale = 6) none  
#> 4 time_2            2 mu    student_t(df = 3, location = 0, scale = 0.06… none  
#> 5 Intercept_3       3 mu    student_t(df = 3, location = 14.2, scale = 6) none  
#> 6 time_3            3 mu    student_t(df = 3, location = 0, scale = 0.06… none  
#> 7 sigma_1           1 sigma student_t(df = 3, location = 0, scale = 6)    [0.00…

# Set priors and re-run
prior = list(
  Intercept_1 = 15,
  time_2 = "dt(0, 2, 1) T(0, )",  # t-dist slope. Truncated to positive.
  cp_2 = "dunif(cp_1, 80)",       # change point to segment 2 > cp_1 and < 80.
  Intercept_3 = "Intercept_1"     # Shared intercept between segment 1 and 3
)

fit3 = mcp(model, data = data, prior = prior, warmup = 2000, iter = 6000, seed = 42)
#> Warning: Some parameters may not have converged well:
#>   * rhat > 1.01 or ess_bulk < 400 or ess_tail < 400: Intercept_3
#> Inspect `summary(fit)` and `plot_pars(fit)`, and consider increasing `iter`/`warmup` or simplifying the model before trusting these results.

# Show the JAGS model
demo_fit$jags_code
#> model {
#>   # mcp helper values
#>   cp_0 = CONST1_
#>   cp_3 = CONST2_
#> 
#>   # Priors for population-level effects
#>   cp_frac_1_ ~ dbeta(1, 2)  # Relative fraction of remaining span (Uniform order statistics)
#>   cp_1 = cp_0 + cp_frac_1_ * (cp_3 - cp_0)  # Ordered change point
#>   cp_frac_2_ ~ dbeta(1, 1)  # Relative fraction of remaining span (Uniform order statistics)
#>   cp_2 = cp_1 + cp_frac_2_ * (cp_3 - cp_1)  # Ordered change point
#>   Intercept_1 ~ dt(14.2, 1/(6)^2, 3)   # Robustly centered mean intercept with a minimum scale of 2.5
#>   time_2 ~ dt(0, 1/(0.06053105)^2, 3)   # Regularizing mean coefficient scaled to a reference predictor change
#>   Intercept_3 ~ dt(14.2, 1/(6)^2, 3)   # Robustly centered mean intercept with a minimum scale of 2.5
#>   time_3 ~ dt(0, 1/(0.06053105)^2, 3)   # Regularizing mean coefficient scaled to a reference predictor change
#>   sigma_1 ~ dt(0, 1/(6)^2, 3) T(0.001,)  # Positive residual SD calibrated on the response scale
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
