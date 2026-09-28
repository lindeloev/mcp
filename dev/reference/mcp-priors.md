# Priors in mcp

The `prior` argument of
[`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) is a
named list. Names are parameter names (`cp_1`, `Intercept_1`, `x_2`,
`sigma_1`, etc.; see
[mcp-formula](https://lindeloev.github.io/mcp/dev/reference/mcp-formula.md)
for how they arise) and the values specify their priors. Use
[`prior_summary()`](https://lindeloev.github.io/mcp/dev/reference/prior_summary.md)
to see the resolved priors of a fitted model, including defaults.

## Details

Each value in the list is either

- A distribution in `mcp`'s JAGS-string syntax (e.g.,
  `Intercept_1 = "dnorm(0, 1) T(0,)"`) indicating a conventional prior
  distribution. Data-calibrated, regularizing defaults are used where
  priors are not specified. These are designed for stable estimation and
  prediction, but should be justified before hypothesis testing. `mcp`
  uses conventional distribution scales rather than JAGS precision: SD
  for [`dnorm()`](https://rdrr.io/r/stats/Normal.html), scale for
  [`dt()`](https://rdrr.io/r/stats/TDist.html), `ddexp()`, and
  [`dlogis()`](https://rdrr.io/r/stats/Logistic.html), and log-SD for
  [`dlnorm()`](https://rdrr.io/r/stats/Lognormal.html).

- A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed
  value.

- A model parameter name (e.g., `Intercept_2 = "Intercept_1"`), equating
  parameters via JAGS. Note that sharing coefficients via `same()` in
  the segment formulas (e.g., `~ same(1)`) is generally preferred. If
  two group-level deviations are shared via the prior, they will need to
  have the same grouping variable.

Default coefficient priors follow the autoscaled defaults of `rstanarm`:
`normal(0, 2.5 * sd(y) / sd(x))` for identity links and
`normal(0, 2.5 / sd(x))` for log links, where `sd(x)` is the standard
deviation of the predictor's model-matrix column (including dummy
variables for factors). Logit and probit links use a narrower base
scale, `normal(0, 1.5 / sd(x))`. For the change-point variable, `sd(x)`
is taken over all data, and powers such as `I(x^2)` use
`sd((x - min(x))^2)`. Intercepts also follow `rstanarm`, e.g.,
`normal(mean(y), 2.5 * sd(y))` for Gaussian models (`normal(0, 1.5)` for
logit and probit), while `sigma` and group-level SDs follow `brms`,
e.g., `student_t(3, 0, max(2.5, mad(y)))`. These are scale parameters,
for which a half-Student-t is the standard weakly informative choice
(Gelman, 2006).

Like `rstanarm` and `brms`, which place intercept priors on centered
predictors, the default prior on the segment-1 intercept applies to the
level at the start of the segment, `min(x)`. The reported `Intercept_1`
is still at `x = 0`, as in [`lm()`](https://rdrr.io/r/stats/lm.html).
Later segments' intercepts are at their change point. User-specified
intercept priors apply to `Intercept_1` as written. Notice that the
autoscaling does not factor in the number of segments or their expected
x-range. Consider scaling by expected segment length on a case-by-case
basis if you have such a priori expectations.

The default prior on change points is `dirichlet(1)` (uniform order
statistics). For a single change point, this is the Beta(1, 1) / Uniform
distribution over `[min(x), max(x)]`. For multiple change points, it
corresponds to a flat Dirichlet distribution over segment lengths. You
can also explicitly set `cp_i = "dirichlet(alpha)"` with the same
positive `alpha` for all change points to regularize spacing
(`alpha > 1` penalizes change points from occurring close together,
while `alpha < 1` favors clustering). Under the hood, this is
parameterized as an exact sequential stick-breaking Beta chain for fast
and robust sampling. [Read
more](https://lindeloev.github.io/mcp/articles/priors.html).

## Notes on priors

- *Ordered change point priors:* Default population-level `cp_i` priors
  are ordered and the ordering is imposed through the priors. For
  user-defined priors, `mcp` adds truncation (e.g., `T(cp_1, )`) only
  when the prior has neither explicit truncation nor an inherently
  bounded form such as [`dunif()`](https://rdrr.io/r/stats/Uniform.html)
  or `dirichlet()`.

- *Data-dependent terms:* If `mcp` encounters a data-dependent term like
  `min(time)`, `max(time)`, `mean(response)`, `median(response)`,
  `sd(response)`, or `mad(response)` in the prior string, they are
  resolved from the model data so a numerical value is passed to JAGS.
  The following terms are also allowed: `n_segments()` and `n_cp()`. The
  older constants `MINX`, `MAXX`, `MEANX`, `SDX`, `MINY`, `MAXY`,
  `MEANY`, `SDY`, and `N_CP` remain accepted with a deprecation warning.

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

## See also

[`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md),
[`prior_summary()`](https://lindeloev.github.io/mcp/dev/reference/prior_summary.md),
[mcp-formula](https://lindeloev.github.io/mcp/dev/reference/mcp-formula.md)

## Examples

``` r
model = list(
  response ~ 1,
  ~ 0 + time,
  ~ 1 + time
)
prior = list(
  cp_1 = "dunif(20, 50)",         # Distribution
  time_2 = "dnorm(0, 2) T(0, )",  # Truncated distribution
  Intercept_3 = 15,               # Fixed value
  time_3 = "time_2"               # Shared with another parameter
)
fit = mcp(model, data = mcp_example_data("demo"), prior = prior, sample = "none")
prior_summary(fit)
#> # A tibble: 7 × 5
#>   parameter   segment dpar  prior                                      bounds   
#>   <chr>         <int> <chr> <chr>                                      <chr>    
#> 1 cp_1              2 cp    uniform(min = 20, max = 50)                [20, 50] 
#> 2 cp_2              3 cp    uniform(min = cp_1, max = 99.74132)        [cp_1, m…
#> 3 Intercept_1       1 mu    normal(mean = 15.32512, sd = 16.86313)     none     
#> 4 time_2            2 mu    normal(mean = 0, sd = 2)                   [0, Inf] 
#> 5 Intercept_3       3 mu    15                                         none     
#> 6 time_3            3 mu    time_2                                     none     
#> 7 sigma_1           1 sigma student_t(df = 3, location = 0, scale = 6) [0.001, …
```
