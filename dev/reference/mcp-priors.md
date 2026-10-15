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

- A distribution string (e.g., `Intercept_1 = "dnorm(0, 1) T(0, )"`).
  Distributions use JAGS names and argument order, including
  [`T()`](https://rdrr.io/r/base/logical.html) truncation, but with
  conventional scales instead of JAGS precision: SD for
  [`dnorm()`](https://rdrr.io/r/stats/Normal.html), SD on the log scale
  for [`dlnorm()`](https://rdrr.io/r/stats/Lognormal.html) (`sdlog` in
  [`stats::dlnorm()`](https://rdrr.io/r/stats/Lognormal.html)), and
  scale for [`dt()`](https://rdrr.io/r/stats/TDist.html), `ddexp()`, and
  [`dlogis()`](https://rdrr.io/r/stats/Logistic.html). Change points
  also accept `dirichlet(alpha)`, and all priors accept data-dependent
  terms such as `sd(y)` (see below).

- A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed
  value.

- A model parameter name (e.g., `Intercept_2 = "Intercept_1"`), sharing
  that parameter's value. Note that sharing coefficients via `same()` in
  the segment formulas (e.g., `~ same(1)`) is generally preferred. If
  two group-level deviations are shared via the prior, they will need to
  have the same grouping variable.

Weakly informative defaults are used where priors are not specified.
These are designed for stable estimation and prediction, but should be
justified before hypothesis testing. Default priors follow `rstanarm`
for intercepts and coefficients and `brms` for scale parameters, with a
few documented deviations:

- *Coefficients* are autoscaled as in `rstanarm`. For the mean (`mu`),
  \$\$\beta \sim \mathrm{Normal}\left(0, \frac{2.5 \cdot
  s_y}{\mathrm{sd}(\mathrm{column})}\right)\$\$ where \\s_y =
  \mathrm{sd}(y)\\ for identity links and \\s_y = 1\\ for other links.
  \\\mathrm{sd}(\mathrm{column})\\ is the SD of the coefficient's
  model-matrix column, including the 0/1 dummies of factors. This makes
  the prior largely invariant to the units of `x` and `y`. Logit and
  probit links use 1.5 and 1 instead of 2.5. For x-terms, `sd(column)`
  is `sd(x)` over all data (powers such as `I(x^2)` use
  `sd((x - min(x))^2)`) times the SD of the term without `x`, e.g.,
  `sd(x) * sd(z)` for `x:z`. 0/1 indicators such as the dummies in
  `x:state` contribute a factor of 1. Unlike `rstanarm`, this does not
  depend on the location of `x` and `z`.

- *Intercepts* are `normal(mean(y), 2.5 * sd(y))` for identity links,
  `normal(0, 1.5)` for logit, and `normal(0, 1)` for probit. The latter
  two are (approximately) uniform on the probability scale. For log
  links (Gaussian, Poisson, and negative binomial), intercepts are
  `normal(median(log(y)), max(2.5, mad(log(y))))`, rounded to one
  decimal and with `y <= 0` replaced by 0.1 so that the log is defined.
  This robustly centers the prior on the typical log-response, as in
  `brms` (which uses a Student-t shape and replaces zeros the same way),
  whereas `rstanarm` uses `normal(0, 2.5)`. For a factor without an
  intercept, e.g., `y ~ 0 + state`, each factor level's coefficient is
  an intercept and gets this prior where the segment starts at zero
  (segment 1 or joined to nothing). With an
  [`offset()`](https://rdrr.io/r/stats/offset.html), the intercept
  location and scale are computed from the offset-adjusted response,
  `y - offset` for identity links and `log(y) - offset` for log links.

- *`sigma` and group-level SDs* follow `brms`, e.g.,
  `student_t(3, 0, max(2.5, mad(y)))`, since half-Student-t priors are a
  standard weakly informative choice for scale parameters (Gelman,
  2006). For `mu` with a log link and for modeled `shape`, group-level
  SDs are half-normal instead, e.g., `normal(0, 2.5) T(0, )`, because
  exponentiated heavy tails imply extreme group multipliers.

Like `rstanarm` and `brms`, which place intercept priors on centered
predictors, the default prior on the segment-1 intercept applies to the
value at the start of the segment, `min(x)`. The reported `Intercept_1`
is still at `x = 0`, as in [`lm()`](https://rdrr.io/r/stats/lm.html), so
its prior is given the x-terms, e.g.,
`dnorm(mean(y) - x_1 * min(x), 2.5 * sd(y))` in `fit$prior`. Later
segments' intercepts are at their change point. Likewise, factor levels
and contrasts with their own x-terms apply at `min(x)`, e.g., `stateb_1`
in `y ~ 1 + state * x` gets `dnorm(-stateb.x_1 * min(x), ...)`. Unlike
in `rstanarm` and `brms`, other predictors are not centered, so the
prior applies where they are zero. If such a predictor is far from zero,
center it before fitting or set a prior on `Intercept_1`. User-specified
intercept priors apply directly to `Intercept_1` (at `x = 0`). The
autoscaling does not account for segment length. If you expect steep
slopes or strong curvature within short segments, consider setting a
wider prior.

The default prior on change points is `dirichlet(1)` (uniform order
statistics). For a single change point, this is the Beta(1, 1) / Uniform
distribution over `[min(x), max(x)]`. For multiple change points, it
corresponds to a flat Dirichlet distribution over segment lengths. You
can also explicitly set `cp_i = "dirichlet(alpha)"` with the same
positive `alpha` for all change points to regularize spacing
(`alpha > 1` penalizes change points from occurring close together,
while `alpha < 1` favors clustering). Under the hood, `dirichlet(1)` is
sampled as sorted uniforms and other `alpha` as an exact sequential
stick-breaking Beta chain for fast and robust sampling. [Read
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
#>   parameter   segment dpar  prior                                        bounds 
#>   <chr>         <int> <chr> <chr>                                        <chr>  
#> 1 cp_1              2 NA    uniform(min = 20, max = 50)                  [20, 5…
#> 2 cp_2              3 NA    uniform(min = 0.61863, max = 99.741)         [0.618…
#> 3 Intercept_1       1 mu    normal(mean = 13.825, sd = 17.383)           none   
#> 4 time_2            2 mu    normal(mean = 0, sd = 2)                     [0, In…
#> 5 Intercept_3       3 mu    15                                           none   
#> 6 time_3            3 mu    time_2                                       none   
#> 7 sigma_1           1 sigma student_t(df = 3, location = 0, scale = 5.7) [0, In…
```
