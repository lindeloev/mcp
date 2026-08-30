# Group-level (random) effects in mcp

`mcp` supports group-level effects (also called random effects) in both
the predictor and change-point parts of a formula. The syntax follows
`lme4` and `brms`: `(1|group)` specifies a group-level intercept,
`(factor||group)` specifies independent group coefficients,
[`fixef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
reports population-level effects, and
[`ranef()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
reports group-level deviations.

This article in brief:

- Predictor and change-point group-level effects
- How group-level effects persist across segments
- How to simulate group-level change-point deviations
- Get posteriors using `ranef(fit)`
- Plot using `plot(fit, facet_by="my_group")` and
  `plot_pars(fit, pars = "group", type = "dens_overlay", ncol = 3)`.
- The default priors bound group-specific change points relative to
  adjacent population-level change points.
- The article on modeling standard deviation via
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) contains [another
  group-level change-point
  example](https://lindeloev.github.io/mcp/dev/articles/dpar.md).

``` r

library(mcp)
future::plan(future::multisession, workers = 3)
set.seed(42)  # Make the script deterministic
```

## Specifying group-level change points

You specify group-level effects using the familiar
[`lmer`](https://www.rdocumentation.org/packages/lme4/versions/1.1-21/topics/lmer)
and `brms` syntax `(1|group)`. In the cp part, this models a
group-specific deviation from the population-level change point. For
example:

``` r

model = list(
  y ~ 1,  # Intercept_1
  1 + (1|id) ~ 0 + x  # cp_1, cp_1_sd, cp_1_id[i]
)
```

You can have multiple group-level change points, but they must use the
same grouping factor so their realized locations can be ordered within
each group:

``` r

model = list(
  y ~ 1,  # Intercept_1
  1 + (1|id) ~ 0 + x,  # cp_1, cp_1_sd, cp_1_id[i]
  1 + (1|id) ~ 0,      # cp_2, cp_2_sd, cp_2_id[i]
  (1|id) ~ 1           # cp_3 (implicit), cp_3_sd, cp_3_id[i]
)
```

Here are some properties of group-level change-point deviations:

**Exactly zero centered:** The group-level deviations sum to zero around
the associated population-level change point. Besides making `cp_i` the
arithmetic mean of the group-specific change points, this avoids a
poorly identified translation between `cp_i` and the average group
deviation.

**Hierarchical:** Consider the population-level change point `cp_1` and
its associated group-level deviations `cp_1_id`. Their spread is
governed by the population-level parameter `cp_1_sd`.

**Constraint:** Realized group-specific change points must remain
strictly ordered within each group.

Unlike predictor group-level effects, a change-point group effect
applies only to the boundary where it is written. There is no carry-over
as there is for group-level effects on RHS.

## Predictor group-level effects

The same syntax in the predictor part gives each group a deviation from
the population-level intercept:

``` r

model = list(
  y ~ 1 + (1|id),
  ~ 0 + x,
  ~ 1 + (0|id)
)
```

Here, the group-level intercept introduced in segment 1 also applies in
segment 2. A later `(1|id)` would replace it from that segment onward,
while `(0|id)` turns it off, as in segment 3 above. Persistence is
tracked separately for each grouping factor and distributional
parameter.

Use `||` for independent group-level slopes and factor coefficients:

``` r

model = list(
  y ~ 1 + condition + (condition||id),
  ~ 0 + x
)
```

With the default treatment coding, `(condition||id)` contains a
group-level intercept and one group-level deviation for each
non-reference contrast of `condition`. `(0 + condition||id)` instead
gives each factor level its own group-level coefficient. Similarly,
`(1 + z||id)` specifies independent group-level intercepts and slopes on
a numeric predictor `z`. Each coefficient has its own population-level
SD.

The population-level and group-level formulas need not contain the same
coefficients. For example, `y ~ 1 + (0 + condition||id)` has a
population-level intercept but only group-specific condition
coefficients.

Group-level intercepts also work inside distributional formulas:

``` r

model = list(
  y ~ 1 + (1|id) + sigma(1 + (1|id))
)
```

This model has group-level deviations in both the conditional mean and
log-SD. The `||` syntax works inside distributional formulas too, for
example `sigma(1 + (condition||id))`.

A later predictor group-level term for the same distributional parameter
and grouping factor replaces the entire earlier block. Thus, if
`(condition||id)` is followed in a later segment by `(1||id)`, only the
new intercept deviation applies from that segment onward. `(0|id)` turns
the block off.

Multi-coefficient terms with `|`, such as `(1 + x|id)`, would imply
correlated group coefficients and are not yet supported; use
`(1 + x||id)` for independent coefficients. Group-level terms inside
[`ar()`](https://rdrr.io/r/stats/ar.html) and `ma()` are also not
supported.

Unlike change-point deviations, predictor deviations are not constrained
to sum exactly to zero. They have an ordinary mean-zero hierarchical
normal distribution, matching multilevel regression in `lme4` and
`brms`. For example, `(1|id)` in segment 1 creates the deviation vector
`Intercept_1_id` and its population-level SD parameter
`Intercept_1_id_sd`. `(condition||id)` additionally creates names such
as `conditionB_1_id` and `conditionB_1_id_sd`; inside
[`sigma()`](https://rdrr.io/r/stats/sigma.html), the corresponding names
start with `sigma`.

## Simulating group-level change points

Let us simulate group-specific deviations in the change point between a
plateau and a slope:

``` r

model = list(
  y ~ 1,  # Intercept_1
  1 + (1|id) ~ 0 + x  # cp_1, cp_1_sd, cp_1_id[i]
)
```

This follows the same interface as predictor group-level effects: supply
the population-level SD and `fit$simulate()` draws one deviation per
group and centers them exactly on zero.

``` r

library(dplyr, warn.conflicts = FALSE)
group_levels = c("Clark", "Louis", "Batman", "Batgirl", "Spiderman", "Jane")
df = data.frame(
  x = runif(length(group_levels) * 30, 0, 100),  # 30 data points for each
  id = rep(group_levels, each = 30),  # the group names
  y = 1
)
empty = mcp(model, data = df, sample = FALSE)

set.seed(42)
df$y = empty$simulate(empty, df,
  # Population-level:
  Intercept_1 = 20, x_2 = 0.5, cp_1 = 50, sigma = 2,
  
  # Generate exactly centered group-level deviations with this SD
  cp_1_sd = 15)

head(df)
```

    ##          x    id        y
    ## 1 91.48060 Clark 36.10665
    ## 2 93.70754 Clark 34.00776
    ## 3 28.61395 Clark 24.03685
    ## 4 83.04476 Clark 28.74026
    ## 5 64.17455 Clark 22.60974
    ## 6 51.90959 Clark 24.57329

For models with multiple group-level change points, all deviations are
drawn once and then checked for ordering. An unordered draw produces an
error; `fit$simulate()` does not redraw, sort, or clip the change
points. The result:

``` r

library(ggplot2)
ggplot(df, aes(x=x, y=y)) + 
  geom_point() +
  facet_wrap(~id)
```

![](group_effects_files/figure-html/unnamed-chunk-9-1.png)

## Summarise and plot group-level effects

Fitting the model is simple:

``` r

fit = mcp(model, data = df, seed = 42)
```

If we just use `plot(fit)`, we would see all points in one plot. We want
to facet by `id`, so:

``` r

set.seed(42)
plot(fit, facet_by = "id")
```

![](group_effects_files/figure-html/unnamed-chunk-11-1.png)

It seems that `mcp` recovered the group-specific change points well.
There is a lot of information in these data because the population-level
intercept and slopes on either side of the change point are shared
across participants (`id`).

`summary(fit)` (or `fixef(fit)`) returns posterior summaries for the
population-level effects. To get the group-level deviations (random
effects), use:

``` r

mcp::ranef(fit)
```

    ##                 name        mean        sd       lower       upper     rhat
    ## 1   cp_1_id[Batgirl]   3.4359904 0.8359116   1.8387702   5.1155898 1.000170
    ## 2    cp_1_id[Batman]  -1.0441880 0.9203563  -2.8374458   0.7860453 1.000344
    ## 3     cp_1_id[Clark]  17.0539585 0.9404897  15.2304185  18.9432569 1.000462
    ## 4      cp_1_id[Jane]  -6.4355915 0.8460597  -8.0639491  -4.7646496 1.000115
    ## 5     cp_1_id[Louis] -13.8447712 0.7595938 -15.3688120 -12.3909721 1.000583
    ## 6 cp_1_id[Spiderman]   0.8346018 0.8040346  -0.7445893   2.4225902 1.000137
    ##   ess_bulk ess_tail         sim match
    ## 1     5460     5770   4.2419513    OK
    ## 2     5871     6189   0.1959384    OK
    ## 3     5216     5815  15.3133890    OK
    ## 4     5854     6000  -6.8428555    OK
    ## 5     4976     6457 -13.7214603    OK
    ## 6     6226     6974   0.8130371    OK

Inspecting the `sim` and `match` columns, we see that they recovered the
simulation parameters well.

Prediction methods include all group-level effects by default. Set
`varying = FALSE` for population-only predictions, `varying = "cp"` for
all group-level effects in the cp part, `varying = "predictor"` for
predictor-side group-level effects, or supply an exact name such as
`varying = "cp_1_id"`. `ranef(fit)` deliberately remains simple and
returns all group-level effects.

Good convergence is not always as obvious as in this example. While
`plot_pars(fit)` shows population-level parameters only, you can select
the group-level deviations with `"group"`:

``` r

set.seed(42)
plot_pars(fit, pars = "group", type = "trace", ncol = 3, nvariables = NULL)
```

![](group_effects_files/figure-html/unnamed-chunk-13-1.png)

The `ncol` argument controls the number of columns. Group-level effects
often have many levels, so this is useful for viewing all deviations.

Using `pars = "group"` plots all group-level deviations. To select one
group-level effect, use a regular expression in `regex_pars`; for
example, `^` anchors the start of a parameter name:

``` r

set.seed(42)
plot_pars(fit, regex_pars = "^cp_1_id", type = "dens_overlay", ncol = 2, nvariables = NULL)
```

![](group_effects_files/figure-html/unnamed-chunk-14-1.png)

You can also do posterior predictive checking with facets. I think that
for the relatively univariate models supported as of `mcp` 0.3, this
does not add much new information over and above
`plot(fit, facet_by = "id")`, but it’s a standard assessment that many
will be acquainted with:

``` r

set.seed(42)
pp_check(fit, facet_by = "id")
```

![](group_effects_files/figure-html/unnamed-chunk-15-1.png)

## Priors for group-level effects

You can see the priors of the model like this:

``` r

prior_summary(fit)
```

    ## # A tibble: 6 × 5
    ##   parameter   segment dpar  prior                                         bounds
    ##   <chr>         <int> <chr> <chr>                                         <chr> 
    ## 1 cp_1              2 cp    uniform(min = 0.02388966, max = 98.88917)     [min(…
    ## 2 cp_1_sd           2 cp    normal(mean = 0, sd = 197.7306)               [0, I…
    ## 3 cp_1_id           2 cp    normal(mean = 0, sd = cp_1_sd)                [min(…
    ## 4 Intercept_1       1 mu    student_t(df = 3, location = 23.4, scale = 7… none  
    ## 5 x_2               2 mu    student_t(df = 3, location = 0, scale = 0.07… none  
    ## 6 sigma_1           1 sigma student_t(df = 3, location = 0, scale = 7.2)  [0, I…

The prior `cp_1_sd` is the population-level standard deviation governing
the `cp_1_id` deviations across levels of `id`. Change-point deviations
are exactly zero-centered, while predictor group-level effects use
ordinary hierarchical mean-zero deviations. Their realized locations are
constrained to remain ordered within each group.

## JAGS code

Here is the JAGS code for the model used in this article:

``` r

fit$jags_code
```

    ## model {
    ##   # mcp helper values
    ##   cp_0 = CONST1_
    ##   cp_2 = CONST2_
    ## 
    ##   # Priors for population-level effects
    ##   cp_1 ~ dunif(CONST1_, CONST2_)  # Within the observed change-point span
    ##   cp_1_sd ~ dnorm(0, 1/(197.7306)^2) T(0,)  # Group-level change-point variation
    ##   Intercept_1 ~ dt(23.4, 1/(7.2)^2, 3)   # Robustly centered mean intercept with a minimum scale of 2.5
    ##   x_2 ~ dt(0, 1/(0.07282637)^2, 3)   # Regularizing mean coefficient scaled to a reference predictor change
    ##   sigma_1 ~ dt(0, 1/(7.2)^2, 3) T(0,)  # Positive residual SD calibrated on the response scale
    ## 
    ##   # Priors for group-level effects
    ##   for (id_ in 1:n_unique_id) {
    ##     cp_1_id_uncentered[id_] ~ dnorm(CONST3_, 1/(cp_1_sd)^2) T(CONST1_-cp_1,CONST2_-cp_1)  # Zero-mean ordered group-level change-point offsets
    ##   }
    ##   cp_1_id = cp_1_id_uncentered - mean(cp_1_id_uncentered)  # vectorized zero-centering
    ## 
    ##   # Model and likelihood
    ##   for (i_ in 1:length(x)) {
    ##     # par_x local to each segment
    ##     x_local_1_[i_] = min(x[i_], (cp_1 + cp_1_id[id[i_]]))
    ##     x_local_2_[i_] = min(x[i_], cp_2) - (cp_1 + cp_1_id[id[i_]])
    ##     
    ##     # Formula for mu
    ##     link_mu_[i_] =
    ##       (x[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(1)], c(Intercept_1)) * 1 + 
    ##       (x[i_] >= (cp_1 + cp_1_id[id[i_]])) * inprod(rhs_matrix_[i_, c(2)], c(x_2)) * x_local_2_[i_]
    ##     
    ##     # Formula for sigma
    ##     link_sigma_[i_] =
    ##       (x[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(3)], c(sigma_1)) * 1
    ## 
    ##     # Likelihood and log-density for family = gaussian()
    ##     mu_[i_] = link_mu_[i_]
    ##     sigma_[i_] = max(1e-03, link_sigma_[i_])
    ##     y[i_] ~ dnorm(mu_[i_], 1 / sigma_[i_]^2)
    ##   }
    ## }
