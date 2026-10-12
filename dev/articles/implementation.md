# From formulas to sampler code

The `mcp` model itself is described in [Understanding mcp
formulas](https://lindeloev.github.io/mcp/dev/articles/formulas.html#modeling-slope-change-points).
The actual implementation is equivalent, but quite different in form.
This article shows how `mcp` turns the segment formulas into code for
the sampler, currently JAGS. You can see it as R code in `fit$simulate`
and as JAGS code in `fit$jags_code`.

## Modeling changes in intercept

A change point is simply like an `ifelse` statement or multiplying with
indicators (0s and 1s):

``` r

# Model parameters
x = 1:20
cp_1 = 12
Intercept_1 = 5
Intercept_2 = 10

# Ifelse version
y_ifelse = ifelse(x <= cp_1, yes = Intercept_1, no = Intercept_2)

# Indicator equivalent using dummy helpers
cp_0 = -Inf
cp_2 = Inf
y_indicator = (x > cp_0) * (x <= cp_1) * Intercept_1 +  # Between cp_0 and cp_1
              (x > cp_1) * (x <= cp_2) * Intercept_2  # Between cp_1 and cp_2

# Show it
par(mfrow = c(1,2))
plot(x, y_ifelse, main = "ifelse(x <= cp_1)")
plot(x, y_indicator, main = "(x > cp_1) * Intercept_2")
```

![](implementation_files/figure-html/unnamed-chunk-2-1.png)

The magic of (Bayesian) MCMC sampling is that it can actually infer the
change point from this simple formulation. We let `mcp` write the JAGS
code for this simple two-plateaus model and see how it uses the
indicator formulation of change points:

``` r

library(mcp)
model = list(y ~ 1, ~ 1)
fit = mcp(model, data = data.frame(y = 1:10, x = 1:10), sample = FALSE, par_x = "x")
fit$jags_code
```

    ## model {
    ##   # Change points
    ##   cp_0 = min(x)
    ##   for (cp_i_ in 1:1) {
    ##     cp_unit_[cp_i_] ~ dunif(0, 1)
    ##   }
    ##   cp_sorted_ = sort(cp_unit_)  # Uniform order statistics
    ##   cp_1 = cp_0 + cp_sorted_[1] * (cp_2 - cp_0)
    ##   cp_2 = max(x)
    ## 
    ##   # Priors for population-level effects
    ##   Intercept_1 ~ dnorm(5.5, 1/(7.5691258852437295)^2)   # Mean intercept (rstanarm default)
    ##   Intercept_2 ~ dnorm(5.5, 1/(7.5691258852437295)^2)   # Mean intercept (rstanarm default)
    ##   sigma_1 ~ dt(0, 1/(3.7)^2, 3) T(0,)  # Positive residual SD calibrated on the response scale
    ## 
    ##   # Model and likelihood
    ##   for (i_ in 1:length(x)) {
    ##     # par_x local to each segment
    ##     x_local_1_[i_] = min(x[i_], cp_1)
    ##     x_local_2_[i_] = min(x[i_], cp_2) - cp_1
    ##     
    ##     # Formula for mu
    ##     link_mu_[i_] =
    ##       (x[i_] >= cp_0) * (x[i_] < cp_1) * inprod(rhs_matrix_[i_, c(1)], c(Intercept_1)) * 1 + 
    ##       (x[i_] >= cp_1) * inprod(rhs_matrix_[i_, c(2)], c(Intercept_2)) * 1
    ##     
    ##     # Formula for sigma
    ##     link_sigma_[i_] =
    ##       (x[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(3)], c(sigma_1)) * 1
    ## 
    ##     # Likelihood for family = gaussian()
    ##     mu_[i_] = link_mu_[i_]
    ##     sigma_[i_] = max(1e-05, link_sigma_[i_])
    ##     y[i_] ~ dnorm(mu_[i_], 1 / sigma_[i_]^2)
    ##   }
    ## }

Look at the section under the comment `# Formula for mu` which is the
(automatically generated) model that was discussed above. Some
unnecessary stuff is added to segment 1 just because it makes the code
easier to generate. (`x[i_] >= cp_0` when `cp_0` is the smallest value
of `x` is, of course, always true).

The section under `# Change points` implements the default
`dirichlet(1)` prior on the change points, which is the distribution of
sorted uniform draws (uniform order statistics). `mcp` samples uniform
draws on the unit interval, `cp_unit_`, sorts them, and scales them to
the observed range of x. Moving one change point therefore leaves the
others in place. The unit interval suits JAGS’s slice sampler, whose
step size starts at 1. Other `dirichlet(alpha)` priors are sampled as a
sequential stick-breaking Beta chain (`cp_1_frac_` etc.).

## Modeling changes in slope

The [mathematical
definition](https://lindeloev.github.io/mcp/dev/articles/formulas.html#modeling-slope-change-points)
of the piecewise model, with segments k \in \\1, \dots, K\\, change
points \tau_k, intercepts \alpha_k, and slopes \beta\_{k,1}, is in
Understanding mcp formulas. Here is how `mcp` writes it in JAGS.

### How mcp implements slope changes in JAGS using indicators

Bayesian MCMC samplers like JAGS cannot evaluate dynamic conditional
statements (`if/else`) on parameters being sampled. To implement this
piecewise model in JAGS, `mcp` uses an **indicator and plateauing
formulation** (evaluated via Iverson brackets \[A\] which equal 1 when
true and 0 when false).

First, `mcp` defines a span X\_{k,i} that measures x from the change
point within segment k and stays constant beyond it:

X\_{1,i} = \min(x_i, \tau_1) X\_{k,i} = \max\bigl(0, \min(x_i, \tau_k) -
\tau\_{k-1}\bigr) \quad (k \ge 2)

When x_i \< \tau\_{k-1}, X\_{k,i} = 0 (before the segment). When
\tau\_{k-1} \le x_i \< \tau_k, X\_{k,i} = x_i - \tau\_{k-1} (linear
growth within the segment). When x_i \ge \tau_k, X\_{k,i} = \tau_k -
\tau\_{k-1} (constant segment width).

For a **purely joined model** (`list(y ~ 0 + x, ~ 0 + x)`), `mcp` writes
\eta_i as a cumulative sum of these segment spans:

\eta_i = \sum\_{k=1}^K \[x_i \ge \tau\_{k-1}\] \cdot \beta\_{k,1}
X\_{k,i}

Because each segment’s contribution stays constant beyond its change
point (\beta\_{k-1,1}(\tau\_{k-1} - \tau\_{k-2})), joining needs no
parameter constraints: earlier segments provide the starting value for
later segments. The expected response on the response scale is \mu_i =
g^{-1}(\eta_i) via link function g(\mu_i) = \eta_i (identity for
Gaussian, so \mu_i = \eta_i).

When a segment is **disjoined** (`~ 1 + x`), it introduces a new
intercept \beta\_{k,0} at \tau\_{k-1}, and all earlier terms get the
indicator \[x_i \< \tau\_{k-1}\] so they end where the disjoined segment
starts. Terms that are not x-terms get the same indicator at the next
change point, which is why they end there.

Let’s demonstrate this equivalence in R:

``` r

# Model parameters
x = 1:20
cp_1 = 12
slope_1 = 2
slope_2 = -1

# Ifelse version
y_ifelse = ifelse(x <= cp_1, 
            yes = slope_1 * x,
            no = cp_1 * slope_1 + slope_2 * (x - cp_1))

# Indicator version with local plateauing coordinates (pmin is vectorized min)
y_local = slope_1 * pmin(x, cp_1) + 
          (x > cp_1) * slope_2 * (x - cp_1)

# Show that both formulations are identical
par(mfrow = c(1,2))
plot(x, y_ifelse, main = "ifelse() version")
plot(x, y_local, main = "Local plateauing coordinate")
```

![](implementation_files/figure-html/unnamed-chunk-4-1.png)

Let us see this in action:

``` r

model = list(y ~ 0 + x, ~ 0 + x)
fit = mcp(model, data = data.frame(y = 1:10, x = 1:10), sample = FALSE, reparam = FALSE)
fit$jags_code
```

    ## model {
    ##   # Change points
    ##   cp_0 = min(x)
    ##   for (cp_i_ in 1:1) {
    ##     cp_unit_[cp_i_] ~ dunif(0, 1)
    ##   }
    ##   cp_sorted_ = sort(cp_unit_)  # Uniform order statistics
    ##   cp_1 = cp_0 + cp_sorted_[1] * (cp_2 - cp_0)
    ##   cp_2 = max(x)
    ## 
    ##   # Priors for population-level effects
    ##   x_1 ~ dnorm(0, 1/(2.5)^2)   # Autoscaled mean coefficient (rstanarm default)
    ##   x_2 ~ dnorm(0, 1/(2.5)^2)   # Autoscaled mean coefficient (rstanarm default)
    ##   sigma_1 ~ dt(0, 1/(3.7)^2, 3) T(0,)  # Positive residual SD calibrated on the response scale
    ## 
    ##   # Model and likelihood
    ##   for (i_ in 1:length(x)) {
    ##     # par_x local to each segment
    ##     x_local_1_[i_] = min(x[i_], cp_1)
    ##     x_local_2_[i_] = min(x[i_], cp_2) - cp_1
    ##     
    ##     # Formula for mu
    ##     link_mu_[i_] =
    ##       (x[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(1)], c(x_1)) * x_local_1_[i_] + 
    ##       (x[i_] >= cp_1) * inprod(rhs_matrix_[i_, c(2)], c(x_2)) * x_local_2_[i_]
    ##     
    ##     # Formula for sigma
    ##     link_sigma_[i_] =
    ##       (x[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(3)], c(sigma_1)) * 1
    ## 
    ##     # Likelihood for family = gaussian()
    ##     mu_[i_] = link_mu_[i_]
    ##     sigma_[i_] = max(1e-05, link_sigma_[i_])
    ##     y[i_] ~ dnorm(mu_[i_], 1 / sigma_[i_]^2)
    ##   }
    ## }

We used `reparam = FALSE` to see the plain model. By default, `mcp`
samples some parameters on other scales, which we will return to
[below](#reparameterizations).

Look at the JAGS code under `# par_x local to each segment` and
`# Formula for mu` to see this exact indicator formulation in action:
`x_local_1_[i_]` and `x_local_2_[i_]` correspond to X\_{1,i} and
X\_{2,i}, and each segment’s slope is multiplied by its activation
indicator (e.g., `(x[i_] >= cp_1)`).

You will find the exact same formula for `y_ = ...` if you do
`print(fit$simulate)`, though this function contains a whole lot of
other stuff too.

## Reparameterizations

By default (`mcp(..., reparam = TRUE)`), `fit$jags_code` samples some
coefficients via auxiliary nodes, ending in `_mid_`, `_end_`, or
`_rise_`. This is a reparameterization, not a new model: the priors, and
hence the posteriors, are exactly those of the plain model. Only the
sampling becomes faster and converges more reliably. `reparam = FALSE`
shows the plain model, except that the segment-1 intercept is still
sampled at the midpoint of segment 1 because sampling is otherwise very
slow when x is far from zero. Here are excerpts of `fit$jags_code` for
`list(y ~ 1 + x, ~ 0 + x, ~ 1 + x)` with `x = 1:100`, first with
`reparam = FALSE`:

``` r

Intercept_1_mid_ ~ dnorm(0 - x_1 * 1 + x_1 * (17.5), 1/(2.5)^2)  # Prior on the value at min(x) = 1
Intercept_1 = Intercept_1_mid_ - (x_1 * (17.5))  # Also with reparam = FALSE
x_1 ~ dnorm(0, 1/(0.087)^2)
x_2 ~ dnorm(0, 1/(0.087)^2)
Intercept_3 ~ dnorm(0, 1/(2.5)^2)
x_3 ~ dnorm(0, 1/(0.087)^2)
```

And with `reparam = TRUE`:

``` r

# Segment 1: intercept sampled as the value at x = 17.5
Intercept_1_mid_ ~ dnorm(0 - x_1 * 1 + x_1 * (17.5), 1/(2.5)^2)
Intercept_1 = Intercept_1_mid_ - (x_1 * (17.5))
x_1 ~ dnorm(0, 1/(0.087)^2)

# Joined segment: slope sampled as the value at the segment end
x_2_end_ ~ dnorm(Intercept_1 + x_1 * cp_1, 1/(0.087*(cp_2-cp_1))^2)
x_2 = (x_2_end_ - (Intercept_1 + x_1 * cp_1)) / (cp_2 - cp_1)

# Disjoined segment: intercept sampled as the value at the segment end
Intercept_3_end_ ~ dnorm(0 + x_3 * (cp_3 - cp_2), 1/(2.5)^2)
Intercept_3 = Intercept_3_end_ - (x_3 * (cp_3 - cp_2))
x_3 ~ dnorm(0, 1/(0.087)^2)
```

### Why this is easier for the sampler

Like most Gibbs samplers, JAGS updates one parameter at a time while
holding the others fixed. In the plain model, an intercept is the value
at the start of its segment, and a joined segment starts where the
previous segment ended. So if JAGS moves a change point \tau_k while
holding the slopes fixed, the lines of segment k+1, and of all later
joined segments, shift up or down. Such moves fit the data badly and are
mostly rejected, so the change points and the coefficients must creep
along together, which mixes slowly.

With `reparam = TRUE`, `mcp` instead samples quantities that stay in
place when a change point moves:

- **Joined segments:** The value at the segment end, e_k. Moving
  \tau\_{k-1} now pivots segment k around its end, so later segments
  stay in place.
- **Disjoined segments:** The value at the segment end, e_k. Then
  segment k’s line, e_k + \beta\_{k,1}(x - \tau_k), does not depend on
  \tau\_{k-1} at all.
- **Segment 1:** The value at the fixed point c = \min(x) + (\max(x) -
  \min(x)) / (2K), which is the midpoint of segment 1 under the
  `dirichlet(1)` prior. Like centering x in ordinary regression, this
  makes the intercept and slope less correlated.

Further x-terms in joined segments (e.g., `I(x^2)`) are sampled as their
rise over the segment, \beta w_k^2. The same rules apply to
distributional parameters such as `sigma(1 + x)`.

### Why this is the same model

All of the above follow from one fact: if a coefficient \beta has a
location-scale prior, e.g., \beta \sim \text{Normal}(m, s), then for any
a and w \> 0,

u = a + w \beta \quad \sim \quad \text{Normal}(a + w m, \\ w s).

So JAGS can sample u with this prior and derive \beta = (u - a) / w, and
\beta has exactly its original prior. This holds when a and w depend on
other parameters, such as change points and other coefficients, because
the prior on \beta given those parameters is unchanged. Here is what
`mcp` samples, with segment widths w_k = \tau_k - \tau\_{k-1} and
v\_{k-1} being the value where segment k starts:

| Coefficient | Sampled u | a | w |
|----|----|----|----|
| Segment 1 intercept \beta\_{1,0} | Value at c | \beta\_{1,1} c | 1 |
| Joined slope \beta\_{k,1} | Value at segment end e_k = v\_{k-1} + \beta\_{k,1} w_k | v\_{k-1} | w_k |
| Joined `I(x^2)` term \beta\_{k,2} | Rise \beta\_{k,2} w_k^2 | 0 | w_k^2 |
| Disjoined intercept \beta\_{k,0} | Value at segment end e_k = \beta\_{k,0} + \beta\_{k,1} w_k | \beta\_{k,1} w_k | 1 |

Let us verify this for the joined slope, where the sampled value at the
segment end depends on both the starting value and the segment width:

``` r

s = 0.5  # Prior SD of the slope
start = rnorm(100000, 10, 3)  # Value where the segment starts (set by earlier segments)
width = runif(100000, 0.5, 50)  # Segment widths (set by the change points)
end = rnorm(100000, mean = start, sd = s * width)  # Sampled value at the segment end
slope = (end - start) / width  # Derived slope

c(mean = mean(slope), sd = sd(slope))  # Normal(0, 0.5), as for the plain slope
```

    ##        mean          sd 
    ## 0.001082462 0.501374861

``` r

c(cor(slope, start), cor(slope, width))  # ... regardless of start and width
```

    ## [1]  0.0044027564 -0.0006086204

`mcp` only reparameterizes its default
[`dnorm()`](https://rdrr.io/r/stats/Normal.html) and
[`dt()`](https://rdrr.io/r/stats/TDist.html) priors, and for the
segment-1 intercept also truncated and
[`dunif()`](https://rdrr.io/r/stats/Uniform.html) priors, since a shift
(w = 1) moves their bounds by a too. User-specified priors are sampled
as written, and so are coefficients in segments bounded by group-level
change points, where the segment widths differ between groups.
