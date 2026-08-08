# Time series change points

Serial dependence is common in time series. Add autoregressive terms
with `ar(p)`, moving-average terms with `ma(q)`, or both in the segment
formulas:

``` r

model = list(y ~ 1 + x + ar(1) + ma(1))
```

The most common use case is still `ar(1)`. Like other `mcp` terms, AR
and MA coefficients carry over to later segments until another
[`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()` term changes them.
Both accept a regression formula, such as `ar(1, 1 + x)` or
`ma(1, 0 + x)`.

### GARMA definition

`mcp` implements AR and MA as a generalized ARMA (GARMA) recurrence on
the response-family link scale. If $`b_t`$ is the ordinary regression
predictor and $`\eta_t`$ is the predictor including serial dependence,
then

``` math
\eta_t = b_t
+ \sum_{j=1}^{p} \phi_{j,t} \left[g(y^*_{t-j}) - b_{t-j}\right]
+ \sum_{k=1}^{q} \theta_{k,t} \left[g(y^*_{t-k}) - \eta_{t-k}\right].
```

Thus [`ar()`](https://rdrr.io/r/stats/ar.html) uses lagged link-scale
residuals relative to the ordinary regression, while `ma()` uses lagged
one-step innovations. Unavailable lags at the beginning of the series
contribute zero. This same recurrence is used by JAGS, fitted values,
predictions, log likelihoods, and fresh-series simulation.

`mcp` evaluates this as a finite, conditional recurrence: the observed
history supplies earlier residuals, and unavailable residuals before the
start of the series are set to zero.

The AR and MA coefficients are modeled directly and are not jointly
constrained to the stationary or invertible regions. With higher-order
terms such as `ar(2)`, independently bounded coefficients can violate
the usual root conditions. Coefficients that change across predictors or
segments, such as `ar(1, 1 + x)` or `model = list(y ~ ar(1), ~ ar(1))`,
instead define a time-varying process to which the usual
constant-coefficient conditions do not directly apply. The same applies
to `ma()`.

GARMA currently supports only the default links for
[`gaussian()`](https://rdrr.io/r/stats/family.html) (identity),
[`binomial()`](https://rdrr.io/r/stats/family.html) (logit),
[`poisson()`](https://rdrr.io/r/stats/family.html) (log), and
[`negbinomial()`](https://lindeloev.github.io/mcp/dev/reference/negbinomial.md)
(log).
[`bernoulli()`](https://lindeloev.github.io/mcp/dev/reference/bernoulli.md)
and non-default links are rejected for now.

#### Observation boundary

The transformed observation $`y^*`$ keeps log and logit residuals finite
when counts lie on a link boundary. The default `boundary = 0.1`
replaces zero counts by 0.1 for Poisson and negative-binomial models.
For binomial models, counts are constrained to the interval from 0.1 to
`trials - 0.1` before conversion to a rate. It has no effect for
Gaussian models.

The default should usually be left unchanged. If needed, set it on
either term, for example `ar(1, boundary = 0.01) + ma(1)`. AR and MA
share one boundary within a segment, and the boundary may differ between
segments. Supplying different AR and MA boundaries in the same segment
is an error.

## Simple example

Let’s simulate some data using the
[`mcp_example()`](https://lindeloev.github.io/mcp/dev/reference/mcp_example.md)
function with the “ar” settings:

``` r

library(mcp)
future::plan(future::multisession, workers = 3)
set.seed(42)  # Make the script deterministic

ex = mcp_example("ar", sample = FALSE)
```

    ## Generating residuals for AR(N) model since the response column/argument was not provided.

``` r

head(ex$data)
```

    ##      price time
    ## 1 26.85479    1
    ## 2 19.91843    2
    ## 3 22.81123    3
    ## 4 24.27657    4
    ## 5 24.15365    5
    ## 6 21.77232    6

See how this was generated in `ex$call`. We model this as a plateau
(`1`) with a second-order autoregressive residual (`ar(2)`) followed by
a joined slope (`0 + time`) with a negative first-order autoregressive
residual (`ar(1)`):

``` r

model = list(
  price ~ 1 + ar(2),  # Intercept_1, ar1_1, ar2_1
  ~ 0 + time + ar(1)  # time_2, ar1_2
)
fit = mcp(model, ex$data)
```

Let’s plot it and we see that AR was strong in the first segment and
weaker-but-negative in the second:

``` r

plot(fit)
```

![](arma_files/figure-html/unnamed-chunk-3-1.png)

We can summarise the inferred coefficients:

``` r

summary(fit)
```

    ## Family: gaussian(link = 'identity')
    ## Iterations: 3000 from 3 chains.
    ## Segments:
    ##   1: price ~ 1 + ar(2)
    ##   2: price ~ 1 ~ 0 + time + ar(1)
    ## 
    ## Population-level parameters:
    ##         name match   sim  mean  lower  upper Rhat ess_bulk ess_tail
    ##         cp_1    OK 75.00 75.00 68.053 83.461    1      458      691
    ##  Intercept_1    OK 20.00 21.19 18.453 24.667    1      633      801
    ##       time_2    OK  0.50  0.45  0.334  0.578    1      854      718
    ##      sigma_1    OK  5.00  5.32  4.679  6.109    1     3444     3331
    ##        ar1_1    OK  0.40  0.43  0.219  0.638    1     4123     6284
    ##        ar1_2    OK -0.30 -0.29 -0.599  0.032    1     2418     1542
    ##        ar2_1    OK  0.15  0.20  0.014  0.392    1     2443     2988
    ## 
    ## Warning: 1 parameter shows poor convergence (Rhat > 1.01 or ESS < 400).

The naming syntax is `[component][order]_[segment]` for intercepts. For
example, `ar1_2` is the first-order autoregressive coefficient and
`ma1_2` the first-order moving-average coefficient in segment 2. Slopes
use the usual `mcp` parameter names, e.g., `ar1_x_3` for a slope on
AR(1) in segment 3.

Comparing the columns `mean` and `sim` we see that the AR coefficients
are reasonably recovered. In fact, the posterior mean is almost always
exactly the same as `arima(data, order = c(N, 0, 0))` (see below), so
the non-perfect fits are due to randomness in the simulation - not in
the fit.

For Gaussian GARMA models, `sigma` describes the *innovations*, i.e.,
the part of the residuals not explained by AR and MA coefficients.
`sd(ex$data$price)` is therefore generally higher. In this case, the SD
of raw data in the plateau is 8.99. As always, it is good to assess
posteriors and convergence more directly:

``` r

plot_pars(fit)
```

![](arma_files/figure-html/unnamed-chunk-5-1.png)![](arma_files/figure-html/unnamed-chunk-5-2.png)

Sometimes, the trace plot shows that the change point (`cp_1`) is not
well identified with this model and data. As discussed in the article on
[tips, tricks, and
debugging](https://lindeloev.github.io/mcp/dev/articles/tips.md), you
could combine a more informative prior with more samples
(`mcp(..., adapt = 10000, iter = 10000)`), if this is a problem.

You can do hypothesis testing with GARMA models using
[`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md).
[Read more
here](https://lindeloev.github.io/mcp/dev/articles/comparison.md) or
scroll down for an applied example.

## Tips, comments, and warnings

The AR and MA terms apply to link-scale *residuals* from the ordinary
regression. In time-series jargon, this is a *dynamical* regression
model where the ordinary regression parameters make up the
*deterministic* structure. See further comments in the section on priors
below.

To control the *direction* of a change, put a one-sided prior on the
corresponding coefficient, e.g., `ar1_x_1 = "dnorm(0, 1) T(0, )"`. See
recommendations in [the section on priors](#ar_prior).

### Regression on AR/MA coefficients

You can specify how the coefficients change with $`x`$ using
`ar(p, formula)` and `ma(q, formula)`. For example, `ar(1, 1 + time)`
models a steady change in AR(1) strength. These regression formulas use
an identity link, so slopes can easily produce non-stationary or
non-invertible coefficients. See [the section on priors](#ar_prior) for
ways to constrain them. The same formula is applied separately to every
order in the term.

Population-level [`mcp` formula
syntax](https://lindeloev.github.io/mcp/dev/articles/formulas.md) is
available inside both terms, so you could use
`ar(3, 1 + I(x^2) + exp(z))` or `ma(1, 0 + group)`. Group-level terms
such as `(1|id)` are not supported inside
[`ar()`](https://rdrr.io/r/stats/ar.html) or `ma()`. Transformed
predictors must stay within the transformation’s domain; for example,
values passed to [`log()`](https://rdrr.io/r/base/Log.html) and
[`sqrt()`](https://rdrr.io/r/base/MathFun.html) must be non-negative.

### Combining ar(), ma(), and sigma()

You can combine [`ar()`](https://rdrr.io/r/stats/ar.html) and `ma()`
with any regression model and with [varying change
points](https://lindeloev.github.io/mcp/dev/articles/varying.md). For
Gaussian models, [`sigma()`](https://rdrr.io/r/stats/sigma.html)
controls the innovation standard deviation:

``` r

model = list(
  y ~ 1 + ar(1) + ma(1) + sigma(1),
  ~ 0 + x + ar(1) + sigma(1),
  ~ 1
)
```

### Order your data

AR and MA lags apply to the *order* of rows in the data frame without
taking into account the distance between values of `x`. This has two
important implications:

- You probably want to sort your data according to your `x`. Just do
  `data = data[order(data$x), ]`.
- Adjacent data points that lie years apart are modeled to be just as
  (auto)correlated as adjacent points lying seconds apart.

## Simulating autocorrelated change point data

Assessing the correctness of autocorrelation is less intuitive than
seeing e.g. a mean fit. However, we can verify `mcp` up against more
tested-and-tried functions such as
[`arima()`](https://rdrr.io/r/stats/arima.html) in base R. Let us
simulate a single AR(3) segment, i.e., without change points, and see if
it fits:

``` r

# Model
model = list(response ~ 1 + ar(3))

# Simulate data
df = data.frame(time = 1:200, response = 0)
empty = mcp(model, df, sample = FALSE, par_x = "time")
set.seed(42)  # For consistent "random" results
df$response = empty$simulate(
    empty,
    df,
    Intercept_1 = 20,
    ar1_1 = 0.7, 
    ar2_1 = 0.2, 
    ar3_1 = -0.4, 
    sigma_1 = 8)
```

    ## Generating residuals for AR(N) model since the response column/argument was not provided.

``` r

# Base arima AR(3)
arima(df$response, order = c(3, 0, 0))
```

    ## 
    ## Call:
    ## arima(x = df$response, order = c(3, 0, 0))
    ## 
    ## Coefficients:
    ##          ar1     ar2      ar3  intercept
    ##       0.6944  0.2297  -0.4078    19.5891
    ## s.e.  0.0643  0.0798   0.0650     1.1356
    ## 
    ## sigma^2 estimated as 60.17:  log likelihood = -694.1,  aic = 1398.19

OK, we can see that the `ar` coefficients and sigma
($`sigma = sqrt(sigma^2)`$) is simulated correctly, if taking
[`arima()`](https://rdrr.io/r/stats/arima.html) as ground truth.
Inferring with `mcp` is straightforward:

``` r

fit = mcp(model, df, par_x = "time")
```

The Bayesian parameter estimates are in perfect correspondence with
[`arima()`](https://rdrr.io/r/stats/arima.html), even where they deviate
a tiny bit from the simulation parameters (due to the inherent
randomness in simulating data):

``` r

fixef(fit)
```

    ##          name match  sim       mean       lower      upper     Rhat ess_bulk
    ## 1 Intercept_1    OK 20.0 19.6587422 17.41107992 21.9667288 1.000260     5345
    ## 2     sigma_1    OK  8.0  7.8992496  7.15026305  8.7461457 1.000606     4712
    ## 3       ar1_1    OK  0.7  0.6890901  0.56164757  0.8178585 1.000811     2790
    ## 4       ar2_1    OK  0.2  0.2297764  0.07191371  0.3885284 1.000854     2090
    ## 5       ar3_1    OK -0.4 -0.3990256 -0.52560089 -0.2703786 1.000546     3218
    ##   ess_tail
    ## 1     5199
    ## 2     4758
    ## 3     5160
    ## 4     3663
    ## 5     5449

## Inferring an autocorrelation-only change

One “side-effect” of the `mcp` implementation of autocorrelation using
[`ar()`](https://rdrr.io/r/stats/ar.html) *in the formulas* is that you
can infer when autocorrelation parameters and structures change.

Let’s simulate a change point in autocorrelation and see if we can infer
it later:

``` r

# The model
model = list(
  y ~ 1 + x + ar(1),  # Slope
  ~ 0 + x + ar(1)  # Slope
)

# Get predictions
df = data.frame(x = seq(0, 100, length.out = 200), y = 0)
empty = mcp(model, df, sample = FALSE)
set.seed(42)
df$y = empty$simulate(
  empty,
  df,
  cp_1 = 60,
  Intercept_1 = 20,
  x_1 = 1, x_2 = 1,  # same slope
  ar1_1 = 0.8, ar1_2 = 0.2,
  sigma_1 = 5)
```

… and we use a prior to equate the slopes of each segment (read more
about [using
priors](https://lindeloev.github.io/mcp/dev/articles/priors.md) to
equate parameters and define constants). Now let’s see if we can recover
these parameters. We use `sample = "both"` because we will do a
Savage-Dickey test later.

``` r

prior = list(x_2 = "x_1")  # Set the two slopes equal
fit = mcp(model, data = df, prior = prior, sample = "both")
```

First, let’s get a visual to see that the posterior is reasonably narrow
and consistent:

``` r

plot(fit)
```

![](arma_files/figure-html/unnamed-chunk-12-1.png)

You could use `plot(fit, geom_data = "line")` for a more classical line
plot of the time series data. Set `plot(fit, arma = FALSE)` to omit all
AR and MA effects. You can also plot `ar1_` or `ma1_` parameters
directly using
[`plot_dpar()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md):

``` r

plot_dpar(fit, dpar = "ar1", lines = 100)
```

![](arma_files/figure-html/unnamed-chunk-13-1.png)

We recovered the parameters, including the change point:

    ## Family: gaussian(link = 'identity')
    ## Iterations: 3000 from 3 chains.
    ## Segments:
    ##   1: y ~ 1 + x + ar(1)
    ##   2: y ~ 1 ~ 0 + x + ar(1)
    ## 
    ## Population-level parameters:
    ##         name match  sim   mean lower upper Rhat ess_bulk ess_tail
    ##         cp_1    OK 60.0 63.696 51.80 74.03    1      479     1378
    ##  Intercept_1    OK 20.0 21.003 16.00 25.74    1      224      424
    ##          x_1    OK  1.0  0.984  0.92  1.05    1      231      439
    ##          x_2    OK  1.0  0.984  0.92  1.05    1      231      439
    ##      sigma_1    OK  5.0  4.888  4.42  5.42    1     4577     4375
    ##        ar1_1    OK  0.8  0.779  0.67  0.89    1     2343     5798
    ##        ar1_2    OK  0.2  0.094 -0.19  0.39    1     2248     2863
    ## 
    ## Warning: 3 parameters show poor convergence (Rhat > 1.01 or ESS < 400).

We can also plot some of the parameters. As usual, we see that the
change point is not well defined by any known distribution. The fact
that the posterior mean is around 60 does not (necessarily) mean that
there is a high credence in this value. Usually, I find that any bi- or
N-modality on the posterior matches well with what you would guess from
looking at the raw data. As they say: Bayesian inference is common sense
applied to data.

``` r

plot_pars(fit, regex_pars = "cp_1|ar_*")
```

![](arma_files/figure-html/unnamed-chunk-15-1.png)

As usual, we can test hypotheses ([read more
here](https://lindeloev.github.io/mcp/dev/articles/comparison.md)). We
can also ask how much more likely it is (relative to the prior) that
there the two autocorrelations are equal compared to them differing.
Because we sampled both the prior and posterior
(`mcp(..., sample = "both")`), we can do a Savage-Dickey density ratio
test:

``` r

hypothesis(fit, "ar1_1 = ar1_2")
```

    ##          hypothesis      mean     lower     upper           p          BF
    ## 1 ar1_1 - ar1_2 = 0 0.6846472 0.3820554 0.9841136 0.002109506 0.002113965

In this case, the evidence for equality is so small that it is rounded
to zero. This means that not even a single MCMC sample visited a state
with noticeable density at zero difference.

Of course, we can also do directional tests. For example, what is the
evidence that `ar1_1` is more than 0.3 greater than `ar1_2`? Answer:
around 100 to one.

``` r

hypothesis(fit, "ar1_1 - 0.3 > ar1_2")
```

    ##                hypothesis      mean      lower     upper     p       BF
    ## 1 ar1_1 - 0.3 - ar1_2 > 0 0.3846472 0.08205536 0.6841136 0.993 301.4464

## Priors on AR and MA coefficients

The default prior on each AR and MA intercept is a truncated normal
distribution:

``` r
dnorm(0, 0.5) T(-1, 1)
```

It is symmetric around zero and gently shrinks away from the boundaries:
its central 95% interval is approximately $`[-0.84, 0.84]`$. This is a
regularizing prior on each direct coefficient, not a joint stationarity
or invertibility constraint for orders above one.

`brms` also models AR and MA coefficients directly within $`[-1, 1]`$,
but uses flat default priors over that range. Thus `mcp` keeps the
familiar parameter interpretation while adding mild regularization
toward zero; fits should generally be similar when the likelihood is
informative.

Categorical contrasts default to `dnorm(0, 0.25)`. Numeric slopes are
also normal and scaled so that one representative predictor change has
prior SD 0.25. Thus approximately 95% of the prior change lies within
$`\pm 0.5`$.

The defaults express no assumption about the sign. For a time series
where positive first-order autocorrelation is expected, alternatives
include `ar1_1 = "dunif(0, 1)"` or `ar1_1 = "dnorm(0.5, 0.5) T(0, 1)"`.
[Read more about specifying and checking
priors](https://lindeloev.github.io/mcp/dev/articles/priors.md).

Here is a complete list of the (default) priors in the model above:

``` r

cbind(fit$prior)
```

    ##             [,1]                    
    ## cp_1        "dunif(0, 100)"         
    ## ar1_1       "dnorm(0, 0.5) T(-1, 1)"
    ## ar1_2       "dnorm(0, 0.5) T(-1, 1)"
    ## Intercept_1 "dt(72.4, 35.2, 3)"     
    ## x_1         "dt(0, 0.704, 3)"       
    ## x_2         "x_1"                   
    ## sigma_1     "dt(0, 35.2, 3) T(0, )"

We can also visualize them because we sampled the prior. `prior = TRUE`
works in most `mcp` functions, including
[`plot()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
and
[`summary()`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md).

``` r

plot_pars(fit, prior = TRUE)
```

![](arma_files/figure-html/unnamed-chunk-19-1.png)![](arma_files/figure-html/unnamed-chunk-19-2.png)

Notice that the posteriors are smoothed at sharp cutoffs, slightly
misrepresenting the true distribution.

Let’s inspect the priors for a more advanced AR model, since you would
often have to inform these:

``` r

model = list(
  y ~ 1 + ar(2, 1 + x),
  ~ 0 + ar(1, 1 + I(x^2))
)
empty = mcp(model, data = data.frame(y = 1:10, x = 1:10), sample = FALSE)
cbind(empty$prior)
```

    ##             [,1]                    
    ## cp_1        "dunif(1, 10)"          
    ## ar1_1       "dnorm(0, 0.5) T(-1, 1)"
    ## ar1_x_1     "dnorm(0, 0.05555556)"  
    ## ar2_1       "dnorm(0, 0.5) T(-1, 1)"
    ## ar2_x_1     "dnorm(0, 0.05555556)"  
    ## ar1_2       "dnorm(0, 0.5) T(-1, 1)"
    ## ar1_xE2_2   "dnorm(0, 0.01234568)"  
    ## Intercept_1 "dt(5.5, 3.7, 3)"       
    ## sigma_1     "dt(0, 3.7, 3) T(0, )"

AR and MA coefficient regressions use an identity link, while their
residual inputs use the supported response-family link described above.
Careful priors become especially important with higher orders or
coefficient slopes, where individual values in \[-1, 1\] do not ensure a
stable recurrence. After fitting,
[`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) checks
up to 500 posterior draws at up to 100 observed predictor rows and warns
if more than 10% of the checked draws violate the root condition at any
checked row. `fit$simulate()` checks the supplied coefficient trajectory
before generating a fresh series. For predictor- or segment-varying
coefficients these are local smoke tests, not proofs of global
stationarity or invertibility.

Here are a few ways in which you may want to inform the AR or MA
parameters:

- For a constant AR(2) process, stationarity requires
  $`-1 < \phi_2 < 1`$ and $`\phi_2 - 1 < \phi_1 < 1 - \phi_2`$. One
  dependent-prior construction that stays in this region is
  `ar2_1 = "dunif(-1, 1)"` together with
  `ar1_1 = "dunif(ar2_1 - 1, 1 - ar2_1)"`. The simpler suggestion
  `ar2_1 = "dunif(0, ar1_1)"` is not sufficient. Higher orders require
  the corresponding joint root condition.
- Slopes can quickly make AR coefficients non-stationary or MA
  coefficients non-invertible. Constrain their magnitude with reference
  to a plausible change over the predictor range. For example, a shallow
  negative change over the observed x-span could use
  `"dnorm(0, 0.1 / (max(x) - min(x))) T(, 0)"`.

## JAGS code

Here is the JAGS code for the second simulation example, i.e., the one
with a single slope going from AR(2) to AR(1). You can print
`fit$simulate` and see that it runs much of the same code.

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
    ##   ar1_1 ~ dnorm(0, 1/(0.5)^2) T(-1,1)  # Zero-centered regularizing dependence coefficient
    ##   ar1_2 ~ dnorm(0, 1/(0.5)^2) T(-1,1)  # Zero-centered regularizing dependence coefficient
    ##   Intercept_1 ~ dt(72.4, 1/(35.2)^2, 3)   # Robustly centered mean intercept with a minimum scale of 2.5
    ##   x_1 ~ dt(0, 1/(0.704)^2, 3)   # Regularizing mean coefficient scaled to a reference predictor change
    ##   x_2 = x_1  # Same value as x_1
    ##   sigma_1 ~ dt(0, 1/(35.2)^2, 3) T(0,)  # Positive residual SD calibrated on the response scale
    ## 
    ##   # Apply GARMA recursion to link-scale residuals
    ##   resid_arma_[1] = 0
    ##   for (i_ in 2:length(x)) {
    ##     resid_arma_[i_] = ar1_[i_] * resid_abs_[i_ - 1]
    ##   }
    ##   # Model and likelihood
    ##   for (i_ in 1:length(x)) {
    ##     # par_x local to each segment
    ##     x_local_1_[i_] = min(x[i_], cp_1)
    ##     x_local_2_[i_] = min(x[i_], cp_2) - cp_1
    ##     
    ##     # GARMA observation boundary
    ##     garma_boundary_[i_] =
    ##       (x[i_] < cp_1) * 0.1 +
    ##       (x[i_] >= cp_1) * 0.1
    ##     
    ##     # Formula for ar1
    ##     ar1_[i_] =
    ##     
    ##       # Segment 1: y ~ 1 + x + ar(1)
    ##       (x[i_] >= cp_0) * (x[i_] < cp_1) * inprod(rhs_matrix_[i_, c(1)], c(ar1_1)) * 1 + 
    ##     
    ##       # Segment 2: y ~ 1 ~ 0 + x + ar(1)
    ##       (x[i_] >= cp_1) * inprod(rhs_matrix_[i_, c(2)], c(ar1_2)) * 1
    ##     
    ##     # Formula for mu
    ##     link_mu_[i_] =
    ##     
    ##       # Segment 1: y ~ 1 + x + ar(1)
    ##       (x[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(3)], c(Intercept_1)) * 1 + 
    ##       (x[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(4)], c(x_1)) * x_local_1_[i_] + 
    ##     
    ##       # Segment 2: y ~ 1 ~ 0 + x + ar(1)
    ##       (x[i_] >= cp_1) * inprod(rhs_matrix_[i_, c(5)], c(x_2)) * x_local_2_[i_]
    ##     
    ##     # Formula for sigma
    ##     link_sigma_[i_] =
    ##     
    ##       # Segment 1: y ~ 1 + x + ar(1)
    ##       (x[i_] >= cp_0) * inprod(rhs_matrix_[i_, c(6)], c(sigma_1)) * 1
    ## 
    ##     # Likelihood and log-density for family = gaussian()
    ##     mu_[i_] = link_mu_[i_] + resid_arma_[i_]
    ##     sigma_[i_] = max(1e-03, link_sigma_[i_])
    ##     y[i_] ~ dnorm(mu_[i_], 1 / sigma_[i_]^2)  # SD as precision
    ##     garma_y_[i_] = y[i_]
    ##     garma_link_y_[i_] = garma_y_[i_]
    ##     resid_abs_[i_] = garma_link_y_[i_] - link_mu_[i_]
    ##   }
    ## }
