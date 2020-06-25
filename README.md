# mcp: Regression with Multiple Change Points

[![mcp Travis-CI status](https://travis-ci.org/lindeloev/mcp.svg?branch=master)](https://travis-ci.org/lindeloev/mcp)
[![mcp Coveralls status](https://codecov.io/gh/lindeloev/mcp/branch/master/graph/badge.svg)](https://coveralls.io/r/lindeloev/mcp)
[![mcp CRAN status](https://www.r-pkg.org/badges/version/mcp)](https://CRAN.R-project.org/package=mcp)
[![mcp CRAN downloads](https://cranlogs.r-pkg.org/badges/mcp)](https://cranlogs.r-pkg.org/badges/mcp)

`mcp` does regression with one or Multiple Change Points (MCP) between Generalized and hierarchical Linear Segments using Bayesian inference. `mcp` aims to provide maximum flexibility for analyses with a priori knowledge about the number of change points and the form of the segments in between.

Change points are also called **switch points**, **break points**, **broken line** regression, **broken stick** regression, **bilinear** regression, **piecewise linear** regression, **local linear** regression, **segmented** regression, and (performance) **discontinuity** models. `mcp` aims to be be useful for all of them. See how `mcp` compares to [other R packages](https://lindeloev.github.io/mcp/articles/packages.html).

Under the hood, `mcp` takes a formula-representation of linear segments and turns it into [JAGS](https://sourceforge.net/projects/mcmc-jags/) code. `mcp` leverages the power of `tidybayes`, `bayesplot`, `coda`, and `loo` to make change point analysis easy and powerful.


# Install

 1. <a href="https://sourceforge.net/projects/mcmc-jags/" target="_blank">Install the latest version of JAGS</a>. Linux users can fetch binaries <a href="http://mcmc-jags.sourceforge.net/" target="_blank">here</a>.
 
 2. Install from CRAN:
 
    ```r
    install.packages("mcp")
    ```
    
    or install the development version from GitHub:
 
    ```r
    if (!requireNamespace("remotes")) install.packages("remotes")
    remotes::install_github("lindeloev/mcp")
    ```


# At a glance
Here are some example `mcp` models. `mcp` takes a list of formulas - one for each segment. The change point(s) are the `x` at which data changes from being better predicted by one formula to the next. The first formula is just `response ~ predictors` and the most common formula for segment 2+ would be ` ~ predictors` (more details [here](https://lindeloev.github.io/mcp/articles/formulas.html)).


![](https://github.com/lindeloev/mcp-paper/raw/master/all_plots_small.png)

Scroll down to see brief introductions to each of these, or browse the website articles for more thorough worked examples and discussions.



# Brief worked example
## Fit a model
The following model infers the two change points between three segments.

```r
library(mcp)

# Define the model
model = list(
  response ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1 + time     # disjoined slope (int_3, time_3) at cp_2
)

# Fit it. The `ex_demo` dataset is included in mcp
fit = mcp(model, data = ex_demo)
```

## Plot and summary

The default plot includes data, fitted lines drawn randomly from the posterior, and change point(s) posterior density for each chain:
```r
plot(fit)
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_demo.png)

Use `summary()` to summarise the posterior distribution as well as sampling diagnostics. They were [simulated with mcp](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_demo.R) so the summary include the "true" values in the column `sim` and the column `match` show whether this true value is within the interval:

```r
summary(fit)
```
```r
Family: gaussian(link = 'identity')
Iterations: 9000 from 3 chains.
Segments:
  1: response ~ 1
  2: response ~ 1 ~ 0 + time
  3: response ~ 1 ~ 1 + time

Population-level parameters:
    name match  sim  mean lower  upper Rhat n.eff
    cp_1    OK 30.0 30.27 23.19 38.760    1   384
    cp_2    OK 70.0 69.78 69.27 70.238    1  5792
   int_1    OK 10.0 10.26  8.82 11.768    1  1480
   int_3    OK  0.0  0.44 -2.49  3.428    1   810
 sigma_1    OK  4.0  4.01  3.43  4.591    1  3852
  time_2    OK  0.5  0.53  0.40  0.662    1   437
  time_3    OK -0.2 -0.22 -0.38 -0.035    1   834
```

`rhat` is the [Gelman-Rubin convergence diagnostic](https://www.rdocumentation.org/packages/coda/versions/0.19-3/topics/gelman.diag), `eff` is the [effective sample size](https://mc-stan.org/docs/2_18/reference-manual/effective-sample-size-section.html).

`plot_pars(fit)` can be used to inspect the posteriors and convergence of all parameters. See the documentation of `plot_pars()` for many other plotting options. Here, we plot just the (population-level) change points. They often have "strange" posterior distributions, highlighting the need for a computational approach:

```r
plot_pars(fit, regex_pars = "cp_")
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_demo_combo.png)


## Tests and model comparison

We can test (joint) probabilities in the model using `hypothesis()` ([see more here](https://lindeloev.github.io/mcp/articles/comparison.html)). For example, what is the evidence (given priors) that the first change point is later than 25 against it being less than 25?

```r
hypothesis(fit, "cp_1 > 25")
```
```r
     hypothesis mean lower upper     p   BF
1 cp_1 - 25 > 0 5.27 -1.81 13.76 0.917 11.1
```


For model comparisons, we can fit a null model and compare the predictive performance of the two models using (approximate) leave-one-out cross-validation ([see more here](https://lindeloev.github.io/mcp/articles/comparison.html)). Our null model omits the first plateau and change point, essentially testing the credence of that change point:

```r
# Define the model
model_null = list(
  response ~ 1 + time,  # intercept (int_1) and slope (time_1)
  ~ 1 + time            # disjoined slope (int_2, time_1)
)

# Fit it
fit_null = mcp(model_null, ex_demo)
```

Leveraging the power of `loo::loo`, we see that the two-change-points model is preferred (it is on top), but the `elpd_diff / se_diff` ratio ratio indicate that this preference is not very strong.
```r
fit$loo = loo(fit)
fit_null$loo = loo(fit_null)

loo::loo_compare(fit$loo, fit_null$loo)
```
```
       elpd_diff se_diff
model1  0.0       0.0
model2 -7.6       4.6
```



# Highlights from in-depth guides
The articles on the [mcp website](https://lindeloev.github.io/mcp) go in-depth with the functionality of `mcp`. Here is an executive summary, to give you a quick sense of what mcp can do.

[About mcp formulas and models](https://lindeloev.github.io/mcp/articles/formulas.html):
 * Parameter names are `int_i` (intercepts), `cp_i` (change points), `x_i` (slopes), `phi_i` (autocorrelation), and `sigma_*` (variance).
 * The change point model is basically an `ifelse` model.
 * Use `rel()` to specify that parameters are relative to those corresponding in the previous segments.
 * Generate data using `fit$simulate()`.

[Using priors](https://lindeloev.github.io/mcp/articles/priors.html):
 * See priors in `fit$prior`.
 * Set priors using `mcp(..., prior = list(cp_1 = "dnorm(0, 1)", cp_1 = "dunif(0, 45)")`.
 * The default prior for change points is fast for estimation but is mathematically "messy". The Dirichlet prior (`cp_i = "dirichlet(1)"`) is slow but beautiful.
 * Fix parameters to specific values using `cp_1 = 45`.
 * Share parameters between segments using `slope_1 = "slope_2"`.
 * Truncate priors using `T(lower, upper)`, e.g., `int_1 = "dnorm(0, 1) T(0, )"`. `mcp` applies this automatically to change point priors to enforce order restriction. This is true for [varying change points](https://lindeloev.github.io/mcp/articles/varying.html) too.
 * Do prior predictive checks using `mcp(model, data, sample = "prior")`.

[Varying change points](https://lindeloev.github.io/mcp/articles/varying.html):
 * Simulate varying change points using `fit$simulate()`.
 * Get posteriors using `ranef(fit)`.
 * Plot using `plot(fit, facet_by="my_group")` and `plot_pars(fit, pars = "varying", type = "dens_overlay", ncol = 3)`.
 * The default priors restrict varying change points to lie between the two adjacent change points.

[Supported families and link functions](https://lindeloev.github.io/mcp/articles/families.html): 
 * `mcp` currently supports specific combinations of families (`gaussian()`, `binomial()`, `bernoulli()`, and `poisson()`) and link functions (`identity`, `logit`, `probit`, and `log`).
 * Use informative priors to avoid issues when using non-default priors.
 * Use `binomial(link = "logit")` for [binomial change points in mcp](https://lindeloev.github.io/mcp/articles/binomial.html). Also relevant for `bernoulli(link = "logit")`.
 * Use `poisson(link = "log")` for [Poisson change points in mcp](https://lindeloev.github.io/mcp/articles/poisson.html).

[Model comparison and hypothesis testing](https://lindeloev.github.io/mcp/articles/comparison.html):
 * Do Leave-One-Out Cross-Validation using `loo(fit)` and `loo::loo_compare(fit1$loo, fit2$loo)`.
 * Compute Savage-Dickey density rations using `hypothesis(fit, "cp_1 = 40")`.
 * Leverage directional and conditional tests to assess interval hypotheses (`hypothesis(fit, "cp_1 > 30 & cp_1 < 50")`), combined other hypotheses (`hypothesis(fit, "cp_1 > 30 & int_1 > int_2")`), etc.

Modeling [variance](https://lindeloev.github.io/mcp/articles/variance.html) and [autoregression](https://lindeloev.github.io/mcp/articles/arma.html):
 * `~ sigma(1)` models an intercept change in variance. `~ sigma(0 + x)` models increasing/decreasing variance.
 * `~ ar(N)` models Nth order autoregression on residuals. `~ar(N, 0 + x)` models increasing/decreasing autocorrelation.
 * You can model anything for `sigma()` and `ar()`. For example, `~ x + sigma(1 + x + I(x^2))` models polynomial change in variance with `x` on top of a slope on the mean.
 * Simulate effects and change points on `sigma()` and `ar()` using `fit$simulate()`

[Tips, tricks, and debugging](https://lindeloev.github.io/mcp/articles/debug.html)
 * Speed up fitting using `mcp(..., cores = 3)` / `options(mcp_cores = 3)`, and/or fewer iterations, `mcp(..., adapt = 500)`.
 * Help convergence along using `mcp(..., inits = list(cp_1 = 20, int_2 = -3))`.
 * Most errors will be caused by circularly defined priors.



# Some examples
`mcp` aims to support a wide variety of models. Here are some example models for inspiration.


## Means
Find the single change point between two plateaus ([see how this data was simulated with mcp](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_plateaus.R)).

```r
model = list(
    y ~ 1,  # plateau (int_1)
    ~ 1     # plateau (int_2)
)
fit = mcp(model, ex_plateaus, par_x = "x")
plot(fit)
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_plateaus.png)


## Varying change points

Here, we find the single change point between two joined slopes. While the slopes are shared by all participants, the change point varies by `id`. Read more about [varying change points in mcp](https://lindeloev.github.io/mcp/articles/varying.html).

```r
model = list(
  y ~ 1 + x,          # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope, varying by id
)
fit = mcp(model, ex_varying)
plot(fit, facet_by = "id")
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_varying.png)

Summarise the varying change points using `ranef()` or plot them using `plot_pars(fit, "varying")`. Again, [this data was simulated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_varying.R) so the columns `match` and `sim` are added to show simulation values and whether they are inside the interval. Set the `width` wider for a more lenient criterion.

```r
ranef(fit, width = 0.98)
```

```r
           name match   sim  mean   lower   upper Rhat n.eff
 cp_1_id[Benny]    OK -17.5 -18.1 -21.970 -14.877    1   895
  cp_1_id[Bill]    OK -10.5  -7.6 -10.658  -4.451    1   420
  cp_1_id[Cath]    OK  -3.5  -2.8  -5.634   0.027    1   888
  cp_1_id[Erin]    OK   3.5   3.1   0.041   5.952    1  3622
  cp_1_id[John]    OK  10.5  11.3   7.577  14.989    1  2321
  cp_1_id[Rose]    OK  17.5  14.1  10.485  18.079    1  5150
```


## Generalized linear models
`mcp` supports Generalized Linear Modeling. See extended examples using [`binomial()`](https://lindeloev.github.io/mcp/articles/binomial.html) and [`poisson()`](https://lindeloev.github.io/mcp/articles/poisson.html). These data were simulated with `mcp` [here](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_binomial.R).

Here is a binomial change point model with three segments. We plot the 95% HDI too:

```r
model = list(
  y | trials(N) ~ 1,  # constant rate
  ~ 0 + x,            # joined changing rate
  ~ 1 + x             # disjoined changing rate
)
fit = mcp(model, ex_binomial, family = binomial())
plot(fit, q_fit = TRUE)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_binomial.png)

Use `plot(fit, rate = FALSE)` if you want the points and fit lines on the original scale of `y` rather than divided by `N`.



## Time series
`mcp` allows for flexible time series analysis with autoregressive residuals of arbitrary order. Below, we model a change from a plateau with strong positive AR(2) residuals to a slope with medium AR(1) residuals. These data were simulated with mcp [here](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_ar.R) and the generating values are in the `sim` column. You can also do regression on the AR coefficients themselves using e.g., `ar(1, 1 + x)`. [Read more here](https://lindeloev.github.io/mcp/articles/arma.html).

```r
model = list(
  price ~ 1 + ar(2),
  ~ 0 + time + ar(1)
)
fit = mcp(model, ex_ar)
summary(fit)
```

The AR(N) parameters on intercepts are named `ar[order]_[segment]`. All parameters, including the change point, are well recovered:

```r
Population-level parameters:
    name match   sim    mean     lower   upper Rhat n.eff
   ar1_1    OK   0.7   0.741  5.86e-01   0.892 1.01   713
   ar1_2    OK  -0.4  -0.478 -6.88e-01  -0.255 1.00  2151
   ar2_1    OK   0.2   0.145 -6.56e-04   0.284 1.01   798
    cp_1       120.0 117.313  1.14e+02 118.963 1.05   241
   int_1        20.0  17.558  1.51e+01  19.831 1.02   293
 sigma_1    OK   5.0   4.829  4.39e+00   5.334 1.00  3750
  time_2    OK   0.5   0.517  4.85e-01   0.553 1.00   661
```

The fit plot shows the inferred autocorrelated nature:

```r
plot(fit_ar)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_ar.png)



## Variance change and prediction intervals
You can model variance by adding a `sigma()` term to the formula. The inside `sigma()` can take everything that the formulas outside do. Read more in [the article on variance](https://lindeloev.github.io/mcp/articles/variance.html). The example below models two change points. The first is variance-only: variance abruptly increases and then declines linearly with `x`. The second change point is the stop of the variance-decline and the onset of a slope on the mean.

Effects on variance is best visualized using *prediction intervals*. See more in the documentation for `plot.mcpfit()`.

```r
model = list(
  y ~ 1,
  ~ 0 + sigma(1 + x),
  ~ 0 + x
)
fit = mcp(model, ex_variance, cores = 3, adapt = 5000, iter = 5000)
plot(fit, q_predict = TRUE)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_variance.png)



## Quadratic and other exponentiations
Write exponents as `I(x^N)`. E.g., quadratic `I(x^2)`, cubic `I(x^3)`, or some other power function `I(x^1.5)`. The example below detects the onset of linear + quadratic growth. This is often called the BLQ model (Broken Line Quadratic) in nutrition research.

```r
model = list(
  y ~ 1,
  ~ 0 + x + I(x^2)
)
fit = mcp(model, ex_quadratic)
plot(fit)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_quadratic.png)




## Trigonometric and others
You can use `sin(x)`, `cos(x)`, and `tan(x)` to do trigonometry. This can be useful for seasonal trends and other periodic data. You can also do `exp(x)`, `abs(x)`, `log(x)`, and `sqrt(x)`, but beware that the two latter will currently fail in segment 2+. Raise an issue if you need this.

```r
model = list(
  y ~ 1 + sin(x),
  ~ 1 + cos(x) + x
)

fit = mcp(model, ex_trig)
plot(fit)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_trig.png)





## Using `rel()` and priors
Read more about [formula options](https://lindeloev.github.io/mcp/articles/formulas.html) and [priors](https://lindeloev.github.io/mcp/articles/priors.html).

Here we find the two change points between three segments. The slope and intercept of segment 2 are parameterized relative to segment 1, i.e., modeling the *change* in intercept and slope since segment 1. So too with the second change point (`cp_2`) which is now the *distance* from `cp_1`. 

Some of the default priors are overwritten. The first intercept (`int_1`) is forced to be 10, the slopes are in segment 1 and 3 is shared. It is easy to see these effects in the `ex_rel_prior` dataset because they violate it somewhat. The first change point has to be at `x = 20` or later.

```r
model = list(
  y ~ 1 + x,
  ~ rel(1) + rel(x),
  rel(1) ~ 0 + x
)

prior = list(
  int_1 = 10,  # Constant, not estimated
  x_3 = "x_1",  # shared slope in segment 1 and 3
  int_2 = "dnorm(0, 20)",
  cp_1 = "dunif(20, 50)"  # has to occur in this interval
)
fit = mcp(model, ex_rel_prior, prior, iter = 10000)
plot(fit, cp_dens = FALSE)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_fix_rel.png)

Comparing the summary to the fitted lines in the plot, we can see that `int_2` and `x_2` are relative values. We also see that the "wrong" priors made it harder to recover the parameters used [to simulate this data](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_rel_prior.R) (`match` and `sim` columns):

```r
summary(fit)
```

```r
Population-level parameters:
    name match   sim  mean  lower upper Rhat n.eff
    cp_1    OK  25.0 23.15  20.00 25.81 1.00   297
    cp_2        40.0 51.85  47.06 56.36 1.02   428
   int_1        25.0 10.00  10.00 10.00  NaN     0
   int_2    OK -10.0 -6.86 -21.57 11.89 1.03   190
 sigma_1         7.0  9.70   8.32 11.18 1.00  7516
     x_1         1.0  1.58   1.24  1.91 1.07   120
     x_2    OK  -3.0 -3.28  -3.61 -2.96 1.04   293
     x_3         0.5  1.58   1.24  1.91 1.07   120
```




# Do much more
Don't be constrained by these simple `mcp` functions. `fit$samples` is a regular `mcmc.list` object and all methods apply. You can work with the MCMC samples just as you would with `brms`, `rstanarm`, `jags`, or other samplers using the always excellent `tidybayes`:

```r
library(tidybayes)

# Extract all parameters:
tidy_draws(fit$samples) %>%
  # tidybayes stuff here

# Extract some parameters:
fit$pars$model  # check out which parameters are inferred.
spread_draws(fit$samples, cp_1, cp_2, int_1, year_1) %>%
 # tidybayes stuff here
```


# Citation
[This preprint](https://osf.io/fzqxv) formally introduces `mcp`. Find citation info at the link, call `citation("mcp")` or copy-paste this into your reference manager:

```
  @Article{,
    title = {mcp: An R Package for Regression With Multiple Change Points},
    author = {Jonas Kristoffer Lindeløv},
    journal = {OSF Preprints},
    year = {2020},
    doi = {10.31219/osf.io/fzqxv},
    encoding = {UTF-8},
  }
```
