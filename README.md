# mcp: Bayesian Inference of Multiple Change Points

[![Travis-CI status](https://img.shields.io/travis/lindeloev/mcp.svg)](https://travis-ci.org/lindeloev/mcp)
[![Coveralls status](https://codecov.io/gh/lindeloev/mcp/branch/master/graph/badge.svg)](https://coveralls.io/r/lindeloev/mcp)


Bayesian inference of Hierarchical Multiple Change Points (MCP) between Generalized Linear Segments - and the coefficients in those segments. 

Change points are also sometimes called **switch points**, **break points**, **broken-line** models, **bilinear** models, (performance) **discontinuity** models, **piecewise linear regression**, and **local linear regression**. `mcp` aims to be be useful for all of them. See how `mcp` compares to [other packages](https://lindeloev.github.io/mcp/articles/packages.html).

Under the hood, `mcp` takes a formula-representation of linear segments and turns it into [JAGS](https://sourceforge.net/projects/mcmc-jags/) code. `mcp` leverages the power of `tidybayes`, `bayesplot`, `coda`, and `loo` to make change point analysis easy and powerful.

# Install

 1. [Install the latest version of JAGS](https://sourceforge.net/projects/mcmc-jags/). Linux users can fetch binaries [here](http://mcmc-jags.sourceforge.net/).
 
 2. Install `mcp` by running this in R:
 
    ```r
    if (!requireNamespace("remotes")) install.packages("remotes")
    remotes::install_github("lindeloev/mcp")
    ```



# Brief example
`mcp` takes a list of formulas for `y` as a function of `x`. The change point(s) are the `x` at which data changes from being better predicted by one formula to the next. The first formula is just `response ~ predictors` and the following formulas typically take the form `changepoint ~ predictors` (more details [here](https://lindeloev.github.io/mcp/articles/formulas.html)).

The following model infers the two change points between three segments.

```r
library(mcp)

# Define the model
segments = list(
  response ~ 1,  # plateau (int_1)
  1 ~ 0 + time,  # joined slope (time_2) at cp_1
  1 ~ 1 + time  # disjoined slope (int_1, time_2) at cp_2
)

# Fit it. The `ex_demo` dataset is included in mcp
fit = mcp(segments, data = ex_demo)
```

Plot lines drawn randomly from the posterior on top of data to inspect the fit:
```r
plot(fit)
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_demo.png)

Use `summary` to summarise the posterior distribution as well as sampling diagnostics. They were simulated to lie at `cp_1 = 30` and `cp_2 = 70` and these values are well recovered, as are the other coefficients ([see how `ex_demo` was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_demo.R)). 

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
   name   mean  lower   upper rhat  eff   ts_se
   cp_1 30.510 22.832 38.6869 1.01  366 378.575
   cp_2 69.776 69.273 70.2350 1.00 5868   0.131
  int_1 10.279  8.839 11.6956 1.00  952   4.274
  int_3  0.576 -2.454  3.6310 1.00  711  31.702
  sigma  4.008  3.464  4.6100 1.00 4495   0.168
 time_2  0.536  0.403  0.6702 1.01  421   0.106
 time_3 -0.223 -0.403 -0.0494 1.00  725   0.105
```

`rhat` is the [Gelman-Rubin convergence diagnostic](https://www.rdocumentation.org/packages/coda/versions/0.19-3/topics/gelman.diag), `eff` is the [effective sample size](https://mc-stan.org/docs/2_18/reference-manual/effective-sample-size-section.html), and `ts_se` is the time-series standard error.

`plot(fit, "combo")` can be used to inspect the posteriors and convergence of all parameters. See the documentation of `plot.mcpfit` for many other plotting options. Here we plot just the (population-level) change points. They often have "strange" posterior distributions, highlighting the need for a computational approach:

```r
plot(fit, "combo", regex_pars = "cp_")
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_demo_combo.png)


For model comparisons, we can fit a null model and compare the predictive performance of the two models using (approximate) leave-one-out cross-validation. Our null model omits the first plateau and change point, essentially testing the credence of that change point:

```r
# Define the model
segments_null = list(
  response ~ 1 + time,  # intercept (int_1) and slope (time_1)
  1 ~ 1 + time  # disjoined slope (int_2, time_1)
)

# Fit it
fit_null = mcp(segments_null, ex_demo)
```

Leveraging the power of `loo::loo`, we see that the two-change-points model is preferred (it is on top), but not very strongly (`elpd_diff / se_diff` ratio). See the documentation of the `loo` package and associated papers for more details, or scroll down for more notes on model comparison using `mcp`.

```r
fit$loo = loo(fit)
fit_null$loo = loo(fit_null)

loo::loo_compare(fit$loo, fit_null$loo)
```
```
       elpd_diff se_diff
model1  0.0       0.0   
model2 -7.8       4.7   
```



# Highlights from in-depth guides
The articles in the web site's menu go in-depth with the functionality of `mcp`. Here is an executive summary, to give you a sense.

[About mcp formulas and models](https://lindeloev.github.io/mcp/articles/formulas.html):
 * Parameter names `int_i` (intercepts), `cp_i` (change points), `x_i` (slopes).
 * The change point model is basically an `ifelse` model.
 * Use `rel()` to specify that parameters are relative to those in the previous  segments.
 * Generate data using `fit$func_y`.

[Using priors in mcp](https://lindeloev.github.io/mcp/articles/priors.html) include:
 * See priors in `fit$prior`.
 * Set priors using `mcp(segments, data, prior = list(cp_1 = "dnorm(0, 1)", cp_1 = "dunif(0, 45))`.
 * Fix parameters to specific values using `cp_1 = 45`.
 * Share parameters between segments using `slope_1 = "slope_2"`.
 * Allows for truncated priors using `T(lower, upper)`, e.g., `int_1 = "dnorm(0, 1) T(0, )"`. `mcp` applies this automatically to change point priors to enforce order restriction. This is true for [varying change points](https://lindeloev.github.io/mcp/articles/varying.html) too.
 * Do prior predictive checks using `mcp(segments, data, sample="prior")`.

[Varying change points in mcp](https://lindeloev.github.io/mcp/articles/varying.html) include:
 * Simulate varying change points using `fit$func_y()`.
 * Get posteriors using `ranef(fit)`.
 * Plot using `plot(fit, facet_by="my_group")` and `plot(fit, "dens_overlay", pars = "varying", ncol = 3)`.
 * The default priors restrict varying change points to lie between the two adjacent change points.

`mcp` currently supports the following GLM:
 * `gaussian(link = "identity")` (default). See examples above and below.
 * `binomial(link = "logit")`. See [binomial change points in mcp](https://lindeloev.github.io/mcp/articles/binomial.html).
 * `bernoulli(link = "logit")`. See [binomial change points in mcp](https://lindeloev.github.io/mcp/articles/binomial.html).
 * `poisson(link = "log")`. See [Poisson change points in mcp](https://lindeloev.github.io/mcp/articles/poisson.html).

[Model comparison and hypothesis testing in mcp](https://lindeloev.github.io/mcp/articles/comparison.html):
 * Do Leave-One-Out Cross-Validation using `loo(fit)` and `loo::loo_compare(fit1$loo, fit2$loo)`.
 * Compute Savage-Dickey density rations using `hypothesis(fit, "cp_1 = 40")`.
 * Leverage directional and conditional tests to assess interval hypotheses (`hypothesis(fit, "cp_1 > 30 & cp_1 < 50")`), combined hypotheses (`hypothesis(fit, "cp_1 > 30 & int_1 > int_2")`), etc.



# Some examples
`mcp` aims to support a wide variety of models. Here are some example models for inspiration.


## Two plateaus
Find the single change point between two plateaus ([see how this data was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_plateaus.R))

```r
segments = list(
    y ~ 1,  # plateau (int_1)
    1 ~ 1  # plateau (int_2)
)
fit = mcp(segments, ex_plateaus, par_x = "x")
plot(fit)
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_plateaus.png)


## Varying change points

Here, we find the single change point between two joined slopes. While the slopes are shared by all participants, the change point varies by `id`.

Read more about [varying change points in mcp](https://lindeloev.github.io/mcp/articles/varying.html) and [see how this data was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_varying.R).

```r
segments = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope, varying by id
)
fit = mcp(segments, ex_varying)
plot(fit, facet_by = "id")
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_varying.png)

Summarise the individual change points using `ranef()` or plot them using `plot(fit, "combo", "varying")`:

```r
ranef(fit)
```

```r
            name       mean      lower      upper     rhat  eff     ts_se
1 cp_1_id[Benny] -18.939028 -21.662758 -16.178849 1.000739 2539  4.801090
2  cp_1_id[Bill] -10.932586 -13.570255  -8.355125 1.001008  709 22.656037
3  cp_1_id[Cath]  -3.868473  -6.161151  -1.534219 1.000035  996 24.755440
4  cp_1_id[Erin]   6.298535   3.689802   9.089925 1.000215 5308  2.047966
5  cp_1_id[John]   9.742310   7.184306  12.258944 1.000451 1854  7.511643
6  cp_1_id[Rose]  17.699242  14.393830  21.048751 1.000152 4019  4.435384
```


## Generalized linear models
`mcp` supports Generalized Linear Modeling. See extended examples using [`binomial()`](https://lindeloev.github.io/mcp/articles/binomial.html) and [`poisson()`](https://lindeloev.github.io/mcp/articles/poisson.html). [See how this data was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_binomial.R).

Here is a binomial change point model with three segments. We plot the 95% HDI too:

```r
segments = list(
  y | trials(N) ~ 1,  # constant rate
  1 ~ 0 + x,  # joined changing rate
  1 ~ 1 + x  # disjoined changing rate
)
fit = mcp(segments, ex_binomial, family = binomial())
plot(fit, quantiles = TRUE)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_binomial.png)

Use `plot(fit, rate = FALSE)` if you want the points and fit lines on the original scale of `y` rather than divided by `N`.


## Using `rel()` and priors
Read more about [formula options](https://lindeloev.github.io/mcp/articles/formulas.html) and [priors](https://lindeloev.github.io/mcp/articles/priors.html) and [see how this data was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_rel_prior.R).

Here we find the two change points between three segments. The slope and intercept of segment 2 are parameterized relative to segment 1, i.e., modeling the *change* in intercept and slope since segment 1. So too with the second change point (`cp_2`) which is now the *distance* from `cp_1`. 

Some of the default priors are overwritten. The first intercept (`int_1`) is forced to be 10, the slopes are in segment 1 and 3 is shared. It is easy to see these effects in the `ex_rel_prior` dataset because they violate it somewhat. The first change point has to be at `x = 20` or later.

```r
segments = list(
  y ~ 1 + x,
  1 ~ rel(1) + rel(x),
  rel(1) ~ 0 + x
)

prior = list(
  int_1 = 10,  # Constant, not estimated
  x_3 = "x_1",  # shared slope in segment 1 and 3
  int_2 = "dnorm(0, 20)",
  cp_1 = "dunif(20, 50)"  # has to occur in this interval
)
fit = mcp(segments, ex_rel_prior, prior)
plot(fit)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_fix_rel.png)

Looking at the summary, we can see that `int_2` and `x_2` are relative values:
```r
summary(fit)
```

```r
Family: gaussian(link = 'identity')
Iterations: 9000 from 3 chains.
Segments:
  1: y ~ 1 + x
  2: y ~ 1 ~ rel(1) + rel(x)
  3: y ~ rel(1) ~ 0 + x

Population-level parameters:
  name  mean lower upper rhat  eff    ts_se
  cp_1 22.46 20.00 26.59 1.21   21 1.68e+03
  cp_2 46.55 41.60 49.80 1.21   29 1.22e+03
 int_1 10.00 10.00 10.00  NaN    0 0.00e+00
 int_2  6.42 -8.84 18.37 1.13   36 1.10e+04
 sigma  7.09  6.15  8.14 1.00 2943 8.49e-01
   x_1  1.85  1.66  2.04 1.03   57 1.37e+00
   x_2 -3.94 -4.17 -3.70 1.01  117 1.16e+00
   x_3  1.85  1.66  2.04 1.03   57 1.37e+00
```

(rhat cannot be estimated when the model contains constant coefficients in `mcp 0.1`. This will be fixed.)


## Quadratic
Write exponents as `I(x^N)`. E.g., `I(x^2)`, `I(x^1.5)`, `I(x^3)`. The following detects the onset of linear + quadratic growth. In nuitrition, this is often called the BLQ model (Broken Line Quadratic).

```r
segments = list(
  y ~ 1,
  1 ~ 0 + x + I(x^2)
)
fit = mcp(segments, ex_quadratic)
plot(fit)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_quadratic.png)


## Trigonometric and other 
You can use `sin(x)`, `cos(x)`, and `tan(x)` to do trigonometry. You can also do `exp(x)`, `abs(x)`, `log(x)`, and `sqrt(x)`, but beware that the two latter will fail in segment 2+ because from the "perspective" of that segment, earlier `x` values are negative.

```r
segments = list(
  y ~ 1 + sin(x),
  1 ~ 0 + cos(x) + x
)

fit = mcp(segments, ex_trig)
plot(fit)
```

![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_trig.png)



# Diagnosing problems
**Convergence:** A common problem when using MCMC is lacking convergence between chains. This will show up as large `rhat` values (> 1.1 is a common criterion) and non-converging lines in `plot(fit, "trace")`. 

 * The first thing to try is always to make the model warm up longer to see if it reaches convergence later: `mcp(fit, data, adapt = 4000, update = 4000)`. 
 
 * It can be a sign of a deeper non-identifiability in the model. This will show up as strong correlations in the joint distribution of any pair of implicated parameters: `plot(fit, "hex", pars = c("int_1", "int_2))`. This may give you ideas how to change the model.
 
 * You can set the initial values for the JAGS sampler using, e.g., `mcp(..., inits = list(cp_1 = 20, int_2 = -20, etc.))`. This will be passed to `jags.fit` and you can see more documentation there.


**Speed:** A lot of data and complicated models will slow down fitting. 

 * Run the chains in parallel using, e.g., `mcp(..., chains=4, cores=4)`.

 * More data usually means better identifiability and faster convergence. Lower the warmup period using, e.g., `mcp(..., adapt=300, update = 200)`.

**Won't run:** Most of these problems should stem from inappropriate priors and such problems may be exacerbated by fragile link functions (e.g., `binomial(link = "identity")`. The article on [priors in mcp](https://lindeloev.github.io/mcp/articles/priors.html) may be helpful, but in particular:

 * Errors on "directed cycle" usually stems from using parameters in priors. For example, if you set `prior = list(int_1 = "dnorm(int_2, 1)"", int_2 = "dnorm(int_1, 1)")` this is an infinite regress.
 
 * Errors on "incompatible with parent nodes" usually stem from impossible values. For example, if you set `prior = list(sigma = "dnorm(0, 1)"")`, this allows for a negative standard deviation, which is impossible. Think about changing the prior distributions and perhaps truncate them using `T(lower, upper)`.


If you encounter these or other problems, don't hesitate to [raise a Github Issue](https://github.com/lindeloev/mcp/issues), asking for help or posting a bug report.




# Do much more
Don't be constrained by these simple `mcp` functions. `fit$samples` is a regular `mcmc.list` object and all methods apply. You can work with the MCMC samples just as you would with `brms`, `stan_glm`, `jags`, or other samplers using the always excellent `tidybayes`:

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
