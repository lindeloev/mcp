# mcp: Bayesian Inference of Multiple Change Points

[![Travis-CI status](https://img.shields.io/travis/lindeloev/mcp.svg)](https://travis-ci.org/lindeloev/mcp)
[![Coveralls status](https://codecov.io/gh/lindeloev/mcp/branch/master/graph/badge.svg)](https://coveralls.io/r/lindeloev/mcp)


Bayesian inference of Hierarchical Multiple Change Points (MCP) between Generalized Linear Segments - and the parameters in those segments. Under the hood, `mcp` takes a formula-representation of linear segments and turns it into [JAGS](https://sourceforge.net/projects/mcmc-jags/) code. The rest of the package leverages the power of `tidybayes`, `bayesplot`, `coda`, and `loo` to make change point analysis easy and powerful.


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

Use `summary` to summarise the posterior distribution as well as sampling diagnostics. They were simulated to lie at `cp_1 = 30` and `cp_2 = 70` and these values are well recovered, as are the other simulation parameters ([see how `ex_demo` was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_demo.R)). 

```r
summary(fit)
```
```r
Family: gaussian(link = 'identity')
Iterations: 9000 from 3 chains.
Segments:
   response ~ 1 
   response ~ 1 ~ 0 + time 
   response ~ 1 ~ 1 + time 

Population-level parameters:
   name   mean   X2.5   X97.5 rhat  eff   ts_se
   cp_1 30.973 23.054 38.5818 1.01  396 316.789
   cp_2 69.783 69.272 70.2324 1.00 5859   0.129
  int_1 10.317  8.913 11.7922 1.01 1683   2.939
  int_3  0.579 -2.483  3.5711 1.00  712  29.767
  sigma  4.003  3.462  4.6293 1.00 4224   0.187
 time_2  0.543  0.412  0.6775 1.01  413   0.103
 time_3 -0.223 -0.404 -0.0493 1.01  714   0.104
```

`rhat` is the [Gelman-Rubin convergence diagnostic](https://www.rdocumentation.org/packages/coda/versions/0.19-3/topics/gelman.diag), `eff` is the [effective sample size](https://mc-stan.org/docs/2_18/reference-manual/effective-sample-size-section.html), and `ts_se` is the time-series standard error.

`plot(fit, "combo")` can be used to inspect the posteriors and convergence of all parameters. See the documentation of `plot.mcpfit` for many other plotting options. Here we plot just the (population-level) change points. They often have "strange" posterior distributions, highlighting the need for a computational approach:

```r
plot(fit, "combo", regex_pars = "cp_")
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_demo_combo.png)


For model comparisons, we can fit a null model and compare the predictive performance of the two models using (approximate) leave-one-out cross-validation. Our null model omits the first plateau and change point, essentially testing the credence of that change point:

```r
# Fit the model
segments_null = list(
  response ~ 1 + time,  # intercept (int_1) and slope (time_1)
  1 ~ 1 + time  # disjoined slope (int_2, time_1)
)
fit_null = mcp(segments_null, ex_demo)
```

Leveraging the power of `loo::loo`, we see that the two-change-points model is preferred (it is on top), but not very strongly (`elpd_diff / se_diff` ratio). See the documentation of the `loo` package and associated papers for more details, or scroll down for more notes on model comparison using `mcp`.

```r
# Compare loos:
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

[varying change points in mcp](https://lindeloev.github.io/mcp/articles/varying.html) include:
 * Simulate varying change points using `fit$func_y()`.
 * Get posteriors using `ranef(fit)`.
 * Plot using `plot(fit, facet_by="my_group")` and `plot(fit, "dens_overlay", pars = "varying", ncol = 3)`.
 * The default priors restrict varying change points to lie between the two adjacent change points.

`mcp` currently supports the following GLM:
 * `gaussian(link = "identity")` (default). See examples above and below.
 * `binomial(link = "logit")`. See [binomial change points in mcp](https://lindeloev.github.io/mcp/articles/binomial.html).
 * `bernoulli(link = "logit")`. See [binomial change points in mcp](https://lindeloev.github.io/mcp/articles/binomial.html).
 * `poisson(link = "log")`. See [Poisson change points in mcp](https://lindeloev.github.io/mcp/articles/poisson.html).



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
            name       mean       X2.5      X97.5     rhat  eff     ts_se
1 cp_1_id[Benny] -18.970413 -21.746111 -16.255543 1.000108 3131  5.317845
2  cp_1_id[Bill] -10.959035 -13.530281  -8.430385 1.002649  733 21.986801
3  cp_1_id[Cath]  -3.894166  -6.327470  -1.576214 1.000920 1156 21.315220
4  cp_1_id[Erin]   6.303391   3.462568   8.951041 1.000212 5684  2.104787
5  cp_1_id[John]   9.767637   7.223594  12.412474 1.000726 2089  7.184478
6  cp_1_id[Rose]  17.752586  14.684662  21.315577 1.001246 2704  4.979728
```


## Generalized linear models
`mcp` supports Generalized Linear Modeling. See extended examples using [`binomial()`](https://lindeloev.github.io/mcp/articles/binomial.html) and [`poisson()`](https://lindeloev.github.io/mcp/articles/poisson.html). [See how this data was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_binomial.R).

Here is a binomial change point model with three segments:

```r
segments = list(
  y | trials(N) ~ 1,  # constant rate
  1 ~ 0 + x,  # joined changing rate
  1 ~ 1 + x  # disjoined changing rate
)
fit = mcp(segments, ex_binomial, family = binomial())
plot(fit)
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
  int_1 = 10,  # fixed value
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
   y ~ 1 + x 
   y ~ 1 ~ rel(1) + rel(x) 
   y ~ rel(1) ~ 0 + x 

Population-level parameters:
  name  mean  X2.5 X97.5 rhat  eff    ts_se
  cp_1 29.52 29.05 30.00   NA 2115    0.373
  cp_2 49.33 45.81 52.51   NA  350   69.387
 int_1 10.00 10.00 10.00   NA    0    0.000
 int_2 33.91 24.56 43.75   NA  131 1066.654
 sigma  9.67  8.33 11.13   NA 4538    1.019
   x_1  1.64  1.42  1.87   NA   78    1.580
   x_2 -3.43 -3.73 -3.14   NA  172    1.308
   x_3  1.64  1.42  1.87   NA   78    1.580
```

(rhat cannot be estimated when the model contains fixed-to-value parameters in `mcp 0.1`. This will be fixed.)




# Diagnosing problems
**Convergence:** A common problem when using MCMC is lacking convergence between chains. This will show up as large `rhat` values (> 1.1 is a common criterion) and non-converging lines in `plot(fit, "trace")`. 

 * The first thing to try is always to make the model warm up longer to see if it reaches convergence later: `mcp(fit, data, adapt = 4000, update = 4000)`. 
 
 * It can be a sign of a deeper non-identifiability in the model. This will show up as strong correlations between parameters in the joint distribution of the implicated parameters: `plot(fit, "hex", pars = c("int_1", "int_2))`. This may give you ideas how to change the model.
 
 * You can set the initial values for the JAGS sampler using, e.g., `mcp(..., inits = list(cp_1 = 20, int_2 = -20, etc.))`. This will be passed to `jags.fit` and you can see more documentation there.


**Speed:** A lot of data and complicated models will slow down fitting. 

 * Run the chains in parallel using, e.g., `mcp(..., chains=4, cores=4)`.

 * More data usually means better identifiability and faster convergence. Lower the warmup period using, e.g., `mcp(..., adapt=300, update = 200)`.

**Won't run:** Most of these problems should stem from inappropriate priors and such problems may be exacerbated by fragile link functions (e.g., `binomial(link = "identity")`. The article on [priors in mcp](https://lindeloev.github.io/mcp/articles/priors.html) may be helpful, but in particular:

 * Errors on "directed cycle" usually stems from using parameters in priors. For example, if you set `prior = list(int_1 = "dnorm(int_2, 1)"", int_2 = "dnorm(int_1, 1)")` this is an infinite regress.
 
 * Errors on "incompatible with parent nodes" usually stem from impossible values. For example, if you set `prior = list(sigma = "dnorm(0, 1)"")`, this allows for a negative standard deviation, which is impossible. Think about changing the prior distributions and perhaps truncate them using `T(lower, upper)`.


If you encounter these or other problems, don't hesitate to [raise a Github Issue](https://github.com/lindeloev/mcp/issues), asking for help or posting a bug report.



# Notes on model comparison
We can use the cross-validation from the `loo` package to compare the predictive performance of various models. We compute loo for each model and then compare them. the `mcpfit` objects have placeholders at `fit$loo` and `fit$waic`. You need not use them, but it's nice to keep things together.

`loo::loo_compare`: The important thing is the `elpd_diff`/`se_diff` ratio. This is almost like a z-score, i.e., a ratio of 1.96 corresponds to 95% probability that `model1` (the first `loo` pass to `loo_compare`) has superior predictive accuracy. BUT, the `loo` developers have advised that it is probably too optimistic in practice, so that one should only begin taking it seriously if this ratio exceeds 5 [Citation needed].


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
