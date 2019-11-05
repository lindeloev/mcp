# mcp: Bayesian Inference of Multiple Change Points

[![Travis-CI status](https://img.shields.io/travis/lindeloev/mcp.svg)](https://travis-ci.org/lindeloev/mcp)
[![Coveralls status](https://codecov.io/gh/lindeloev/mcp/branch/master/graph/badge.svg)](https://coveralls.io/r/lindeloev/mcp)


Hierarchical Bayesian inference of Multiple Change Points (MCP) between Generalized Linear Segments - and the parameters in those segments. Under the hood, `mcp` takes an abstract representation of linear segments and turns it into [JAGS](https://sourceforge.net/projects/mcmc-jags/) code. The rest of the package ensures seamless compatibility with your favorite Bayesian packages, including `tidybayes`, `bayesplot`, `coda`, and `loo`.

You can see the roadmap for the immediate future under [issues](https://github.com/lindeloev/mcp/issues). Expect breaking changes in the API until version 1.0. All feedback is welcome, positive or negative or constructive. Please raise issues here or contact me elsewhere, e.g., [@jonaslindeloev](https://twitter.com/jonaslindeloev) on Twitter or via mail at firstname at lindeloev.dk.



# Install

 1. [Install the latest version of JAGS](https://sourceforge.net/projects/mcmc-jags/). Linux users can fetch binaries [here](http://mcmc-jags.sourceforge.net/).
 
 2. Install `mcp` by running this in R: `devtools::install_github("lindeloev/mcp")`. If you don't have `devtools`, install it first using `install.packages("devtools")`.





# Brief example
This model infers the change point (`cp_1`) on `time`, the intercept for segment 1 (`int_1`), and the intercept for segment 2 (`int_2`). Also the standard deviation of the residuals (`sigma`).


Find the two change points between three segments (three formulas):
```r
library(mcp)
segments = list(
  response ~ 1,  # plateau (int_1)
  1 ~ 0 + time,  # joined slope (time_2) at cp_1
  1 ~ 1 + time  # disjoined slope (int_1, time_2) at cp_2
)
fit = mcp(segments, data = ex_demo)  # dataset included in mcp
```

Use `summary` to see the estimates. They were simulated to lie at `cp_1 = 30` and `cp_2 = 70` ([see how `ex_demo` was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_demo.R)) and these values are well recovered, as are the other simulation parameters:

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

Plot lines drawn randomly from the posterior on top of data to inspect the fit:
```r
plot(fit)
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_demo.png)

Inspect the posteriors and convergence of the change points:

```r
plot(fit, "combo", regex_pars = "cp_")
```
![](https://github.com/lindeloev/mcp/raw/master/man/figures/ex_demo_combo.png)

Fit a null model with just one change point, consisting of two disjoined slopes:

```r
# Fit the model
segments_null = list(
  response ~ 1 + time,  # intercept (int_1) and slope (time_1)
  1 ~ 1 + time  # disjoined slope (int_2, time_1)
)
fit_null = mcp(segments_null, ex_demo)
```

Using leave-one-out cross-validation, we see that the two-change-points model is preferred (it is on top), but not very strongly. See more about this in the `loo` package.

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
The article on [understanding mcp formulas and models](https://lindeloev.github.io/mcp/articles/formulas.html):
 * Parameter names `int_i` (intercepts), `cp_i` (change points), `x_i` (slopes).
 * The change point model is basically an `ifelse` model.
 * Use `rel()` to specify that parameters are relative to those in the previous  segments.

The article on [using priors in mcp](https://lindeloev.github.io/mcp/articles/priors.html) include:
 * See priors in `fit$prior`
 * Set priors using `mcp(segments, data, prior = list(cp_1 = "dnorm(0, 1)", cp_1 = "dunif(0, 45)`)
 * Fix parameters to specific values using `cp_1 = 45`
 * Share parameters between segments using `slope_1 = "slope_2"`
 * Allows for truncated priors using `T(lower, upper)`, e.g., `int_1 = "dnorm(0, 1) T(0, )"`. `mcp` applies this automatically to change point priors to enforce order restriction. This is true for [varying change points](https://lindeloev.github.io/mcp/articles/varying.html) too.
 * Sample priors for prior predictive checks using `mcp(segments, data, sample="prior")`.

The article on [varying change points in mcp](https://lindeloev.github.io/mcp/articles/varying.html) include:
 * How to simulate varying change points using `fit$func_y()`.
 * Get posteriors using `ranef(fit)`
 * Plot using `plot(fit, facet_by="my_group")` and `plot(fit, "dens_overlay", pars = "varying", ncol = 3)`.
 * The default priors restrict varying change points to lie between the two adjecent change points.

`mcp` currently supports the following GLM:
 * `gaussian(link = "identity")`. See worked example below.
 * `binomial(link = "logit")`. See [binomial change points in mcp](https://lindeloev.github.io/mcp/articles/binomial.html).
 * `bernoulli(link = "logit")`. See [binomial change points in mcp](https://lindeloev.github.io/mcp/articles/binomial.html).
 * `poisson(link = "log")`. See [Poisson change points in mcp](https://lindeloev.github.io/mcp/articles/poisson.html).


# More examples


## Two plateaus
Find the single change point between two plateaus ([see how `ex_plateaus` was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_plateaus.R))

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

Here, we find the single change point between two joined slopes. While the slopes are shared by all participants (they are population-level), their change point varies.

Read more about [varying change points](https://lindeloev.github.io/mcp/articles/varying.html) and [see how ex_varying was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_plateaus.R).

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

```
            name       mean       X2.5      X97.5     rhat  eff     ts_se
1 cp_1_id[Benny] -18.970413 -21.746111 -16.255543 1.000108 3131  5.317845
2  cp_1_id[Bill] -10.959035 -13.530281  -8.430385 1.002649  733 21.986801
3  cp_1_id[Cath]  -3.894166  -6.327470  -1.576214 1.000920 1156 21.315220
4  cp_1_id[Erin]   6.303391   3.462568   8.951041 1.000212 5684  2.104787
5  cp_1_id[John]   9.767637   7.223594  12.412474 1.000726 2089  7.184478
6  cp_1_id[Rose]  17.752586  14.684662  21.315577 1.001246 2704  4.979728
```


## Generalized linear models
`mcp` supports Generalized Linear Modeling. See extended examples using [`binomial()`](https://lindeloev.github.io/mcp/articles/binomial.html) and [`poisson()`](https://lindeloev.github.io/mcp/articles/poisson.html). [See how `ex_binomial` was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_binomial.R).

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
Read more about [formula options](https://lindeloev.github.io/mcp/articles/formulas.html) and [priors](https://lindeloev.github.io/mcp/articles/priors.html) and [see how `ex_rel_prior` was generated](https://github.com/lindeloev/mcp/tree/master/data-raw/ex_rel_prior.R).

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

```
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
**Convergence:** A common problem when using MCMC is lacking convergence between chains. This will show up as large `rhat` values (> 1.1 is a common criterion) and non-converging lines in `plot(fit, "trace")`. The first thing to try is always to make the model warm up longer to see if it reaches convergence later: `mcp(fit, data, adapt = 4000, update = 4000)`. It can be a sign of a deeper non-identifiability in the model. This will show up as strong correlations between parameters in the joint distribution of the implicated parameters: `plot(fit, "hex", pars = c("int_1", "int_2))`.

**Speed:** A lot of data will slow down fitting. You can try combinations of lowerering the warmup period and running the chains in parallel: `mcp(fit, data, chains = 4, cores = 4, adapt=300, update = 200)`. Convergence issues will also slow down sampling.

If you encounter other problems, don't hesitate to [raise a Github Issue](https://github.com/lindeloev/mcp/issues), asking for help or posting a bug report.



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
