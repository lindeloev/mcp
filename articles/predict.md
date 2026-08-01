# Fits and predictions

This article introduces
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
[`residuals()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
for in-sample and out-of-sample data. I will also show how to get
creative with `mcp`, including how to make predictions around future
change points.

## Preparation: an example model

We need an `mcpfit` to get started. We take the “demo” dataset:

``` r

library(mcp)
future::plan(future::multisession, workers = 3)
set.seed(42)  # Make the script deterministic

data = mcp_example_data("demo")
head(data)
```

    ##     response     time
    ## 1 -3.0084198 91.48060
    ## 2 -7.8768640 93.70754
    ## 3 16.3029101 28.61395
    ## 4 -0.0373553 83.04476
    ## 5 27.4463185 64.17455
    ## 6 22.0610004 51.90959

… and model it as three segments, i.e., two change points:

``` r

# Define the model
model = list(
  response ~ 1,  # plateau (Intercept_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1 + time     # disjoined slope (Intercept_3, time_3) at cp_2
)

# Fit it
fit = mcp(model, data = data)
```

    ## Warning: MultisessionFuture ('future_lapply-1') added, removed, or modified
    ## connections. A future expression must close any opened connections and must not
    ## close connections it did not open. Details: 1 connection added ([index=6,
    ## description=jags_code, class=textConnection, mode=r, text=text, opened=opened,
    ## can.read=yes, can.write=no]), 0 connection removed (<none>), 0 connection
    ## replaced (<none>). See also help("future.options", package = "future") [future
    ## 'future_lapply-1' (c53bd8ca5dc1e6bdfe82a518e9ff3d77-1); on
    ## c53bd8ca5dc1e6bdfe82a518e9ff3d77@runnervmvrwv9<8560>]

    ## Warning: MultisessionFuture ('future_lapply-2') added, removed, or modified
    ## connections. A future expression must close any opened connections and must not
    ## close connections it did not open. Details: 1 connection added ([index=6,
    ## description=jags_code, class=textConnection, mode=r, text=text, opened=opened,
    ## can.read=yes, can.write=no]), 0 connection removed (<none>), 0 connection
    ## replaced (<none>). See also help("future.options", package = "future") [future
    ## 'future_lapply-2' (c53bd8ca5dc1e6bdfe82a518e9ff3d77-2); on
    ## c53bd8ca5dc1e6bdfe82a518e9ff3d77@runnervmvrwv9<8560>]

    ## Warning: MultisessionFuture ('future_lapply-3') added, removed, or modified
    ## connections. A future expression must close any opened connections and must not
    ## close connections it did not open. Details: 1 connection added ([index=6,
    ## description=jags_code, class=textConnection, mode=r, text=text, opened=opened,
    ## can.read=yes, can.write=no]), 0 connection removed (<none>), 0 connection
    ## replaced (<none>). See also help("future.options", package = "future") [future
    ## 'future_lapply-3' (c53bd8ca5dc1e6bdfe82a518e9ff3d77-3); on
    ## c53bd8ca5dc1e6bdfe82a518e9ff3d77@runnervmvrwv9<8560>]

    ## Warning: Some parameters may not have converged well:
    ##   * ess_bulk or ess_tail < 400: cp_1 and time_2
    ## Inspect `summary(fit)` and `plot_pars(fit)`, and consider increasing `iter`/`adapt` or simplifying the model before trusting these results.

This is what the data and the inferred fit looks like with 95% credible
interval and a 80% prediction interval:

``` r

plot(fit, q_fit = TRUE, q_predict = c(0.1, 0.9))
```

![](predict_files/figure-html/unnamed-chunk-3-1.png)

To review what we see here:

- The black dots is the data from `data`.
- The gray lines are 25 samples from the posterior (control using
  `plot(fit, lines = 100)`).
- The dashed red lines are the 2.5% and 97.% quantiles of the fitted
  (expected) values.
- The green lines are the 10% and 90% quantiles of the predicted values.
- The blue curves on the x-axis are the posterior distributions of the
  change point locations (better viewed using
  `plot_pars(fit, pars = c("cp_1", "cp_2"))`).

Behind the scenes,
[`plot()`](https://lindeloev.github.io/mcp/reference/plot.mcpfit.md)
merely calls
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
and
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
to show these inferences.

## Extracting `fitted()` values

To get the fitted values for each data point, simply do `fitted(fit)`:

``` r

head(fitted(fit))
```

    ##       time    fitted     error      Q2.5       Q97.5
    ## 1 91.48060 -3.590968 0.7388680 -5.044832 -2.12438624
    ## 2 93.70754 -4.195168 0.8256646 -5.820004 -2.54570655
    ## 3 28.61395 10.984224 1.0622749  9.001740 12.96695246
    ## 4 83.04476 -1.302205 0.6528397 -2.564974 -0.01641898
    ## 5 64.17455 25.069951 0.8849647 23.366532 26.80580830
    ## 6 51.90959 20.143706 0.6052964 18.951794 21.32189370

In general, this output will include:

- A column for each predictor column in the data. Here, `time` is the
  only predictor and you see the values in the same order as in `data`
  (which is copied to `fit$data`). Models with [group-level
  effects](https://lindeloev.github.io/mcp/articles/varying.md)
  additionally include the relevant grouping columns,
  [`binomial()`](https://rdrr.io/r/stats/family.html) models include the
  number of trials, etc.

- **fitted:** The fitted value (posterior mean). When `summary = TRUE`
  (default), the column is called `fitted`. When `summary = FALSE`, the
  column is named `.epred`, matching
  [`tidybayes::add_epred_draws()`](https://mjskay.github.io/tidybayes/reference/add_predicted_draws.html)
  conventions.

- **error:**: The standard error corresponding to `fitted`, i.e.,
  `diff(fitted + c(-1, 1) * error)` is the 68% credible interval.

- **Q\[some number\]:** The quantiles of the fitted distribution. You
  can set the quantiles using `fitted(fit, probs = c(0.1, 0.5, 0.9))`.

If you compare these values to the plot, you will see that they
correspond.
[`plot()`](https://lindeloev.github.io/mcp/reference/plot.mcpfit.md)
merely calls
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
behind the scenes.

### Expected values for out-of-sample data

To predict out-of-sample data, you can simply use the `newdata`
argument.

``` r

newdata = data.frame(time = c(data$time[1], 25, -20, 200))
fitted(fit, newdata = newdata)
```

    ##       time     fitted     error       Q2.5      Q97.5
    ## 1  91.4806  -3.590968 0.7388680  -5.044832  -2.124386
    ## 2  25.0000  10.003976 0.8637372   8.439159  11.798917
    ## 3 -20.0000   8.986798 0.9194766   7.108671  10.667567
    ## 4 200.0000 -33.033820 7.6010447 -48.073497 -18.580933

Note that:

- We get one row per value of `time`.
- The first value for `time` is in the dataset. The values correspond to
  the same row in `fitted(fit)` because that’s merely a shortcut to do
  `fitted(fit, newdata = fit$data)`.
- The second value (`time = 20`) is within the observed region, but not
  in the dataset.
- The third value (`time = - 20`) is outside the observed region, but
  `mcp` merely extends the first segment backwards in time. Because it’s
  a plateau, we see approximately the same values as for `time = 20`.
- The fourth value (`time = 200`) is way outside the observed region.
  Because it is the extrapolation of the slope in the third segment of
  which we’ve only observed the first tiny bit, the posterior
  distribution is very wide because even a small uncertainty in the
  slope results in very large differences further out.

### Arguments

If you look at the documentation for
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
and
[`residuals()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
you’ll see that they are quite versatile, taking many different
arguments. To mention a few, you can set `fitted(fit, dpar = "sigma")`
to get fitted values for `sigma` [more on modeling
sigma](https://lindeloev.github.io/mcp/articles/variance.md),
`prior = TRUE` to predict using only the prior, and `arma = FALSE` to
exclude AR/MA effects. For group-level effects, use `varying = TRUE`
(all), `FALSE` (none), `"cp"` or `"predictor"` (a formula part), or an
exact group-level parameter name.

## Predictions

[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
is the posterior predictive and it takes exactly the same arguments as
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md).
This means that you can make predictions for in-sample and out-of-sample
data as well. As with
[`fitted()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md),
[`plot()`](https://lindeloev.github.io/mcp/reference/plot.mcpfit.md)
uses
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
under the hood to plot prediction intervals. You can see that the values
correspond to dashed green lines in the plot (the 80% prediction
interval):

``` r

set.seed(42)
head(predict(fit, probs = c(0.1, 0.9)))
```

    ##       time   predict    error       Q10        Q90
    ## 1 91.48060 -3.662650 3.765566 -8.395376  1.1901654
    ## 2 93.70754 -4.218413 3.778896 -9.017438  0.6050123
    ## 3 28.61395 11.036668 3.818941  6.154950 15.8488343
    ## 4 83.04476 -1.333384 3.730039 -6.023077  3.4001579
    ## 5 64.17455 24.966597 3.783374 20.109031 29.8503995
    ## 6 51.90959 20.154743 3.697339 15.410159 24.8859457

Note that
[`predict()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
uses random sampling under the hood, so these values will differ
slightly from call to call. You can make it replicable using
[`set.seed()`](https://rdrr.io/r/base/Random.html) as above. In general,
the more posterior draws, the less the call-to-call variance will be.
Conversely, fewer draws means more call-to-call variation, e.g., if you
do `predict(fit, ndraws = 10)`).

## Residuals

[`residuals()`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
is simply `data$response - fitted()`. It may be useful for model
checking, but the typical needs are covered using posterior predictive
checking (`pp_check(fit)`) and visual inspection of
`plot(fit, q_fit = TRUE, q_predict = TRUE)`.

## Forecasting with future change points

Bayesian inference is the principled updating of prior knowledge using
data. Where there is little or now data, the prior speaks louder.
Sometimes, we can learn surprising stuff simply by inspecting the prior
predictive, e.g., how the priors combine when “put through” the model.
In `mcp`, most functions come with a `prior = FALSE` default, but you
can simply do `plot(fit, prior = TRUE)`, `fitted(fit, prior = TRUE)`, or
`predict(fit, prior = TRUE)`.

Say you want to forecast at `time = 125` and you know that a changepoint
to the baseline level (an intercept change to `Intercept_1`) will occur
approximately after the same interval as between `cp_1` and `cp_2`
(i.e., at `cp_2 + (cp_2 - cp_1)`). Here is a way to “hack” `mcp` to do
this. *(NOTE: I plan on implementing this in a much more user-friendly
way in a future release; see the discussion in [this github
issue](https://github.com/lindeloev/mcp/issues/78) and current status in
[this github issue](https://github.com/lindeloev/mcp/issues/82))*.

### Step 1: run the model for observed data

We already did that above, resulting in our `fit`. But we only do it to
get the default priors that are suitable for inferring change point in
this region, so you could’ve just run it without sampling:

``` r

fit = mcp(model, data = data, sample = FALSE)
```

### Step 2: add the unobserved segment(s) and fit

Now we extend the model with the future segment of which we have prior
knowledge:

``` r

model_forecast = c(fit$model, list(
  ~ 1     # intercept (Intercept_4) after cp_3
))
```

And finally, we extend the list of priors with the two new parameters
(`time_4` and `cp_3`). It may be helpful to review [the article on
priors in mcp](https://lindeloev.github.io/mcp/articles/priors.md).

``` r

prior_forecast = c(fit$prior, list(
  Intercept_4 = "Intercept_1",  # Return to this value
  cp_3 = "dnorm(cp_2 + (cp_2 - cp_1), 20) T(max(time), )"  # In the future at the same interval
))
```

Now let’s fit it:

``` r

fit_forecast = mcp(model_forecast, data = data, prior = prior_forecast)
```

    ## Warning: MultisessionFuture ('future_lapply-1') added, removed, or modified
    ## connections. A future expression must close any opened connections and must not
    ## close connections it did not open. Details: 1 connection added ([index=7,
    ## description=jags_code, class=textConnection, mode=r, text=text, opened=opened,
    ## can.read=yes, can.write=no]), 0 connection removed (<none>), 0 connection
    ## replaced (<none>). See also help("future.options", package = "future") [future
    ## 'future_lapply-1' (c53bd8ca5dc1e6bdfe82a518e9ff3d77-4); on
    ## c53bd8ca5dc1e6bdfe82a518e9ff3d77@runnervmvrwv9<8560>]

    ## Warning: MultisessionFuture ('future_lapply-2') added, removed, or modified
    ## connections. A future expression must close any opened connections and must not
    ## close connections it did not open. Details: 1 connection added ([index=7,
    ## description=jags_code, class=textConnection, mode=r, text=text, opened=opened,
    ## can.read=yes, can.write=no]), 0 connection removed (<none>), 0 connection
    ## replaced (<none>). See also help("future.options", package = "future") [future
    ## 'future_lapply-2' (c53bd8ca5dc1e6bdfe82a518e9ff3d77-5); on
    ## c53bd8ca5dc1e6bdfe82a518e9ff3d77@runnervmvrwv9<8560>]
    ## Warning: MultisessionFuture ('future_lapply-2') added, removed, or modified
    ## connections. A future expression must close any opened connections and must not
    ## close connections it did not open. Details: 1 connection added ([index=7,
    ## description=jags_code, class=textConnection, mode=r, text=text, opened=opened,
    ## can.read=yes, can.write=no]), 0 connection removed (<none>), 0 connection
    ## replaced (<none>). See also help("future.options", package = "future") [future
    ## 'future_lapply-2' (c53bd8ca5dc1e6bdfe82a518e9ff3d77-5); on
    ## c53bd8ca5dc1e6bdfe82a518e9ff3d77@runnervmvrwv9<8560>]

    ## Warning: MultisessionFuture ('future_lapply-3') added, removed, or modified
    ## connections. A future expression must close any opened connections and must not
    ## close connections it did not open. Details: 1 connection added ([index=7,
    ## description=jags_code, class=textConnection, mode=r, text=text, opened=opened,
    ## can.read=yes, can.write=no]), 0 connection removed (<none>), 0 connection
    ## replaced (<none>). See also help("future.options", package = "future") [future
    ## 'future_lapply-3' (c53bd8ca5dc1e6bdfe82a518e9ff3d77-6); on
    ## c53bd8ca5dc1e6bdfe82a518e9ff3d77@runnervmvrwv9<8560>]
    ## Warning: MultisessionFuture ('future_lapply-3') added, removed, or modified
    ## connections. A future expression must close any opened connections and must not
    ## close connections it did not open. Details: 1 connection added ([index=7,
    ## description=jags_code, class=textConnection, mode=r, text=text, opened=opened,
    ## can.read=yes, can.write=no]), 0 connection removed (<none>), 0 connection
    ## replaced (<none>). See also help("future.options", package = "future") [future
    ## 'future_lapply-3' (c53bd8ca5dc1e6bdfe82a518e9ff3d77-6); on
    ## c53bd8ca5dc1e6bdfe82a518e9ff3d77@runnervmvrwv9<8560>]

    ## Warning: Some parameters may not have converged well:
    ##   * ess_bulk or ess_tail < 400: cp_1
    ## Inspect `summary(fit)` and `plot_pars(fit)`, and consider increasing `iter`/`adapt` or simplifying the model before trusting these results.

### Step 3: predict!

We can go right ahead and compute our 50% and 80% prediction intervals
at `time = 125`:

``` r

predict(fit_forecast, newdata = data.frame(time = 125), probs = c(0.1, 0.25, 0.75, 0.9))
```

    ##   time    predict    error       Q10       Q25      Q75      Q90
    ## 1  125 0.04553432 11.47093 -15.85141 -11.52191 9.717135 12.58547

To really understand what’s going on here, it may be helpful to
visualize the model. For now, we will have to hack this a bit too,
manually doing our plot:

``` r

# Get posterior and posterior predictive "predictions"
newdata = data.frame(time = 1:170)
fitted_forecast = fitted(fit_forecast, newdata = newdata, summary = FALSE, ndraws = 50)
predict_forecast = predict(fit_forecast, newdata = newdata, summary = FALSE)

# Plot it
library(ggplot2)
ggplot(predict_forecast, aes(x = time, y = .prediction)) +
  # Prediction intervals and line at x = 125
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.8), geom = "ribbon", alpha = 0.2) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), geom = "ribbon", alpha = 0.3) +
  geom_vline(xintercept = 125, lty = 2, lwd = 1) +

  # Lines for fitted draws
  geom_line(aes(y = .epred, group = .draw), data = fitted_forecast, alpha = 0.2) +

  # Observed data
  geom_point(aes(x = time, y = response), data = data) +
  labs(title = "Predicting with future change points")
```

    ## Warning: Computation failed in `stat_summary()`.
    ## Caused by error in `fun.data()`:
    ## ! The package "Hmisc" is required.

    ## Warning: Computation failed in `stat_summary()`.
    ## Caused by error in `fun.data()`:
    ## ! The package "Hmisc" is required.

![](predict_files/figure-html/unnamed-chunk-12-1.png)

You can read the predicted values from above at `x = 125` off this
graph. We literally just predicted for all values between 1 and 170, and
visualized it using a ribbon. This means that you can also predict
further into the future, if you’d like.

You can extend this approach to an arbitrary number of future segments,
even using the posterior from the “unobserved” segment 4 in the priors
for parameters in future segments. In Bayesian inference, it really does
not make much of a difference whether credence in some parameter values
have been updated using data or not - it’s all credence.

Without doing this formal model of the future change point, one may have
thought that the change point should occur around `time = 110` since
that’s the expected value of `cp_2 + (cp_2 - cp_1)`. However, we
truncated the prior for the future change point (`cp_3`) so that it
occurs *after* the last data point (`max(time)`), i.e., at `time > 100`.
This is knowledge that the third change point had not yet been observed
at `time = 100`, and this pushes the distribution further into the
future (actually around 118; see `summary(fit_forecast)`).
