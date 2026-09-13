# mcp: Multiple Change Point Regression in R

Flexible and informed regression with Multiple Change Points. `mcp` can
infer change points in models on means, variances (and other
distributional parameters), autocorrelation structure, and any
combination of these, as well as the parameters of the segments in
between. All parameters are estimated with uncertainty, and posterior
predictive intervals are supported - also near the change points. `mcp`
supports hypothesis testing via Savage-Dickey density ratios, posterior
contrasts, and PSIS-LOO/WAIC model comparison.

## Extended formulas and model features

- *Distributional regression:* Model residual variance and other
  distributional parameters across segments, e.g., `~ sigma(1 + x)` or
  `~ shape(1)`.

- *Time-series residuals:* Model serial dependence using `ar(p)` and
  `ma(q)` terms with a generalized link-scale recurrence that spans
  continuously across change points.

- *Group-level effects:* Add hierarchical (random) intercepts, slopes,
  and change points, e.g., `~ 1 + (1|id)` or `1 + (1|id) ~ 0 + x`.

- *Model comparison and hypothesis testing:* Compare models using
  PSIS-LOO
  ([`loo()`](https://lindeloev.github.io/mcp/reference/loo.mcpfit.md))
  or WAIC
  ([`waic()`](https://lindeloev.github.io/mcp/reference/loo.mcpfit.md)),
  and test hypotheses
  ([`hypothesis()`](https://lindeloev.github.io/mcp/reference/hypothesis.md))
  using Savage-Dickey density ratios for point-null tests or
  posterior-to-prior odds for directional tests.

See [the mcp website](https://lindeloev.github.io/mcp/) for worked
examples.

## The mcp model

Consider the following model which you can find in `demo_fit` and
`mcp_example("demo")`:

    model = list(
      response ~ 1,  # Plateau in the first segment (Intercept_1)
      ~ 0 + time,    # Joined slope (time_2) in segment 2 which starts at cp_1
      ~ 1 + time     # Disjoined slope (Intercept_3, time_3) at cp_2
    )

![Fitted 3-segment mcp model with a plateau, joined slope, and disjoined
slope](figures/mcp_demo.png)

This model has \\K=3\\ segments separated by \\K-1=2\\ change points:
\\\tau_1\\ and \\\tau_2\\.

More generally, an `mcp` model divides a continuous predictor \\x\\ into
\\K\\ segments separated by ordered change points \\\tau_1 \< \dots \<
\tau\_{K-1}\\. In each segment \\k \in \\1, \dots, K\\\\, the linear
predictor \\\eta_i\\ is evaluated directly from the segment-local
distance \\(x_i - \tau\_{k-1})\\:

\$\$\eta_i = \alpha_k + \beta\_{k,1} (x_i - \tau\_{k-1}) \quad
(\text{with } \tau_0 = 0)\$\$

where the segment-start level \\\alpha_k\\ is freely estimated for the
first and disjoined segments, as in non-segmented regression, and
determined by continuity for joined segments:

\$\$\alpha_k = \begin{cases} \beta\_{k,0}, & \text{Disjoined segments }
(\sim \texttt{1 + x}, \text{ including } k = 1) \\ \alpha\_{k-1} +
\beta\_{k-1,1} (\tau\_{k-1} - \tau\_{k-2}), & \text{Joined segments } (k
\ge 2, \sim \texttt{0 + x}) \end{cases}\$\$

Here, \\\beta\_{k,0}\\ is the segment-start intercept, and
\\\beta\_{k,1}\\ is the slope on \\x\\. In all segments, estimated slope
and intercept parameters are absolute values (not changes relative to
the preceding segment).

If additional continuous covariates or categorical factors are included
(e.g., `+ z + group`), they enter additively on their original scale
(\\\dots + \sum \gamma\_{k,j} z\_{j,i}\\ for covariate \\j\\); only bare
change-point predictor terms \\x\\ and `I(x^k)` are converted to
segment-local coordinates.

The inverse-linked parameter is \\\mu_i = g^{-1}(\eta_i)\\ via link
function \\g(\mu_i) = \eta_i\\, representing the expected response for
most families (or the success probability for binomial models, where the
expected count is \\n_i \mu_i\\). Distributional parameters
([`sigma()`](https://rdrr.io/r/stats/sigma.html), `shape()`, etc.) and
autoregressive terms ([`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`)
follow this exact same segmented structure on their respective link
scales. See more details on the `mcp` model in mcp-package and on the
[mcp website](https://lindeloev.github.io/mcp/articles/formulas.html).

## Time-series residuals (link-scale observation-driven GARMA)

Autoregressive (`ar(p)`) and moving-average (`ma(q)`) terms define a
finite conditional recurrence on the link scale (generalized
autoregressive moving-average, GARMA). They support Gaussian
(`identity`), binomial (`logit`), Bernoulli (`logit`), Poisson (`log`),
and negative-binomial (`log`) families. If \\\eta^{\text{reg}}\_t\\ is
the ordinary regression predictor from the segment formulas and
\\\eta_t\\ is the predictor including serial dependence, the recurrence
decomposes into components:

\$\$\begin{aligned} \text{AR}\_t &= \sum\_{j=1}^{p} \phi\_{j,t}
\left\[g(y^\*\_{t-j}) - \eta^{\text{reg}}\_{t-j}\right\] \\ \text{MA}\_t
&= \sum\_{k=1}^{q} \theta\_{k,t} \left\[g(y^\*\_{t-k}) -
\eta\_{t-k}\right\] \\ \eta_t &= \eta^{\text{reg}}\_t + \text{AR}\_t +
\text{MA}\_t \end{aligned}\$\$

where \\\phi\_{j,t}\\ is the lag-\\j\\ autoregressive (AR) coefficient
at time \\t\\, \\\theta\_{k,t}\\ is the lag-\\k\\ moving-average (MA)
coefficient at time \\t\\, \\g(\cdot)\\ is the link function, and
\\y^\*\_t\\ is the boundary-constrained observation with pseudo-count
\\b\\ (set via argument `boundary = 0.1` in
[`ar()`](https://rdrr.io/r/stats/ar.html) / `ma()`) to keep residuals
finite on the link scale:

- **Gaussian:** \\y^\*\_t = y_t\\.

- **Poisson / Negative Binomial:** \\y^\*\_t = \max(y_t, b)\\ to prevent
  \\\log(0)\\. Here \\b\\ replaces zero counts with a small positive
  count.

- **Binomial / Bernoulli:** \\y^\*\_t = \min(\max(y_t, b), n_t - b) /
  n_t\\, where \\y_t\\ is observed successes, \\n_t\\ is the number of
  trials (\\n_t = 1\\ for Bernoulli), and \\b\\ clamps counts to the
  interval \\\[b, n_t - b\]\\ before converting to a rate, preventing
  \\\text{logit}(0)\\ and \\\text{logit}(1)\\.

Implications:

- For an \\N\\-order component, the last \\N\\ values *before* the
  segment onset are input to the first \\\eta_t\\ in the segment.

- AR and MA components persist into later segments until replaced or
  turned off via `ar(0)` or `ma(0)`.

- AR coefficients are not jointly constrained to stationarity; nor MA
  coefficients to invertibility.

- See [the arma
  article](https://lindeloev.github.io/mcp/articles/arma.html) for more
  details.

## References

- Lindeløv, J. K. (2020). mcp: An R Package for Regression With Multiple
  Change Points. *OSF Preprints*.
  [doi:10.31219/osf.io/fzqxv](https://doi.org/10.31219/osf.io/fzqxv)
  Introduces the `mcp` package, formula syntax, default priors, and
  workflow for regression with multiple change points across generalized
  linear and time-series models. Newer-than-2020 versions of the paper
  may be available at that link. Please cite the newest version.

- Carlin, B. P., Gelfand, A. E., & Smith, A. F. (1992). Hierarchical
  Bayesian Analysis of Changepoint Problems. *Applied Statistics*,
  41(2), 389–405. [doi:10.2307/2347570](https://doi.org/10.2307/2347570)
  Introduced the Gibbs sampling approach for continuous change points in
  segmented regression and hierarchical models, providing the
  computational foundation used by BUGS and JAGS.

## See also

[`mcp`](https://lindeloev.github.io/mcp/reference/mcp.md)

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
