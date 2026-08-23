# mcp: Multiple Change Point Regression in R

Flexible and informed regression with Multiple Change Points. `mcp` can
infer change points in models on means, variances (and other
distributional parameters), autocorrelation structure, and any
combination of these, as well as the parameters of the segments in
between. All parameters are estimated with uncertainty, and prediction
intervals are supported - also near the change points. `mcp` supports
hypothesis testing via Savage-Dickey density ratios, posterior
contrasts, and PSIS-LOO/WAIC model comparison.

## Details

**The mcp model**

For a continuous predictor \\x\\, segment 1 (\\x \le \Delta_1\\) is a
standard regression model with the intercept defined at \\x = 0\\:

\$\$\eta\_{1, i} = f_1(x_i, \mathbf{\beta}\_1) = \beta\_{1,0} +
\beta\_{1,1} x_i + \dots\$\$

For subsequent segments \\k \in \\2, \dots, K\\\\, each change point
\\\Delta\_{k-1}\\ marks the onset of a new segment. The model defines a
segment-local predictor \\X\_{k,i}\\ measured from the segment onset
\\\Delta\_{k-1}\\ and capped at the segment end \\\Delta_k\\:

\$\$X\_{k,i} = \min(x_i, \Delta_k) - \Delta\_{k-1}\$\$

with \\\Delta_K = \max(\mathbf{x})\\. The local predictor \\X\_{k,i}\\
is zero at segment onset (\\x_i \le \Delta\_{k-1}\\) and plateaus at
\\\Delta_k - \Delta\_{k-1}\\ for \\x_i \ge \Delta_k\\.

**Joined (continuous) segments (`~ 0 + x`)**

In joined models, segment \\k \ge 2\\ has no intercept and specifies a
new segment slope \\\beta\_{k,1}\\:

\$\$\eta_i = \eta\_{1,i}(X\_{1,i}) + \sum\_{k=2}^K \[x_i \ge
\Delta\_{k-1}\] \cdot \beta\_{k,1} X\_{k,i}\$\$

where \\X\_{1,i} = \min(x_i, \Delta_1)\\. Because earlier segments
plateau at their boundary values, the expectation is automatically
continuous across change points without requiring parameter constraints.
In both joined and disjoined segments, slope and intercept parameters
are absolute values (not differences relative to the previous segment).

**Disjoined (discontinuous) segments (`~ 1 + x`)**

When segment \\k\\ contains an explicit intercept (\\\beta\_{k,0}\\), it
introduces a step change at \\x_i = \Delta\_{k-1}\\ and previous
segments are truncated at \\\[x_i \< \Delta\_{k-1}\]\\.

**Extended formulas and model features**

- *Distributional regression:* Model residual variance and other
  distributional parameters across segments, e.g., `~ sigma(1 + x)` or
  `~ shape(1)`.

- *Time-series residuals:* Model serial dependence using `ar(p)` and
  `ma(q)` terms with a generalized link-scale recurrence that spans
  continuously across change points.

- *Group-level effects:* Add hierarchical (random) intercepts, slopes,
  and change points, e.g., `~ 1 + (1|id)` or `1 + (1|id) ~ 0 + x`.

- *Model comparison and hypothesis testing:* Compare models using LOO-CV
  ([`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md),
  [`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md))
  and test directional or point-null hypotheses
  ([`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md))
  using Savage-Dickey density ratios.

See [the mcp website](https://lindeloev.github.io/mcp/) for worked
examples and vignettes.

## References

- Lindeløv, J. K. (2020). mcp: An R Package for Regression With Multiple
  Change Points. *OSF Preprints*.
  [doi:10.31219/osf.io/fzqxv](https://doi.org/10.31219/osf.io/fzqxv)

- Carlin, B. P., Gelfand, A. E., & Smith, A. F. (1992). Hierarchical
  Bayesian Analysis of Changepoint Problems. *Applied Statistics*,
  41(2), 389–405. [doi:10.2307/2347570](https://doi.org/10.2307/2347570)

- Stephens, D. A. (1994). Bayesian retrospective multiple-changepoint
  identification. *Applied Statistics*, 43(1), 159–178.
  [doi:10.2307/2986119](https://doi.org/10.2307/2986119)

## See also

[`mcp`](https://lindeloev.github.io/mcp/dev/reference/mcp.md)

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
