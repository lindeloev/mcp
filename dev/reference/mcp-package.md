# mcp: Multiple Change Point Regression in R

Flexible and informed regression with Multiple Change Points. `mcp` can
infer change points in models on means, variances (and other
distributional parameters), autocorrelation structure, and any
combination of these, as well as the parameters of the segments in
between. All parameters are estimated with uncertainty, and posterior
predictive intervals are supported - also near the change points. `mcp`
supports hypothesis testing via Savage-Dickey density ratios, posterior
contrasts, and PSIS-LOO/WAIC model comparison.

## Details

**The mcp model**

An `mcp` model divides a continuous predictor \\x\\ into \\K\\ segments
separated by ordered change points \\\Delta_1 \< \dots \<
\Delta\_{K-1}\\. In each segment \\k \in \\1, \dots, K\\\\, the linear
predictor \\\eta_i\\ is evaluated directly from the segment-local
distance \\(x_i - \Delta\_{k-1})\\:

\$\$\eta_i = \alpha_k + \beta\_{k,1} (x_i - \Delta\_{k-1}) \quad
(\text{with } \Delta_0 = 0)\$\$

where the segment-start level \\\alpha_k\\ is freely estimated for the
first and disjoined segments, as in non-segmented regression, and
inherited continuously for joined segments:

\$\$\alpha_k = \begin{cases} \beta\_{k,0}, & \text{Disjoined segments }
(\sim \texttt{1 + x}, \text{ including } k = 1) \\ \alpha\_{k-1} +
\beta\_{k-1,1} (\Delta\_{k-1} - \Delta\_{k-2}), & \text{Joined segments
} (k \ge 2, \sim \texttt{0 + x}) \end{cases}\$\$

In all segments, estimated slope and intercept parameters are absolute
values (not changes relative to the preceding segment).

If additional continuous covariates or categorical factors are included
(e.g., `+ z + group`), they enter additively on their original scale
(\\\dots + \sum \gamma\_{k,j} z\_{j,i}\\); only the change-point
predictor \\x\\ is converted to segment-local coordinates.

Distributional parameters
([`sigma()`](https://rdrr.io/r/stats/sigma.html), `shape()`, etc.) and
autoregressive terms ([`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`)
follow this exact same segmented structure on their respective link
scales.

**Extended formulas and model features**

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
  ([`loo()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md))
  or WAIC
  ([`waic()`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)),
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
