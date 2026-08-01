# Learnings
This files contains learnings that took a lot of work to reach and some of which meant ideas were abandoned, so it is not "documented" in code/articles.

## JAGS sensitivity to changepoint-parameter values
Even small rounding-to-7-decimal-points can catastrophically collapse JAGS sampling of cp-parameters.
See commit 9b887baacc018b5dd55016366e31e84ee548ed24 and docs for register_jags_constant() and jagsify_constants().

## Truncated-t vs dirichlet vs. stick breaking
Truncated-t yields 3-5 times ESS/second than dirichlet and stick-breaking (The Dirichlet(1,…,1) has a well-known stick-breaking representation using independent, untruncated Betas). So it is still kept default.

## Prediction quantile precision
It added ~100 lines of code to improce precision 4-fold per second using rao-blackwell. I opted not to do it.

## JAGS model optimization
There are no obvious optimizations for the currently generated JAGS code, in terms of increasing effective sample size per second.

## Parameter-plot pagination
`brms::plot.brmsfit()` uses five parameters per page by default. `plot_pars()`
follows its `nvariables` and `ask` conventions while retaining the existing
customizable ggplot return when only one page is needed. Multi-page calls return
every page in a list so no plots are silently discarded.

## Bridge sampling
Requires a lot of work on tracking priors, including their truncation etc. Since bridge sampling is very sensitive to priors anyway, and I doubt users will put a lot of thought into priors, I have opted not to do it.
