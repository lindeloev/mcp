# Package index

## Fit and show mcpfit

mcp() returns and `mcpfit`. These functions illuminate what this fit
contains structurally.

- [`mcp()`](https://lindeloev.github.io/mcp/reference/mcp.md) : Fit
  Multiple Linear Segments And Their Change Points

- [`mcp-package`](https://lindeloev.github.io/mcp/reference/mcp-package.md)
  : mcp: Multiple Change Point Regression in R

- [`mcp_pars()`](https://lindeloev.github.io/mcp/reference/mcp_pars.md)
  : Model parameters

- [`mcp_columns()`](https://lindeloev.github.io/mcp/reference/mcp_columns.md)
  : Model data columns

- [`prior_summary()`](https://lindeloev.github.io/mcp/reference/prior_summary.md)
  : Summarise priors used by an mcp model

- [`print(`*`<mcplist>`*`)`](https://lindeloev.github.io/mcp/reference/print.mcplist.md)
  : Print mcplist

- [`print(`*`<mcptext>`*`)`](https://lindeloev.github.io/mcp/reference/print.mcptext.md)
  : Nice Printing of Multiline Texts

- [`is.mcpfit()`](https://lindeloev.github.io/mcp/reference/is.mcpfit.md)
  :

  Checks if the Argument is an `mcpfit` Object

- [`family(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/model-accessors-mcpfit.md)
  [`nobs(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/model-accessors-mcpfit.md)
  [`model.frame(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/model-accessors-mcpfit.md)
  [`formula(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/model-accessors-mcpfit.md)
  :

  Extract Model Information from an `mcpfit`

- [`mcpfit-class`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  :

  Class `mcpfit` of Models Fitted with the mcp Package

## Visual summaries

Plotting of fits, parameters, and diagnostics.

- [`plot(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/plot.mcpfit.md)
  [`plot_dpar()`](https://lindeloev.github.io/mcp/reference/plot.mcpfit.md)
  : Plot full fits
- [`plot_pars()`](https://lindeloev.github.io/mcp/reference/plot_pars.md)
  : Plot individual parameters
- [`pp_check(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/pp_check.md)
  : Posterior Predictive Checks For Mcpfit Objects
- [`interpolate_newdata()`](https://lindeloev.github.io/mcp/reference/interpolate_newdata.md)
  : Returns a data.frame with all combos of predictors

## Textual summaries and statistics

Summaries that pertain to the whole `mcpfit` or to specific parameters.

- [`summary(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md)
  [`fixef(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md)
  [`ranef(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md)
  [`print(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md)
  : Summarise mcpfit objects

- [`vcov(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/posterior-uncertainty-mcpfit.md)
  [`confint(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/posterior-uncertainty-mcpfit.md)
  :

  Prior and Posterior Covariance and Central Intervals for `mcpfit`
  Objects

- [`hypothesis()`](https://lindeloev.github.io/mcp/reference/hypothesis.md)
  : Test Hypotheses Concerning Individual Parameters

- [`loo(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/loo.mcpfit.md)
  [`waic(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/loo.mcpfit.md)
  : Information Criteria for Model Comparison

## Extract draws or data

Per-observation summaries or draws

- [`predict(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  [`fitted(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  [`log_lik(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  [`residuals(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  :

  Fitted and predicted values of `mcp` models fits

- [`posterior_epred.mcpfit()`](https://lindeloev.github.io/mcp/reference/posterior_epred.mcpfit.md)
  [`posterior_predict.mcpfit()`](https://lindeloev.github.io/mcp/reference/posterior_epred.mcpfit.md)
  [`posterior_linpred.mcpfit()`](https://lindeloev.github.io/mcp/reference/posterior_epred.mcpfit.md)
  :

  Posterior prediction draws for `mcpfit` objects

- [`as_draws()`](https://lindeloev.github.io/mcp/reference/as_draws.mcpfit.md)
  [`as_draws_df()`](https://lindeloev.github.io/mcp/reference/as_draws.mcpfit.md)
  [`as_draws_array()`](https://lindeloev.github.io/mcp/reference/as_draws.mcpfit.md)
  [`as_draws_matrix()`](https://lindeloev.github.io/mcp/reference/as_draws.mcpfit.md)
  [`as_draws_rvars()`](https://lindeloev.github.io/mcp/reference/as_draws.mcpfit.md)
  :

  Extract MCMC Draws from `mcpfit` Objects

- [`ndraws()`](https://lindeloev.github.io/mcp/reference/draws-index-mcp.md)
  [`nchains()`](https://lindeloev.github.io/mcp/reference/draws-index-mcp.md)
  [`niterations()`](https://lindeloev.github.io/mcp/reference/draws-index-mcp.md)
  :

  Index `mcpfit` objects

## Families

Distributional families that are not available in base R.

- [`bernoulli()`](https://lindeloev.github.io/mcp/reference/bernoulli.md)
  : Bernoulli Family for mcp
- [`negbinomial()`](https://lindeloev.github.io/mcp/reference/negbinomial.md)
  : Negative Binomial for mcp
- [`mcpfamily()`](https://lindeloev.github.io/mcp/reference/mcpfamily.md)
  [`format(`*`<mcpfamily>`*`)`](https://lindeloev.github.io/mcp/reference/mcpfamily.md)
  [`print(`*`<mcpfamily>`*`)`](https://lindeloev.github.io/mcp/reference/mcpfamily.md)
  [`is.mcpfamily()`](https://lindeloev.github.io/mcp/reference/mcpfamily.md)
  : Create or Test Objects of Class "mcpfamily"

## Help and demos

Some showcases of typical mcp analyses. Most of these are discussed on
the [front page](https://lindeloev.github.io/mcp).

- [`mcp_example()`](https://lindeloev.github.io/mcp/reference/mcp_example.md)
  [`mcp_example_data()`](https://lindeloev.github.io/mcp/reference/mcp_example.md)
  : Run example models
- [`demo_fit`](https://lindeloev.github.io/mcp/reference/demo_fit.md) :
  Pre-fitted example mcp model

## Deprecated

Functions retained for backward compatibility.

- [`sd_to_prec()`](https://lindeloev.github.io/mcp/reference/sd_to_prec.md)
  : Transform an mcp prior to the parameterization used by JAGS.
