# Package index

## Using mcp

Functions for everyday use of mcp.

- [`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) : Fit
  Multiple Linear Segments And Their Change Points

- [`plot(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
  [`plot_dpar()`](https://lindeloev.github.io/mcp/dev/reference/plot.mcpfit.md)
  : Plot full fits

- [`plot_pars()`](https://lindeloev.github.io/mcp/dev/reference/plot_pars.md)
  : Plot individual parameters

- [`pp_check()`](https://lindeloev.github.io/mcp/dev/reference/pp_check.md)
  : Posterior Predictive Checks For Mcpfit Objects

- [`summary(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
  [`fixef(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
  [`ranef(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
  [`print(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/summary.mcpfit.md)
  : Summarise mcpfit objects

- [`predict(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  [`fitted(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  [`log_lik(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  [`residuals(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md)
  :

  Fitted and predicted values of `mcp` models fits

- [`loo(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
  [`waic(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/dev/reference/loo.mcpfit.md)
  : Information Criteria for Model Comparison

- [`add_loglik()`](https://lindeloev.github.io/mcp/dev/reference/add_loglik.md)
  : Add Log-Likelihood to an mcpfit Object.

- [`hypothesis()`](https://lindeloev.github.io/mcp/dev/reference/hypothesis.md)
  : Test Hypotheses Concerning Individual Parameters

- [`as_draws()`](https://lindeloev.github.io/mcp/dev/reference/as_draws.mcpfit.md)
  [`as_draws_df()`](https://lindeloev.github.io/mcp/dev/reference/as_draws.mcpfit.md)
  [`as_draws_array()`](https://lindeloev.github.io/mcp/dev/reference/as_draws.mcpfit.md)
  [`as_draws_matrix()`](https://lindeloev.github.io/mcp/dev/reference/as_draws.mcpfit.md)
  [`as_draws_rvars()`](https://lindeloev.github.io/mcp/dev/reference/as_draws.mcpfit.md)
  :

  Extract MCMC Draws from `mcpfit` Objects

- [`mcp-package`](https://lindeloev.github.io/mcp/dev/reference/mcp-package.md)
  : mcp: Regression with Multiple Change Points

## Axillary functions

These are used internally by mcp, but are exposed here since they may be
useful for other purposes. Most other useful internal functions deliver
the result already in `mcp(segments, sample = FALSE)`, so
[`mcp()`](https://lindeloev.github.io/mcp/dev/reference/mcp.md) will be
their API.

- [`sd_to_prec()`](https://lindeloev.github.io/mcp/dev/reference/sd_to_prec.md)
  : Transform a JAGS Prior from SD to Precision.

- [`logit()`](https://lindeloev.github.io/mcp/dev/reference/logit.md) :
  Logit function

- [`ilogit()`](https://lindeloev.github.io/mcp/dev/reference/ilogit.md)
  : Inverse logit function

- [`probit()`](https://lindeloev.github.io/mcp/dev/reference/probit.md)
  : Probit function

- [`phi()`](https://lindeloev.github.io/mcp/dev/reference/phi.md) :
  Inverse probit function

- [`is.mcpfit()`](https://lindeloev.github.io/mcp/dev/reference/is.mcpfit.md)
  :

  Checks if the Argument is an `mcpfit` Object

## Families

Distributional families that are not available in base R.

- [`bernoulli()`](https://lindeloev.github.io/mcp/dev/reference/bernoulli.md)
  : Bernoulli Family for mcp
- [`negbinomial()`](https://lindeloev.github.io/mcp/dev/reference/negbinomial.md)
  : Negative Binomial for mcp

## Help and demos

Some showcases of typical mcp analyses. Most of these are discussed on
the [front page](https://lindeloev.github.io/mcp).

- [`mcp_example()`](https://lindeloev.github.io/mcp/dev/reference/mcp_example.md)
  [`mcp_example_data()`](https://lindeloev.github.io/mcp/dev/reference/mcp_example.md)
  : Run example models
- [`demo_fit`](https://lindeloev.github.io/mcp/dev/reference/demo_fit.md)
  : Pre-fitted example mcp model

## Miscellaneous

Stuff you would not usually consult directly.

- [`mcpfit-class`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  :

  Class `mcpfit` of Models Fitted with the mcp Package

- [`ndraws()`](https://lindeloev.github.io/mcp/dev/reference/draws-index-mcp.md)
  [`nchains()`](https://lindeloev.github.io/mcp/dev/reference/draws-index-mcp.md)
  [`niterations()`](https://lindeloev.github.io/mcp/dev/reference/draws-index-mcp.md)
  :

  Index `mcpfit` objects

- [`interpolate_newdata()`](https://lindeloev.github.io/mcp/dev/reference/interpolate_newdata.md)
  : Returns a data.frame with all combos of predictors

- [`mcpfamily()`](https://lindeloev.github.io/mcp/dev/reference/mcpfamily.md)
  [`is.mcpfamily()`](https://lindeloev.github.io/mcp/dev/reference/mcpfamily.md)
  : Create or Test Objects of Class "mcpfamily"

- [`prior_summary()`](https://lindeloev.github.io/mcp/dev/reference/prior_summary.md)
  : Summarise priors used by an mcp model

- [`print(`*`<mcplist>`*`)`](https://lindeloev.github.io/mcp/dev/reference/print.mcplist.md)
  : Print mcplist

- [`print(`*`<mcptext>`*`)`](https://lindeloev.github.io/mcp/dev/reference/print.mcptext.md)
  : Nice Printing of Multiline Texts
