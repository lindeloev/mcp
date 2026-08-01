# Package index

## Using mcp

Functions for everyday use of mcp.

- [`mcp()`](https://lindeloev.github.io/mcp/reference/mcp.md) : Fit
  Multiple Linear Segments And Their Change Points

- [`plot(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/plot.mcpfit.md)
  [`plot_dpar()`](https://lindeloev.github.io/mcp/reference/plot.mcpfit.md)
  : Plot full fits

- [`plot_pars()`](https://lindeloev.github.io/mcp/reference/plot_pars.md)
  : Plot individual parameters

- [`pp_check()`](https://lindeloev.github.io/mcp/reference/pp_check.md)
  : Posterior Predictive Checks For Mcpfit Objects

- [`summary(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md)
  [`fixef(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md)
  [`ranef(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md)
  [`print(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/summary.mcpfit.md)
  : Summarise mcpfit objects

- [`predict(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  [`fitted(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  [`log_lik(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  [`residuals(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/execute-mcp-model.md)
  :

  Fitted and predicted values of `mcp` models fits

- [`loo(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/loo.mcpfit.md)
  [`waic(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/loo.mcpfit.md)
  : Information Criteria for Model Comparison

- [`add_loglik()`](https://lindeloev.github.io/mcp/reference/add_loglik.md)
  : Add Log-Likelihood to an mcpfit Object.

- [`hypothesis()`](https://lindeloev.github.io/mcp/reference/hypothesis.md)
  : Test Hypotheses Concerning Individual Parameters

- [`as_draws(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/as_draws.mcpfit.md)
  :

  Extract MCMC Draws from `mcpfit` Objects

- [`mcp-package`](https://lindeloev.github.io/mcp/reference/mcp-package.md)
  : mcp: Regression with Multiple Change Points

## Axillary functions

These are used internally by mcp, but are exposed here since they may be
useful for other purposes. Most other useful internal functions deliver
the result already in `mcp(segments, sample = FALSE)`, so
[`mcp()`](https://lindeloev.github.io/mcp/reference/mcp.md) will be
their API.

- [`sd_to_prec()`](https://lindeloev.github.io/mcp/reference/sd_to_prec.md)
  : Transform a JAGS Prior from SD to Precision.

- [`logit()`](https://lindeloev.github.io/mcp/reference/logit.md) :
  Logit function

- [`ilogit()`](https://lindeloev.github.io/mcp/reference/ilogit.md) :
  Inverse logit function

- [`probit()`](https://lindeloev.github.io/mcp/reference/probit.md) :
  Probit function

- [`phi()`](https://lindeloev.github.io/mcp/reference/phi.md) : Inverse
  probit function

- [`is.mcpfit()`](https://lindeloev.github.io/mcp/reference/is.mcpfit.md)
  :

  Checks if the Argument is an `mcpfit` Object

## Families

Distributional families that are not available in base R.

- [`bernoulli()`](https://lindeloev.github.io/mcp/reference/bernoulli.md)
  : Bernoulli Family for mcp
- [`negbinomial()`](https://lindeloev.github.io/mcp/reference/negbinomial.md)
  : Negative Binomial for mcp

## Help and demos

Some showcases of typical mcp analyses. Most of these are discussed on
the [front page](https://lindeloev.github.io/mcp).

- [`mcp_example()`](https://lindeloev.github.io/mcp/reference/mcp_example.md)
  [`mcp_example_data()`](https://lindeloev.github.io/mcp/reference/mcp_example.md)
  : Get example models and data

- [`demo_fit`](https://lindeloev.github.io/mcp/reference/demo_fit.md) :

  Example `mcpfit`

## Miscellaneous

Stuff you would not usually consult directly.

- [`mcpfit-class`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  :

  Class `mcpfit` of Models Fitted with the mcp Package

- [`niterations(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/draws-index-mcp.md)
  [`nchains(`*`<mcpfit>`*`)`](https://lindeloev.github.io/mcp/reference/draws-index-mcp.md)
  :

  Index `mcpfit` objects

- [`interpolate_newdata()`](https://lindeloev.github.io/mcp/reference/interpolate_newdata.md)
  : Returns a data.frame with all combos of predictors

- [`mcpfamily()`](https://lindeloev.github.io/mcp/reference/mcpfamily.md)
  [`is.mcpfamily()`](https://lindeloev.github.io/mcp/reference/mcpfamily.md)
  : Create or Test Objects of Class "mcpfamily"

- [`prior_summary()`](https://lindeloev.github.io/mcp/reference/prior_summary.md)
  : Summarise priors used by an mcp model

- [`print(`*`<mcplist>`*`)`](https://lindeloev.github.io/mcp/reference/print.mcplist.md)
  : Print mcplist

- [`print(`*`<mcptext>`*`)`](https://lindeloev.github.io/mcp/reference/print.mcptext.md)
  : Nice Printing of Multiline Texts
