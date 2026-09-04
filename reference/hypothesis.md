# Test Hypotheses Concerning Individual Parameters

Returns posterior probabilities and Bayes Factors for flexible
hypotheses involving model parameters. The documentation for the
argument `hypotheses` below shows examples of how to specify hypotheses,
and [read worked examples on the mcp
website](https://lindeloev.github.io/mcp/articles/comparison.html). For
directional hypotheses, `hypothesis` executes the hypothesis string in a
data-frame environment and summarises the proportion of posterior and
prior draws where the expression evaluates to TRUE. The Bayes factor is
the posterior odds divided by the prior odds. For equality hypotheses, a
Savage-Dickey ratio is computed. Both kinds of Bayes factor require
prior draws, so remember `mcp(..., sample = "both")`. This function is
heavily inspired by the `hypothesis` function from the `brms` package.
When `prior = TRUE`, the summary is based on prior draws and `BF` is
`NA`.

## Usage

``` r
hypothesis(fit, ...)

# S3 method for class 'mcpfit'
hypothesis(fit, hypotheses, width = 0.95, prior = FALSE, ...)
```

## Arguments

- fit:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object.

- ...:

  Must be empty. Reserved for future use.

- hypotheses:

  String representation of a logical test involving model parameters.
  Takes R code that evaluates to TRUE or FALSE in a vectorized way.

  **Directional hypotheses** are specified using \<, \>, \<=, or \>=.
  `hypothesis` returns the posterior probability \\P(H \mid
  \text{data})\\ and the Bayes factor in favor of the stated hypothesis
  \\H\\: \$\$\text{BF}\_{10} = \frac{P(H \mid \text{data}) / (1 - P(H
  \mid \text{data}))}{P(H) / (1 - P(H))}\$\$ where \\P(H)\\ is the prior
  probability of \\H\\. The Bayes factor requires both prior and
  posterior draws from `mcp(sample = "both")`. For example:

  - `"cp_1 > 30"`: the first change point is above 30.

  - `"Intercept_1 > Intercept_2"`: the intercept is greater in segment 1
    than 2.

  - `"x_2 - x_1 <= 3"`: the difference between slope 1 and 2 is less
    than or equal to 3.

  - `"Intercept_1 > -2 & Intercept_1 < 2"`: Intercept_1 is between -2
    and 2 (an interval hypothesis). This can be useful as a Region Of
    Practical Equivalence test (ROPE).

  - `"cp_1^2 < 30 | (log(x_1) + log(x_2)) > 5"`: be creative.

  - `` "`cp_1_id[1]` > `cp_1_id[2]`" ``: id1 is greater than id2, as
    estimated through the group-level change-point deviation for `id` at
    change point 1 (starting segment 2). Note that ``` `` ``` are
    required when using `[i]`.

  **Equality hypotheses** use the equal sign (=) and a Savage-Dickey
  density ratio: posterior density divided by prior density at the
  tested point equality

  \\\theta = \theta_0\\: \$\$\text{BF}\_{01} = \frac{p(\theta = \theta_0
  \mid \text{data})}{p(\theta = \theta_0)}\$\$

  where \\\theta\\ is the evaluated parameter (or affine contrast),
  \\\theta_0\\ is the hypothesized null value, \\p(\theta = \theta_0
  \mid \text{data})\\ is the posterior density, and \\p(\theta =
  \theta_0)\\ is the prior density. This is a Bayes factor for a nested
  point-null model against the fitted continuous model
  (\\\text{BF}\_{10} = 1 / \text{BF}\_{01}\\). Prior and posterior draws
  are required, using `mcp(sample = "both")`.

  The point-null model's nuisance prior is the fitted model's
  conditional prior at the equality. Equality tests are limited to named
  scalar parameters and affine contrasts.

  Examples:

  - `"cp_1 = 30"`: compare the point-null model `cp_1 = 30` with the
    continuous alternative.

  - `"Intercept_1 - Intercept_2 = 0"`: compare equal segment intercepts
    with the continuous alternative.

- width:

  Float. The width of the central posterior interval (between 0 and 1).

- prior:

  Logical. Summarise prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)?

## Value

A data.frame with a row per hypothesis and the following columns:

- `hypothesis` is the hypothesis; often re-arranged to test against
  zero.

- `mean` is the posterior mean of the left-hand side of the hypothesis,
  or the prior mean when `prior = TRUE`.

- `lower` is the lower bound of the central posterior interval of width
  `width`, or the corresponding prior interval when `prior = TRUE`.

- `upper` is the upper bound of ditto.

- `prob` is the posterior probability of a directional hypothesis, or
  the prior probability when `prior = TRUE`. It is `NA` for equality
  hypotheses, which compare models rather than an event within the
  fitted model.

- `BF` Bayes Factor in favor of the hypothesis. For "=" it is the
  Savage-Dickey density ratio. For directional hypotheses, it is the
  posterior odds divided by the prior odds. It is `NA` when
  `prior = TRUE` or when prior draws are not available (e.g.,
  `sample = "post"`).

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
# demo_fit contains both posterior and prior draws
# A directional hypothesis returns its posterior probability and Bayes factor
hypothesis(demo_fit, "cp_1 > 30")
#>      hypothesis     mean     lower    upper  prob       BF
#> 1 cp_1 - 30 > 0 1.459879 -1.938938 4.998374 0.813 4.691063

# Combine directional statements for an interval (a ROPE-style hypothesis)
hypothesis(demo_fit, "cp_1 > 20 & cp_1 < 30")
#>              hypothesis mean lower upper  prob       BF
#> 1 cp_1 > 20 & cp_1 < 30   NA    NA    NA 0.187 1.189817

# Evaluate several directional hypotheses at once
hypothesis(demo_fit, c("cp_1 > 20", "cp_2 > 70"))
#>      hypothesis     mean      lower     upper  prob       BF
#> 1 cp_1 - 20 > 0 11.45988  8.0610620 14.998374 1.000      Inf
#> 2 cp_2 - 70 > 0  1.11715 -0.5594432  2.765253 0.818 4.440893

# Directional hypotheses can be used for a focused posterior question
hypothesis(demo_fit, "cp_1 > 30")
#>      hypothesis     mean     lower    upper  prob       BF
#> 1 cp_1 - 30 > 0 1.459879 -1.938938 4.998374 0.813 4.691063

# Inspect the corresponding prior probability without a Bayes factor
hypothesis(demo_fit, "cp_1 > 30", prior = TRUE)
#>      hypothesis     mean     lower    upper  prob BF
#> 1 cp_1 - 30 > 0 3.827508 -27.98519 55.25295 0.481 NA
```
