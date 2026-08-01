# Test Hypotheses Concerning Individual Parameters

Returns posterior probabilities and Bayes Factors for flexible
hypotheses involving model parameters. The documentation for the
argument `hypotheses` below shows examples of how to specify hypotheses,
and [read worked examples on the mcp
website](https://lindeloev.github.io/mcp/articles/comparison.html). For
directional hypotheses, `hypothesis` executes the hypothesis string in a
data-frame environment and summarises the proportion of posterior and
prior samples where the expression evaluates to TRUE. The Bayes factor
is the posterior odds divided by the prior odds. For equals-hypotheses,
a Savage-Dickey ratio is computed. Both kinds of Bayes factor require
prior samples, so remember `mcp(..., sample = "both")`. This function is
heavily inspired by the `hypothesis` function from the `brms` package.

## Usage

``` r
hypothesis(fit, hypotheses, width = 0.95, digits = 3, prior = FALSE)
```

## Arguments

- fit:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object.

- hypotheses:

  String representation of a logical test involving model parameters.
  Takes R code that evaluates to TRUE or FALSE in a vectorized way.

  Directional hypotheses are specified using \<, \>, \<=, or \>=.
  `hypothesis` returns the posterior probability and the Bayes factor in
  favor of the stated hypothesis. The Bayes factor requires both prior
  and posterior samples from `mcp(sample = "both")`. For example:

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
    estimated through the varying-by-"id" change point in segment 1.
    Note that ``` `` ``` required for varying effects.

  Hypotheses can also test equality using the equal sign (=). This runs
  a Savage-Dickey test, i.e., the proportion by which the probability
  density has increased from the prior to the posterior at a given
  value. Therefore, it requires `mcp(sample = "both")`. There are two
  requirements: First, there can only be one equal sign, so don't use
  and (&) or or (\|). Second, the point to test has to be on the right,
  and the variables on the left.

  - `"cp_1 = 30"`: is the first change point at 30? Or to be more
    precise: by what factor has the credence in cp_1 = 30 risen/fallen
    when conditioning on the data, relative to the prior credence?

  - `"Intercept_1 + Intercept_2 = 0"`: Is the sum of two intercepts
    zero?

  - `` "`cp_1_id[John]`/`cp_1_id[Erin]` = 2" ``: is the varying change
    point for John (which is relative to \`cp_1“) double that of Erin?

- width:

  Float. The width of the central posterior interval (between 0 and 1).

- digits:

  a non-null value for digits specifies the minimum number of
  significant digits to be printed in values. The default, NULL, uses
  getOption("digits"). (For the interpretation for complex numbers see
  signif.) Non-integer values will be rounded down, and only values
  greater than or equal to 1 and no greater than 22 are accepted.

- prior:

  TRUE/FALSE. Summarise prior instead of posterior?

## Value

A data.frame with a row per hypothesis and the following columns:

- `hypothesis` is the hypothesis; often re-arranged to test against
  zero.

- `mean` is the posterior mean of the left-hand side of the hypothesis.

- `lower` is the lower bound of the central posterior interval of width
  `width`.

- `upper` is the upper bound of ditto.

- `p` Posterior probability. For "=" (Savage-Dickey), it is the BF
  converted to p. For directional hypotheses, it is the proportion of
  samples that returns TRUE.

- `BF` Bayes Factor in favor of the hypothesis. For "=" it is the
  Savage-Dickey density ratio. For directional hypotheses, it is the
  posterior odds divided by the prior odds.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
