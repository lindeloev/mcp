# Transform an mcp prior to the parameterization used by JAGS.

**\[deprecated\]**

## Usage

``` r
sd_to_prec(prior_str)
```

## Arguments

- prior_str:

  String. A JAGS prior. Can be truncated, e.g.
  `dt(3, 2, 1) T(my_var, )`.

## Value

A string

## Details

This function is deprecated. JAGS uses precision rather than SD for some
distributions. For example, this function converts `dnorm(4.2, 1.3)`
into `dnorm(4.2, 1/1.3^2)`. It allows users to specify priors using
conventional scale parameters before they are translated to JAGS code.
Users normally do not need to call this function.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
