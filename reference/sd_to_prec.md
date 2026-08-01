# Transform a JAGS Prior from SD to Precision.

JAGS uses precision rather than SD. This function converts
`dnorm(4.2, 1.3)` into `dnorm(4.2, 1/1.3^2)`. It allows users to specify
priors using SD and then it's transformed for the JAGS code. It works
for the following distributions:
dnorm\|dt\|dcauchy\|ddexp\|dlogis\|dlnorm. In all of these, tau/sd is
the second parameter.

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

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
