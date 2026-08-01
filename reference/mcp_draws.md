# Get tidy samples with or without varying effects

Extract posterior or prior draws formatted as tidy data frames

## Usage

``` r
mcp_draws(
  fit,
  population = TRUE,
  varying = TRUE,
  absolute = FALSE,
  prior = FALSE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated()
)
```

## Arguments

- fit:

  An
  [`mcpfit`](https://lindeloev.github.io/mcp/reference/mcpfit-class.md)
  object

- population:

  - `TRUE` All population effects. Same as `fit$pars$population`.

  - `FALSE` No population effects. Same as
    [`c()`](https://rdrr.io/r/base/c.html).

  - Character vector: Only include specified population parameters - see
    `fit$pars$population`.

- varying:

  One of:

  - `TRUE` All varying effects (`fit$pars$varying`).

  - `FALSE` No varying effects ([`c()`](https://rdrr.io/r/base/c.html)).

  - `"cp"` or `"predictor"`: All varying effects belonging to that part
    of the model.

  - Character vector: Only include specified varying parameters - see
    `fit$pars$varying`.

- absolute:

  - `TRUE` Returns the absolute location of all varying change points.

  - `FALSE` Just returns the varying effects.

  - Character vector: Only do absolute transform for these varying
    parameters - see `fit$pars$varying`.

- prior:

  TRUE/FALSE. Summarise prior instead of posterior?

- ndraws:

  Integer or `NULL`. Number of posterior draws to return/summarise. If
  there are varying effects, this is the number of draws from each
  varying group. `NULL` means "all". Ignored if both are `FALSE`. More
  samples trade speed for accuracy.

- nsamples:

  Deprecated. Use `ndraws` instead.

## Value

`tibble` of posterior draws in `tidybayes` format.

## Details

Returns in a format useful for `fit$simulate()` with population
parameters in wide format and varying effects in long format (the number
of rows will be `ndraws * n_levels_in_varying`).

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
