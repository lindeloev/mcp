# Get tidy draws with or without group-level effects

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
  [`mcpfit`](https://lindeloev.github.io/mcp/dev/reference/mcpfit-class.md)
  object

- population:

  - `TRUE` All population-level model parameters.

    - `FALSE` No population-level effects. Same as
      [`c()`](https://rdrr.io/r/base/c.html).

    - Character vector: Only include specified population-level
      parameters.

- varying:

  Group-level effects. One of:

  - `TRUE` All group-level deviations.

  - `FALSE` No group-level deviations
    ([`c()`](https://rdrr.io/r/base/c.html)).

  - `"cp"` or `"predictor"`: All group-level deviations belonging to
    that part of the model.

  - Character vector: Only include specified group-level parameters.

- absolute:

  - `TRUE` Returns the absolute location of all group-specific change
    points.

    - `FALSE` Return the group-level deviations.

    - Character vector: Apply the absolute transform only to these
      group-level parameters.

- prior:

  Logical. Summarise prior draws (`TRUE`) instead of posterior draws
  (`FALSE`, default)?

- ndraws:

  Integer or `NULL`. Number of posterior draws to return/summarise. If
  there are group-level effects, this is the number of draws from each
  group. `NULL` means "all". More draws trade speed for accuracy.

- nsamples:

  Deprecated. Use `ndraws` instead.

## Value

`tibble` of posterior draws in `tidybayes` format.

## Details

Returns in a format useful for `fit$simulate()` with population-level
parameters in wide format and group-level deviations in long format (the
number of rows is multiplied by the number of selected group levels).

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>
