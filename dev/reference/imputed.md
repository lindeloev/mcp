# Imputed Responses

Get posterior draws of the missing responses in `fit$data`. Use
`complete = TRUE` to include the observed responses in each draw. For
predictions of fresh responses, including at the missing rows, use
[`predict()`](https://lindeloev.github.io/mcp/dev/reference/execute-mcp-model.md).
See the [missing data
article](https://lindeloev.github.io/mcp/articles/missing.html) for
plotting imputations and completing your original data.

## Usage

``` r
imputed(
  fit,
  complete = FALSE,
  summary = TRUE,
  width = 0.95,
  rate = FALSE,
  ndraws = NULL,
  draws_format = "tidy",
  ...
)
```

## Arguments

- fit:

  An `mcpfit` object with posterior draws.

- complete:

  Logical. Include observed responses (`TRUE`) or only missing responses
  (`FALSE`)?

- summary:

  Logical. Summarise the draws for each returned row (`TRUE`) or return
  individual draws (`FALSE`)?

- width:

  Width of the central interval (strictly between 0 and 1), returned in
  the `.lower` and `.upper` columns. `FALSE` omits the interval. Used
  when `summary = TRUE`.

- rate:

  Logical. For binomial models, return counts (`FALSE`, the default) or
  proportions (`TRUE`).

- ndraws:

  Integer or `NULL`. Number of draws to select at random. `NULL` uses
  all draws.

- draws_format:

  One of "tidy" or "matrix". Controls the output format when
  `summary == FALSE`.

- ...:

  Must be empty. Reserved for future use.

## Value

With `summary = TRUE`, a tibble with `.mean`, `.sd`, and, by default,
`.lower` and `.upper`. With `summary = FALSE`, a tidy tibble with
`.imputed`, `.chain`, `.iteration`, and `.draw`, or a
draws-by-observation matrix when `draws_format = "matrix"`.

Tibbles include the original columns of `fit$data` and `.row`, the row
number in `fit$data`. Matrix columns are named by `.row`. `.draw`
matches the draw in
[`as_draws_df()`](https://lindeloev.github.io/mcp/dev/reference/as_draws.mcpfit.md).

## Examples

``` r
# \donttest{
fit = mcp_example("missing", plot = FALSE)
#> NA values detected in 'y'. JAGS will treat them as latent responses and impute them during sampling.
imputed(fit, width = 0.8)
#> # A tibble: 10 × 8
#>        y     x state  .row .mean   .sd .lower .upper
#>    <dbl> <int> <fct> <int> <dbl> <dbl>  <dbl>  <dbl>
#>  1    NA     8 B         8  35.4  4.45  29.7    41.2
#>  2    NA    19 A        19  15.5  4.30   9.99   21.0
#>  3    NA    27 A        27  17.2  4.34  11.6    22.7
#>  4    NA    28 B        28  39.4  4.30  33.9    44.7
#>  5    NA    29 A        29  17.6  4.40  11.9    23.3
#>  6    NA    30 B        30  39.8  4.39  34.2    45.3
#>  7    NA    31 A        31  17.9  4.43  12.3    23.6
#>  8    NA    68 B        68  41.2  4.45  35.6    46.9
#>  9    NA    84 B        84  34.6  4.37  29.0    40.2
#> 10    NA    96 B        96  29.6  4.44  23.9    35.3
head(imputed(fit, summary = FALSE, ndraws = 3))
#> # A tibble: 6 × 8
#>       y     x state  .row .chain .iteration .draw .imputed
#>   <dbl> <int> <fct> <int>  <int>      <int> <int>    <dbl>
#> 1    NA     8 B         8      1       1122  1122     32.7
#> 2    NA     8 B         8      3        681  6681     37.4
#> 3    NA     8 B         8      3       2856  8856     32.9
#> 4    NA    19 A        19      1       1122  1122     15.4
#> 5    NA    19 A        19      3        681  6681     16.8
#> 6    NA    19 A        19      3       2856  8856     15.6
# }
```
