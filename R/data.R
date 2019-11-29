#' A change point in a time series
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_ar.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 300 rows and 2 variables:
#' \describe{
#'   \item{time}{The x-axis.}
#'   \item{price}{The y-axis.}
#' }
"ex_ar"


#' Two change points between three binomial segments
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_binomial.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{x}{The x-axis.}
#'   \item{y}{The y-axis.}
#'   \item{N}{The number of trials for for `y`.}
#' }
"ex_binomial"


#' Two change points between three linear segments
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_demo.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{time}{The x-axis.}
#'   \item{response}{The y-axis.}
#' }
"ex_demo"


#' A change point between two plateaus
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_plateaus.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{x}{The x-axis.}
#'   \item{y}{The y-axis.}
#' }
"ex_plateaus"


#' A change point from plateau to quadratic
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_quadratic.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 81 rows and 2 variables:
#' \describe{
#'   \item{x}{The x-axis.}
#'   \item{y}{The y-axis.}
#' }
"ex_quadratic"


#' Two change points between three linear segments
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_rel_prior.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{x}{The x-axis.}
#'   \item{y}{The y-axis.}
#' }
"ex_rel_prior"


#' A change point between two trigonometric segments
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_trig.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 234 rows and 2 variables:
#' \describe{
#'   \item{x}{The x-axis.}
#'   \item{y}{The y-axis.}
#' }
"ex_trig"


#' Two change points between three heteroskedastic segments
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_variance.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{x}{The x-axis.}
#'   \item{y}{The y-axis.}
#' }
"ex_variance"


#' One change point varying by participant
#'
#' See how this data was simulated
#' [here](https://github.com/lindeloev/mcp/blob/master/data-raw/ex_varying.R)
#' using `fit$simulate()`, including which parameters were used. See an analysis
#' [here](https://lindeloev.github.io/mcp/).
#'
#' @format A data frame with 150 rows and 4 variables:
#' \describe{
#'   \item{x}{The x-axis.}
#'   \item{y}{The y-axis.}
#'   \item{id}{The participant id (character).}
#'   \item{id_numeric}{The participant id (numeric).}
#' }
"ex_varying"
