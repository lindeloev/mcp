#' @keywords internal
#' @import patchwork
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @importFrom loo loo waic loo_compare
#' @importFrom stats gaussian binomial mad median
#' @importFrom rlang !! := .env
"_PACKAGE"

utils::globalVariables(c(".data", ".env"))
