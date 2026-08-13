#' @keywords internal
#' @import patchwork
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @importFrom loo loo waic loo_compare
#' @importFrom stats ave binomial coef family formula gaussian mad median model.frame nobs
#' @importFrom utils capture.output
#' @importFrom rlang !! := .env
"_PACKAGE"

utils::globalVariables(c(".data", ".env"))
