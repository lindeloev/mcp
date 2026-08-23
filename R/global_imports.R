#' @keywords internal
#' @import patchwork
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @importFrom loo loo waic loo_compare
#' @importFrom stats ave binomial coef confint family formula gaussian mad median model.frame nobs vcov
#' @importFrom utils capture.output
#' @importFrom rlang !! := .env
NULL

utils::globalVariables(c(".data", ".env"))
