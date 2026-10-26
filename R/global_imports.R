# Importing from coda loads it, so `[` dispatches to its mcmc.list method even when
# JAGS ran in future workers and coda was never loaded in this session.
#' @keywords internal
#' @import patchwork
#' @importFrom dplyr .data
#' @importFrom bayesplot pp_check
#' @importFrom coda varnames
#' @importFrom nlme fixef ranef
#' @importFrom loo loo waic loo_compare
#' @importFrom stats ave binomial coef confint family formula gaussian mad median model.frame nobs vcov
#' @importFrom utils capture.output
#' @importFrom rlang !! := .env
NULL

utils::globalVariables(c(".data", ".env"))
