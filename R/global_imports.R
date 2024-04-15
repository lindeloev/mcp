#' @keywords internal
#' @import patchwork
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @importFrom loo loo waic loo_compare
#' @importFrom stats gaussian binomial
#' @importFrom rlang !! :=
"_PACKAGE"


# Hack to make R CMD pass for function geom_cp_density()
utils::globalVariables(c("value", ".chain", "cp_name", "."))
