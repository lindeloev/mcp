#' Bernoulli family for mcp
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
#'
bernoulli = function(link = "logit") {
  out = list(
    family = "bernoulli",
    link = "logit"
  )
  class(out) = "family"
  return(out)
}
