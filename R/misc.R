#' Bernoulli family for mcp
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
#'

bernoulli = function(link = "logit") {
  get_family("bernoulli", link)
}


#' A family object to store link functions between R and JAGS.
#'
#' This will make more sense once more link functions / families are added.
#'
#' @aliases get_family
#' @param name Name of the family
#' @param link Link function name. Accepts "logit", "log", and "identity"
get_family = function(name, link) {
  # Exciting!
  if (link == "logit") {
    out = list(
      family = name,
      link = "logit",
      link_jags = "logit",
      link_r = "logit",  # included in mcp
      linkinv_jags = "ilogit",
      linkinv_r = "ilogit"  # included in mcp
    )
  }

  # Pretty boring
  if (link == "log") {
    out = list(
      family = name,
      link = "log",
      link_jags = "log",
      link_r = "log",
      linkinv_jags = "exp",
      linkinv_r = "exp"
    )
  }

  # Pretty boring
  if (link == "identity") {
    out = list(
      family = name,
      link = "",
      link_jags = "",
      link_r = "",
      linkinv_jags = "",
      linkinv_r = ""
    )
  }

  class(out) = "family"
  return(out)
}





#' Inverse logit function
#'
#' @aliases ilogit
#' @param x A vector of logits
ilogit = stats::binomial(link = "logit")$linkinv

#' Logit function
#'
#' @aliases logit
#' @param x A vector of probabilities (0-1)
logit = stats::binomial(link = "logit")$linkfun
