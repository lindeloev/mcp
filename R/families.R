#' Bernoulli family for mcp
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
#'
bernoulli = function(link = "logit") {
  assert_value(link, allowed = c("identity", "logit", "probit
                                 "))
  # Just copy binomial()
  family = binomial(link = link)
  family$family = "bernoulli"
  mcpfamily(family)
}

#' Exponential family for mcp
#'
#' @aliases exponential
#' @param link Link function
#' @export
#'
exponential = function(link = "identity") {
  assert_value(link, allowed = c("identity"))

  family = list(
    family = "exponential",
    link = link  # on lambda
  )
  class(family) = "family"
  family = mcpfamily(family)
}


#' Add A family object to store link functions between R and JAGS.
#'
#' This will make more sense once more link functions / families are added.
#'
#' @aliases mcpfamily
#' @keywords internal
#' @param family A family object, e.g., `binomial(link = "identity")`.
mcpfamily = function(family) {
  # Set linkfun_str
  if (family$link == "identity") {
    family$linkfun_str = ""
  } else {
    family$linkfun_str = family$link
  }

  # Set linkinv_str
  family$linkinv_str = switch(
    family$link,
    logit = "ilogit",
    probit = "phi",
    log = "exp",
    identity = ""
  )

  if (rlang::has_name(family, "linkfun") == FALSE)
    family$linkfun = eval(parse(text = family$link))
  if (rlang::has_name(family, "linkinv") == FALSE)
    family$linkinv = eval(parse(text = family$linkinv_str))

  return(family)
}



#' Logit function
#'
#' @aliases logit
#' @param mu A vector of probabilities (0.0 to 1.0)
#' @return A vector with same length as `mu`
#' @export
logit = stats::binomial(link = "logit")$linkfun

#' Inverse logit function
#'
#' @aliases ilogit
#' @param eta A vector of logits
#' @return A vector with same length as `eta`
#' @export
ilogit = stats::binomial(link = "logit")$linkinv


#' Probit function
#'
#' @aliases probit
#' @param mu A vector of probabilities (0.0 to 1.0)
#' @return A vector with same length as `mu`
#' @export
probit = stats::binomial(link = "probit")$linkfun


#' Inverse probit function
#'
#' @aliases phi
#' @param eta A vector of probits
#' @return A vector with same length as `mu`
#' @export
phi = stats::binomial(link = "probit")$linkinv
