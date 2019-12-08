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
#' @keywords internal
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
      link = "identity",
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
#' @param eta A vector of logits
ilogit = stats::binomial(link = "logit")$linkinv

#' Logit function
#'
#' @aliases logit
#' @param mu A vector of probabilities (0-1)
logit = stats::binomial(link = "logit")$linkfun



#' Converts logical(0) to null. Returns x otherwise
#'
#'@aliases logical0_to_null
#'@keywords internal
#'@param x Anything
#'@return NULL or x

logical0_to_null = function(x) {
  if (length(x) > 0)
    return(x)
  else return(NULL)
}

#' Extracts the order from ARMA parameter name(s)
#'
#' If several names are provided (vector), it returns the maximum. If `pars_arma`
#' is an empty string, it returns `0`.
#'
#' @aliases get_arma_order
#' @keywords internal
#' @param pars_arma Character vector
#' @return integer
get_arma_order = function(pars_arma) {
  if (length(pars_arma) > 0) {
    order_str = sub("(ma|ar)([0-9]+).*", "\\2", pars_arma)
    order_max = max(as.numeric(order_str))
    return(order_max)
  } else {
    return(0)
  }
}

#' Throws an error if a number/vector contains non-numeric, decimal, or less-than-lower
#'
#' The expected behavior of is.integer, with informative error messages.
#'
#' @aliases check_integer
#' @keywords internal
#' @param x Numeric value or vector
#' @param name Name to show in error message.
#' @param lower the smallest allowed value. lower = 1 checks for positive integers.
#'
check_integer = function(x, name, lower = -Inf) {
  greater_than = ifelse(lower == -Inf, " ", paste0(" >= ", lower, " "))
  if (!is.numeric(x))
    stop("Only integers", greater_than, "allowed for '", name, "'")

  if (!all(x == floor(x)) | !all(x >= lower))
    stop("Only integers", greater_than, "allowed for '", name, "'")

  TRUE
}



# Ask reminder questions for CRAN export
release_questions = function() {
  c(
    "Have you run the test of fits?",
    "Have you built the README plots and checked them? source('man/figures/make_README_plots.R')",
    "Have you re-built the site using pkgdown::build_site() AFTER deleting caches of articles?",
    "Have you checked all articles and plots after re-building the site?",
    "Have you run the script to insert the correct logo.png in the HTML meta?"
  )
}



