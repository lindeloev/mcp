#' Bernoulli family for mcp
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
#'
bernoulli = function(link = "logit") {
  assert_value(link, allowed = c("identity", "logit", "probit"))

  # Just copy binomial()
  family = binomial(link = link)
  family$family = "bernoulli"
  family$likfun = function(x, prob, log = FALSE) dbinom(x, 1, prob, log)
  mcpfamily(family)
}

#' Exponential family for mcp
#'
#' @aliases exponential
#' @param link Link function (Character).
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


#' Negative binomial for mcp
#'
#' Parameterized as `mu` (mean; poisson lambda) and `size` (a shape parameter),
#' so you can do `rnbinom(10, mu = 10, size = 1)`. Read more in the doc for `rnbinom`,
#'
#' @aliases negbinomial
#' @param link Link function (Character).
#' @export
negbinomial = function(link = "log") {
  assert_value(link, allowed = c("log", "identity"))

  family = list(
    family = "negbinomial",
    link = link  # on lambda
  )
  class(family) = "family"
  family = mcpfamily(family)
}


#' Create or test objects of type class "mcpfamily"
#'
#' @aliases mcpfamily
#' @param x A family object, e.g., `binomial(link = "identity")`.
#' @seealso \code{\link{family}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
mcpfamily = function(x) {
  family = x
  tmp = x  # Hack: the name "family" causes errors in subset()
  assert_types(family, "family")

  # Get default priors for RHS
  dpar_prior = subset(default_dpar_priors, family == tmp$family & link == tmp$link)
  if (nrow(dpar_prior) == 0)
    stop("mcp has no default priors for ", family$family, "(link = \"", family$link, "\") so it's likely not supported. See `mcpfamily()` on how to create a custom family.")

  # Add priors and list parameters
  dpar_prior = rbind(dpar_prior, subset(default_common_priors, dpar == "ar"))
  family$default_prior = dpar_prior
  family$dpars = unique(dpar_prior$dpar)

  # Add likelihood function for stats families
  if (is.null(family$likfun)) {
    loglik_funs = list(
      gaussian = dnorm,
      binomial = dbinom,
      poisson = dpois
    )
    family$likfun = loglik_funs[[family$family]]
  }

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

  # Add link functions as R functions, if they are not already present
  if (rlang::has_name(family, "linkfun") == FALSE)
    family$linkfun = get(family$link)
  if (rlang::has_name(family, "linkinv") == FALSE)
    family$linkinv = get(family$linkinv_str)

  class(family) = c("mcpfamily", "family")
  family
}


#' @aliases is.mcpfamily
#' @describeIn mcpfamily
#' @export
is.mcpfamily = function(x) {
  if (inherits(x, "mcpfamily") == FALSE)
    return(FALSE)

  assert_types(x, "family")
  assert_types(x$default_prior, "tibble", "data.frame")
  assert_types(x$dpars, "character")
  assert_types(x$linkfun_str, "character", len = 1)
  assert_types(x$linkinv_str, "character", len = 1)
  assert_types(x$linkfun, "function")
  assert_types(x$linkinv, "function")

  TRUE
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
