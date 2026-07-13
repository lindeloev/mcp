#' Bernoulli Family for mcp
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
bernoulli = function(link = "logit") {
  assert_value(link, allowed = c("identity", "logit", "probit"))

  # Just copy binomial()
  family = binomial(link = link)
  family$family = "bernoulli"
  family$likfun = function(x, prob, log = FALSE) stats::dbinom(x, 1, prob, log)
  mcpfamily(family)
}

#' Negative Binomial for mcp
#'
#' Parameterized as `mu` (the conditional mean) and `shape` (the same quantity
#' as `size` in `rnbinom()`). Thus `Var(y) = mu + mu^2 / shape`, which approaches
#' the Poisson variance as `shape` approaches infinity.
#'
#' @aliases negbinomial
#' @param link Link function for `mu`.
#' @param link_shape Link function for `shape`.
#' @details `shape(1)` is added implicitly and is constant across segments unless
#'   a `shape()` formula is supplied. For example, `y ~ 1 + x + shape(1 + x)`
#'   models both the mean and shape. Regression coefficients for both dpars are
#'   on their link scales.
#' @export
negbinomial = function(link = "log", link_shape = "log") {
  assert_value(link, allowed = "log")
  assert_value(link_shape, allowed = "log")

  family = list(
    family = "negbinomial",
    link = link,
    link_shape = link_shape,
    likfun = function(x, mu, shape, log = FALSE) {
      stats::dnbinom(x, mu = mu, size = shape, log = log)
    }
  )
  class(family) = "family"
  family = mcpfamily(family)
}


#' Create or Test Objects of Class "mcpfamily"
#'
#' @aliases mcpfamily
#' @param x A family object, e.g., `binomial(link = "identity")`.
#' @seealso \code{\link{family}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
mcpfamily = function(x) {
  family = x
  tmp = x  # Hack: the name "family" causes errors in filter()
  assert_types(family, "family")

  # Get default priors for RHS
  dpar_prior = default_dpar_priors %>% dplyr::filter(.data$family == tmp$family & .data$link == tmp$link)
  if (nrow(dpar_prior) == 0)
    stop("mcp has no default priors for ", family$family, "(link = \"", family$link, "\") so it's likely not supported. See `mcpfamily()` on how to create a custom family.")

  # Add priors
  dpar_prior = rbind(dpar_prior, default_common_priors %>% dplyr::filter(.data$dpar == "ar"))
  family$default_prior = dpar_prior

  # Distributional parameters are properties of the family, not of the prior
  # table. Keeping this metadata separate prevents a missing/default prior from
  # silently removing a parameter from the likelihood.
  family = add_dpar_specs(family)

  # Add likelihood function for stats families
  if (is.null(family$likfun)) {
    loglik_funs = list(
      gaussian = stats::dnorm,
      binomial = stats::dbinom,
      poisson = stats::dpois
    )
    family$likfun = loglik_funs[[family$family]]
  }

  # Legacy mu-link fields remain available while all dpar links are also stored
  # in `family$dpar_specs` and `family$links`.
  family$linkfun_str = get_link_str(family$link)
  family$linkinv_str = get_link_str(family$link, inverse = TRUE)

  # Add link functions as R functions, if they are not already present
  if (rlang::has_name(family, "linkfun") == FALSE)
    family$linkfun = get_link_function(family$link)
  if (rlang::has_name(family, "linkinv") == FALSE)
    family$linkinv = get_link_function(family$link, inverse = TRUE)

  class(family) = c("mcpfamily", "family")
  family
}


#' Describe a distributional parameter
#'
#' This metadata is deliberately small. Family-specific likelihood code remains
#' explicit in `get_jags_code()` and `simulate_vectorized()`.
#'
#' @keywords internal
#' @noRd
new_dpar_spec = function(dpar, link, implicit = FALSE, lower = NA_real_) {
  tibble::tibble(
    dpar = dpar,
    link = link,
    implicit = implicit,
    lower = lower
  )
}


#' Default distributional parameters for built-in families
#'
#' Additional families can provide their own `dpar_specs` before calling
#' `mcpfamily()`.
#'
#' @keywords internal
#' @noRd
get_default_dpar_specs = function(family) {
  if (family$family == "gaussian") {
    specs = dplyr::bind_rows(
      new_dpar_spec("mu", family$link),
      new_dpar_spec("sigma", "identity", implicit = TRUE, lower = 1e-9)
    )
  } else if (family$family == "negbinomial") {
    specs = dplyr::bind_rows(
      new_dpar_spec("mu", family$link),
      new_dpar_spec("shape", family$link_shape, implicit = TRUE)
    )
  } else {
    # Poisson, Bernoulli, binomial, and custom mean-only families.
    specs = new_dpar_spec("mu", family$link)
  }

  specs
}


#' Names reserved for distributional-parameter formula wrappers
#'
#' Availability is still determined exclusively by `family$dpar_specs`. This
#' list only lets the formula parser distinguish an unsupported dpar wrapper
#' from an ordinary user-defined transformation function.
#'
#' @keywords internal
#' @noRd
known_dpar_wrappers = function() {
  c("sigma", "shape")
}


#' Add distributional-parameter metadata when absent
#'
#' The fallback keeps mcpfit objects saved by older mcp versions usable without
#' mutating their serialized family object.
#'
#' @keywords internal
#' @noRd
add_dpar_specs = function(family) {
  if (is.null(family$dpar_specs))
    family$dpar_specs = get_default_dpar_specs(family)
  assert_dpar_specs(family$dpar_specs)
  family$dpars = c(family$dpar_specs$dpar, "ar")
  family$links = stats::setNames(family$dpar_specs$link, family$dpar_specs$dpar)
  family
}


#' Validate distributional-parameter metadata
#'
#' @keywords internal
#' @noRd
assert_dpar_specs = function(x) {
  assert_types(x, "data.frame", "tibble")
  required = c("dpar", "link", "implicit", "lower")
  assert_data_cols(x, required)

  if (!is.character(x$dpar) || anyNA(x$dpar) || any(x$dpar == ""))
    stop("`family$dpar_specs$dpar` must contain non-empty parameter names.")
  if (!is.character(x$link) || anyNA(x$link))
    stop("`family$dpar_specs$link` must contain link names without missing values.")
  if (anyDuplicated(x$dpar))
    stop("Each distributional parameter must occur exactly once in `family$dpar_specs`.")
  if ("mu" %notin% x$dpar)
    stop("`family$dpar_specs` must contain the response-mean parameter 'mu'.")
  if (any(x$dpar %in% c("epred", "ar")))
    stop("'epred' and 'ar' are reserved and cannot be family distributional parameters.")

  supported_links = c("identity", "log", "logit", "probit")
  if (any(x$link %notin% supported_links))
    stop("Unsupported dpar link(s): ", and_collapse(unique(x$link[x$link %notin% supported_links])))
  if (!is.logical(x$implicit) || anyNA(x$implicit))
    stop("`family$dpar_specs$implicit` must be logical without missing values.")
  if (!is.numeric(x$lower))
    stop("`family$dpar_specs$lower` must be numeric.")

  TRUE
}


#' Get metadata for one distributional parameter
#'
#' @keywords internal
#' @noRd
get_dpar_spec = function(family, dpar) {
  spec = family$dpar_specs[family$dpar_specs$dpar == dpar, , drop = FALSE]
  if (nrow(spec) != 1)
    stop_github("Expected exactly one dpar specification for '", dpar, "'.")
  spec
}


#' JAGS function names for links and inverse links
#'
#' Empty strings represent the identity function, matching the historical mcp
#' convention used by generated code.
#'
#' @keywords internal
#' @noRd
get_link_str = function(link, inverse = FALSE) {
  if (!inverse)
    return(ifelse(link == "identity", "", link))

  switch(
    link,
    logit = "ilogit",
    probit = "phi",
    log = "exp",
    identity = ""
  )
}


#' R link or inverse-link function
#'
#' @keywords internal
#' @noRd
get_link_function = function(link, inverse = FALSE) {
  link_object = stats::make.link(link)
  if (inverse) link_object$linkinv else link_object$linkfun
}


#' @aliases is.mcpfamily
#' @describeIn mcpfamily Checks whether x is an `mcpfamily`.
#' @export
is.mcpfamily = function(x) {
  if (inherits(x, "mcpfamily") == FALSE)
    return(FALSE)

  assert_types(x, "family")
  assert_types(x$default_prior, "tibble", "data.frame")
  assert_types(x$dpars, "character")
  if (!is.null(x$dpar_specs)) {
    assert_dpar_specs(x$dpar_specs)
    assert_types(x$links, "character")
  }
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
