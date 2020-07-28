#' Bernoulli family for mcp
#'
#' @aliases bernoulli
#' @param link Link function.
#' @export
#'
bernoulli = function(link = "logit") {
  # Just copy binomial()
  family = binomial(link = link)
  family$family = "bernoulli"
  mcp_family(family)
}

#' Exponential family for mcp
#'
#' @aliases exponential
#' @param link Link function
#' @export
#'
exponential = function(link = "identity") {
  if (link != "identity")
    stop("Only link = 'identity' is currently supported for the exponential() family in mcp.")

  family = list(
    family = "exponential",
    link = "identity",  # on lambda
    linkfun = identity,  # on lambda
    linkinv = identity  # on lambda
  )
  class(family) = "family"
  family = mcp_family(family)
}


#' Add A family object to store link functions between R and JAGS.
#'
#' This will make more sense once more link functions / families are added.
#'
#' @aliases mcp_family
#' @keywords internal
#' @param family A family object, e.g., `binomial(link = "identity")`.
mcp_family = function(family) {
  if (family$link == "logit") {
    family$link_jags = "logit"
    family$link_r = "logit"  # included in mcp
    family$linkinv_jags = "ilogit"
    family$linkinv_r = "ilogit"  # included in mcp
  }

  if (family$link == "probit") {
    family$link_jags = "probit"
    family$link_r = "probit"  # included in mcp
    family$linkinv_jags = "phi"
    family$linkinv_r = "iprobit"  # included in mcp
  }

  if (family$link == "log") {
    family$link_jags = "log"
    family$link_r = "log"
    family$linkinv_jags = "exp"
    family$linkinv_r = "exp"
  }

  # Identity is just the absence of a function
  if (family$link == "identity") {
    family$link_jags = ""
    family$link_r = ""
    family$linkinv_jags = ""
    family$linkinv_r = ""
  }

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
#' @aliases iprobit
#' @param eta A vector of probits
#' @return A vector with same length as `mu`
#' @export
iprobit = stats::binomial(link = "probit")$linkinv


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

#' Asserts whether x contains non-numeric, decimal, or less-than-lower
#'
#' This is basically an extended `is.integer` with informative error messages.
#'
#' @aliases assert_integer
#' @keywords internal
#' @param x Value or vector
#' @param name Name to show in error message.
#' @param lower the smallest allowed value. lower = 1 checks for positive integers.
#'
assert_integer = function(x, name = NULL, lower = -Inf) {
  # Recode
  x = stats::na.omit(x)
  if (is.null(name))
    name = substitute(x)

  # Do checks
  greater_than = ifelse(lower == -Inf, " ", paste0(" >= ", lower, " "))
  if (!is.numeric(x))
    stop("Only integers", greater_than, "allowed for '", name, "'. Got ", x)

  if (!all(x == floor(x)) | !all(x >= lower))
    stop("Only integers", greater_than, "allowed for '", name, "'. Got ", x)

  TRUE
}

#' Asserts whether x is logical
#'
#' This is basically `is.logical` with informative error messages
#'
#' @aliases assert_logical
#' @keywords internal
#' @inheritParams assert_integer
#' @param max_length Maximum length of vector. Throws error if this is exceeded
#'
assert_logical = function(x, name = NULL, max_length = 1) {
  if (is.null(name))
    name = substitute(x)

  if (!is.logical(x))
    stop("`", name, "` must be TRUE or FALSE. Got ", x)

  if (length(x) > max_length)
    stop("`", name, "` must be at most length ", max_length, ". Got ", x)
}

#' Asserts whether x is one of a set of allowed values
#'
#' @aliases assert_charvec
#' @keywords internal
#' @inheritParams assert_integer
#' @param allowed Vector of allowed values
#'
assert_value = function(x, name = NULL, allowed = c()) {
  if (is.null(name))
    name = substitute(x)

  if (!(x %in% allowed)) {
    allowed[is.character(allowed)] = paste0("'", allowed[is.character(allowed)], "'")  # Add quotes for character values
    stop("`", name, "` must be one of ", paste0(allowed, collapse = ", "), ". Got ", x)
  }
}


#' Asserts whether x is an `mcpfit`
#'
#' @aliases assert_mcpfit
#' @keywords internal
#' @param x Object to be tested
#'
assert_mcpfit = function(x) {
  if (!is.mcpfit(x))
    stop("Expected `mcpfit` but got: ", class(x))
}


# Ask reminder questions for CRAN export
release_questions = function() {
  c(
    "Have you run the test of fits? options(test_mcp_fits = TRUE)",
    "Have you built the README plots and checked them? source('man/figures/make_README_plots.R')",
    "Have you re-built the site using pkgdown::build_site() AFTER deleting caches of articles?",
    "Have you checked all articles and plots after re-building the site?",
    "Have you run the script to insert the correct logo.png in the HTML meta?"
  )
}


#' Remove varying or population terms from a formula
#'
#' WARNING: removes response side from the formula
#'
#' @aliases remove_terms
#' @keywords internal
#' @param form A formula
#' @param remove Either "varying" or "population". These are removed.
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
remove_terms = function(form, remove) {
  assert_value(remove, allowed = c("varying", "population"))

  # Find terms with "|"
  attrs = attributes(stats::terms(form))
  term.labels = attrs$term.labels
  varying_bool = stringr::str_detect(term.labels, "\\|")

  # Add parenthesis back to them
  term.labels[varying_bool] = paste0("(", term.labels[varying_bool], ")")

  # Remove non-matching types
  if (remove == "varying") {
    term.labels = term.labels[!varying_bool]
    term.labels = c(attrs$intercept, term.labels)  # Add intercept indicator
  } else if (remove == "population") {
    term.labels = term.labels[varying_bool]
  }

  # Build formula from terms and return
  if (length(term.labels) == 0) {
    return(NULL)
  } else {
    formula_terms = paste0(term.labels, collapse = " + ")
    formula_str = paste0("~", formula_terms)
    return(stats::as.formula(formula_str, env=globalenv()))
  }
}


#' Takes any formula-like input and returns a formula
#' @aliases to_formula
#' @keywords internal
#' @param form Formula or character (with or without initial tilde/"~")
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
to_formula = function(form) {
  if (is.character(form)) {
    # Add tilde
    if (!stringr::str_detect(form, "^(\\s|)~")) {
      form = paste0("~", form)
    }
    form = stats::as.formula(form)
  } else if (!rlang::is_formula(form)) {
    stop("`form` must character or formula but got ", form)
  }

  return(form)
}

#' Expand samples with quantiles
#'
#' TO DO: implement using `fitted()` and `predict()` but avoid double-computing the samples? E.g.:
#' `get_quantiles2 = function(fit, quantiles, facet_by = NULL) {`
#'   `fitted(fit, probs = c(0.1, 0.5, 0.9), newdata = data.frame(x = c(11, 50, 100))) %>%`
#'   `tidyr::pivot_longer(tidyselect::starts_with("Q")) %>%`
#'   `dplyr::mutate(quantile = stringr::str_remove(name, "Q") %>% as.numeric() / 100)`
#' `}`
#'
#' @aliases get_quantiles
#' @keywords internal
#' @inheritParams plot.mcpfit
#' @param samples A tidybayes tibble
#' @param quantiles Vector of quantiles (0.0 to 1.0)
#' @param xvar An rlang::sym() with the name of the x-col in `samples`
#' @param yvar An rlang::sym() with the name of the response col in `samples`
#' @return A tidybayes long format tibble with the column "quantile"
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_quantiles = function(samples, quantiles, xvar, yvar, facet_by = NULL) {
  # Trick to declare no facet = common group for all
  if (is.null(facet_by))
    facet_by = xvar

  # Return data with added quantiles
  samples %>%
    tidyr::expand_grid(quantile = quantiles) %>%

    # Now compute the quantile for each parameter, quantile, and (optionally) facet:
    dplyr::group_by(!!xvar, .data$quantile) %>%
    dplyr::group_by(!!rlang::sym(facet_by), .add = TRUE) %>%
    dplyr::summarise(
      y = stats::quantile(!!yvar, probs = .data$quantile[1])
    )
}

# Hack to make R CMD pass for function geom_cp_density()
utils::globalVariables(c("value", "..scaled..", ".chain", "cp_name", "."))
