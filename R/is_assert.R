# ABOUT: These functions are used internally for for defensive programming.
# -----------------

# Synonym so that assert_types(x, "tibble", "formula") works
#' @keywords internal
is.tibble = tibble::is_tibble

#' @keywords internal
is.formula = rlang::is_formula

#' @keywords internal
is.family = function(x) {
  if (inherits(x, "family") == FALSE)
    return(FALSE)

  assert_types(x$family, "character", len = 1)
  assert_types(x$link, "character", len = 1)

  return(TRUE)
}


#' @keywords internal
stop_github = function(...) {
  stop(..., ". This looks like an internal error in mcp. Please raise an issue on GitHub: https://github.com/lindeloev/mcp/issues")
}


# Asserts whether x contains non-numeric, decimal, or less-than-lower
assert_integer = function(x, name = NULL, lower = -Inf, len = NULL) {
  assert_length(x, len, name = substitute(x))

  # Recode
  if (is.null(name))
    name = substitute(x)
  x = stats::na.omit(x)

  # Do checks
  greater_than = ifelse(lower == -Inf, " ", paste0(" >= ", lower, " "))
  if (!is.numeric(x))
    stop("Only integers", greater_than, "allowed for '", name, "'. Got ", x)

  if (!all(x == floor(x)) || !all(x >= lower))
    stop("Only integers", greater_than, "allowed for '", name, "'. Got ", x)

  return(TRUE)
}


# Asserts whether x is logical
assert_logical = function(x, len = 1) {
  assert_length(x, len, name = substitute(x))

  if (!is.logical(x))
    stop("`", substitute(x), "` must be logical (TRUE or FALSE). Got ", x)

  return(TRUE)
}


# Asserts whether x is one of a set of allowed values
assert_value = function(x, allowed = c(), len = 1) {
  assert_length(x, len, name = substitute(x))

  if (!(x %in% allowed)) {
    allowed[is.character(allowed)] = paste0("'", allowed[is.character(allowed)], "'")  # Add quotes for character values
    stop("`", substitute(x), "` must be one of ", paste0(allowed, collapse = " or "), ". Got ", x)
  }

  return(TRUE)
}


# Asserts whether x is one of a set of allowed types.
# e.g., `assert_types(vec, "numeric", "character", "foo")`
assert_types = function(x, ..., len = NULL) {
  assert_length(x, len, name = substitute(x))
  types = list(...)

  # Test each function on x
  passed = logical(length(types))
  for (i in seq_along(types)) {
    is.type = eval(parse(text = paste0("is.", types[[i]])))  # From character to is.foo() function
    passed[i] = is.type(x)
  }

  # Return helpful error
  if (!any(passed == TRUE))
    stop("`", substitute(x), "` must be one of ", paste0(types, collapse = " or "), ". Got ", class(x))

  return(TRUE)
}


# Asserts whether x is numeric in range
assert_numeric = function(x, lower = -Inf, upper = Inf, len = NULL) {
  assert_length(x, len, name = substitute(x))
  if (!is.numeric(x))
    stop("`", substitute(x), "` must be numeric. Got ", class(x))
  if (any(x < lower) || any(x > upper))
    stop("`", substitute(x), "` contained value(s) outside the interval (", lower, ", ", upper, ").")

  return(TRUE)
}


# Asserts ellipsis. `ellipsis` is a list and `allowed` is a character vector
assert_ellipsis = function(..., allowed = NULL) {
  assert_types(allowed, "null", "character")

  illegal_names = dplyr::setdiff(names(list(...)), allowed)
  if (length(illegal_names) > 0)
    stop("The following arguments are not accepted for this function: ", and_collapse(illegal_names))

  return(TRUE)
}


# Asserts whether the length of x matches len or len = c(lower, upper)
assert_length = function(x, len = NULL, name = NULL) {
  if (is.null(name))
    name = substitute(x)

  if (is.null(len) == FALSE) {
    if (length(len) == 1 && length(x) != len)
      stop("`", name, "` should have length ", len, " but has length ", length(x))

    if (length(len) == 2 && (length(x) < len[1] | length(x) > len[2]))
      stop("`", name, "` should have a length between ", len[1], " and ", len[2], " but has length ", length(x))
  }

  return(TRUE)
}


# Asserts whether matrix x is rank deficient.
assert_rank = function(x, segment) {
  QR = qr(x)
  if (QR$rank < ncol(x)) {
    bad_cols = colnames(x)[QR$pivot[(QR$rank+1):ncol(x)]]
    stop("These terms are perfectly colinear with other terms in segment ", segment, ": ", and_collapse(bad_cols), " (the design matix is rank deficient). Consider checking the data and/or the model.")
    return(bad_cols)
  }

  TRUE
}
