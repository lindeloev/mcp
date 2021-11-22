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

  TRUE
}

#' @keywords internal
is.mcpmodel = function(x) {
  assert_types(x, "list", len = c(1, Inf))

  for (segment in x) {
    if (!inherits(segment, "formula"))
      stop("All segments must be formulas. This segment is not a formula: ", segment)
  }

  TRUE
}


# Is this a continuous vector?
is_continuous = function(x) {
  is.numeric(x) &
    length(unique(stats::na.omit(x))) > 2
}


#' @keywords internal
stop_github = function(...) {
  stop("This looks like an internal error in mcp. To fix this for you and other users, please raise an issue on https://github.com/lindeloev/mcp/issues with the minimum data/code that reproduces this error:\n", ...)
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

  TRUE
}


# Asserts whether x is logical
assert_logical = function(x, len = 1) {
  assert_length(x, len, name = substitute(x))

  if (!is.logical(x))
    stop("`", substitute(x), "` must be logical (TRUE or FALSE). Got ", x)

  TRUE
}


# Asserts whether x is one of a set of allowed values
assert_value = function(x, allowed = c(), len = 1) {
  assert_length(x, len, name = substitute(x))

  if (x %notin% allowed) {
    allowed[is.character(allowed)] = paste0("'", allowed[is.character(allowed)], "'")  # Add quotes for character values
    if (length(allowed) == 1) {
      stop("`", substitute(x), "` must be ", allowed, ". Got ", paste0(x, collapse = ", "))
    } else {
      stop("`", substitute(x), "` must be one of ", paste0(allowed, collapse = " or "), ". Got ", paste0(x, collapse = ", "))
    }
  }

  TRUE
}


# Asserts whether x is one of a set of allowed types.
# e.g., `assert_types(vec, "numeric", "character", "foo")`
assert_types = function(x, ..., len = NULL) {
  assert_length(x, len, name = substitute(x))
  types = list(...)

  # Test each function on x
  passed = logical(length(types))
  for (i in seq_along(types)) {
    is.type = get(paste0("is.", types[[i]]))  # From character to is.foo() function
    passed[i] = is.type(x)
  }

  # Return helpful error
  if (!any(passed == TRUE))
    if (length(types) == 1) {
      stop("`", substitute(x), "` must be ", types, ". Got ", paste0(class(x), collapse = ", "))
    } else {
      stop("`", substitute(x), "` must be one of ", paste0(types, collapse = " or "), ". Got ", paste0(class(x), collapse = ", "))
    }

  TRUE
}


# Asserts whether x is numeric in range
assert_numeric = function(x, lower = -Inf, upper = Inf, len = NULL) {
  assert_length(x, len, name = substitute(x))
  if (!is.numeric(x))
    stop("`", substitute(x), "` must be numeric. Got ", class(x))
  if (any(x < lower) || any(x > upper))
    stop("`", substitute(x), "` contained value(s) outside the interval (", lower, ", ", upper, ").")

  TRUE
}


# Asserts ellipsis. `ellipsis` is a list and `allowed` is a character vector
assert_ellipsis = function(..., allowed = NULL) {
  assert_types(allowed, "null", "character")

  illegal_names = dplyr::setdiff(names(list(...)), allowed)
  if (length(illegal_names) > 0)
    stop("The following arguments are not accepted for this function: ", and_collapse(illegal_names))

  TRUE
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

  TRUE
}


# Asserts whether matrix x is rank deficient.
assert_rank = function(x, segment, dpar) {
  QR = qr(x)
  if (QR$rank < ncol(x)) {
    bad_cols = colnames(x)[QR$pivot[(QR$rank+1):ncol(x)]]
    stop("These terms are perfectly colinear with other terms for ", dpar, " in segment ", segment, ": ", and_collapse(bad_cols), " (the design matrix is rank deficient). Consider checking the data and/or the model.")
    return(bad_cols)
  }

  TRUE
}


# Asserts whether the data contains these cols and that all of them does not contained fail_funcs values
# This is like assert_types(), but applied to multiple columns in data
# Typical fail_funcs would be c(is.na, is.nan, is.infinite)
assert_data_cols = function(data, cols, fail_funcs = c()) {
  missing_cols = cols[(cols %in% colnames(data)) == FALSE]
  if (length(missing_cols) > 0)
    stop("These model terms are missing from the data: ", and_collapse(missing_cols))

  # Only work with the specified columns now
  data = data[, cols]
  for (fail_func in fail_funcs) {
    failed_cols = colnames(data)[unlist(lapply(data, function(x) any(fail_func(x))))]  # Character vector of columns that
    if (length(failed_cols) > 0)
      stop("The column(s) ", and_collapse(failed_cols), " had values where ", as.character(substitute(fail_func)), " was TRUE.")
  }

  TRUE
}


# Asserts whether `rel` is in a model
assert_rel = function(model) {
  has_rel = model %>%
    sapply(formula_to_char) %>%
    stringr::str_detect("rel\\(") %>%
    any()

  if (has_rel)
    stop("rel() for model terms was deprecated in mcp 0.4.0. Relative parameter estimates can be computed by subtracting posterior samples. There is no replacement wrt setting priors.")

  TRUE
}
