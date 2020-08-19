# ABOUT: These functions are used internally for for defensive programming.
# -----------------

# Synonym so that assert_types(x, "tibble", "formula") works
is.tibble = tibble::is_tibble
is.formula = rlang::is_formula

# Asserts whether x is an `mcpfit`
assert_mcpfit = function(x) {
  if (!is.mcpfit(x))
    stop("Expected `mcpfit` but got: ", class(x))
}

# Asserts whether x contains non-numeric, decimal, or less-than-lower
assert_integer = function(x, name = NULL, lower = -Inf) {
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
assert_logical = function(x, max_length = 1) {
  if (!is.logical(x))
    stop("`", substitute(x), "` must be logical (TRUE or FALSE). Got ", x)

  if (length(x) > max_length)
    stop("`", substitute(x), "` must be at most length ", max_length, ". Got length ", length(x))
}

# Asserts whether x is one of a set of allowed values
assert_value = function(x, allowed = c()) {
  if (!(x %in% allowed)) {
    allowed[is.character(allowed)] = paste0("'", allowed[is.character(allowed)], "'")  # Add quotes for character values
    stop("`", substitute(x), "` must be one of ", paste0(allowed, collapse = ", "), ". Got ", x)
  }
}


# Asserts whether x is one of a set of allowed types.
# e.g., `assert_types(vec, "numeric", "character", "foo")`
assert_types = function(x, ...) {
  types = list(...)

  # Test each function on x
  passed = logical(length(types))
  for (i in seq_along(types)) {
    is.type = eval(parse(text = paste0("is.", types[[i]])))  # From character to is.foo() function
    passed[i] = is.type(x)
  }

  # Return helpful error
  if (!any(passed == TRUE))
    stop("`", substitute(x), "` must be ", paste0(types, collapse = " or "), ". Got ", class(x))
}

# Asserts whether x is numeric in range
assert_numeric = function(x, lower = -Inf, upper = Inf) {
  if (!is.numeric(x))
    stop("`", substitute(x), "` must be numeric. Got ", class(x))
  if (any(x < lower) || any(x > upper))
    stop("`", substitute(x), "` contained value(s) outside the interval (", lower, ", ", upper, ").")
}

# Asserts ellipsis. `ellipsis` is a list and `allowed` is a character vector
assert_ellipsis = function(..., allowed = NULL) {
  assert_types(allowed, "null", "character")
  illegal_names = dplyr::setdiff(names(list(...)), allowed)
  if (length(illegal_names) > 0)
    stop("The following arguments are not accepted for this function: '", paste0(illegal_names, collapse = "', '"), "'")
}
