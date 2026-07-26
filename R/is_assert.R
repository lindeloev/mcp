# ABOUT: These functions are used internally for for defensive programming.
# -----------------

#' @keywords internal
is.family = function(x) {
  if (inherits(x, "family") == FALSE)
    return(FALSE)

  checkmate::assert_string(x$family)
  checkmate::assert_string(x$link)

  TRUE
}

#' @keywords internal
is.mcpmodel = function(x) {
  checkmate::assert_list(x, min.len = 1)

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
# This applies checks to multiple columns in data.
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


assert_typescale = function(type, scale) {
  rlang::arg_match0(type, c("predict", "fitted", "residuals", "loglik"))
  rlang::arg_match0(scale, c("response", "linear"))
  if (scale == "linear" & type != "fitted")
    stop("scale = 'linear' is only allowed when type = 'fitted'.")

  TRUE
}


# Check dpar with helpful errors and substitute NULL
assert_dpar = function(dpar, fit, type) {
  if (is.null(dpar))
    dpar = "epred"

  rhs_table = fit$.internal$rhs_table
  allowed_dpars = unique(c(
    paste0(rhs_table$dpar, tidyr::replace_na(as.character(rhs_table$order), "")),   # "mu" "ar1" "ar2" "sigma", etc.
    fit$family$dpars[fit$family$dpars %notin% c("ar", "ma")],  # any model parameters that have no regression terms (~0)
    "epred"
  ))
  checkmate::assert_string(dpar)
  dpar = rlang::arg_match0(dpar, allowed_dpars)

  if (type != "fitted" & dpar != "epred")
    stop("dpar = '", dpar, "' is only allowed when type = 'fitted'.")

  dpar
}
