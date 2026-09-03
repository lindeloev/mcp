# ABOUT: These functions are used internally for defensive programming.
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
  is.numeric(x) &&
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
  }

  TRUE
}


# Asserts whether the data contains these cols and that all of them does not contained fail_funcs values
# This applies checks to multiple columns in data.
# Typical fail_funcs would be c(is.na, is.nan, is.infinite)
assert_data_cols = function(data, cols, fail_funcs = c()) {
  missing_cols = setdiff(cols, colnames(data))
  if (length(missing_cols) > 0)
    stop("These model terms are missing from the data: ", and_collapse(missing_cols))

  # Only work with the specified columns now
  data = data[, cols, drop = FALSE]
  for (fail_func in fail_funcs) {
    failed_cols = colnames(data)[unlist(lapply(data, function(x) any(fail_func(x))))]  # Character vector of columns that
    if (length(failed_cols) > 0) {
      # Fail with informative message
      fail_name = if (identical(fail_func, is.na)) "is.na" else
        if (identical(fail_func, is.nan)) "is.nan" else
          if (identical(fail_func, is.infinite)) "is.infinite" else "the validation check"
      stop("The column(s) ", and_collapse(failed_cols), " had values where ", fail_name, " was TRUE.")
    }
  }

  TRUE
}


# Validate data that define model structure rather than an individual response.
# Response columns are intentionally excluded because some families support
# missing responses as latent values during posterior sampling.
assert_model_data = function(data, par_x, rhs_vars = character(), group_cols = character()) {
  assert_data_cols(data, rhs_vars, fail_funcs = c(is.na, is.nan))
  assert_data_cols(
    data, unique(c(par_x, group_cols)),
    fail_funcs = c(is.na, is.nan, is.infinite)
  )

  invisible(TRUE)
}


# Validate response auxiliaries independently of whether responses themselves
# are observed. This makes the same family-level invariants available to model
# construction and R-side simulation.
assert_response_data = function(family, segments, data) {
  aux_columns = get_family_aux_columns(family, segments)
  response_columns = c(list(y = segments$y[1]), as.list(aux_columns))
  response_data = get_family_response_data(family, segments, data)
  y = if (segments$y[1] %in% names(data)) data[[segments$y[1]]] else rep(NA_real_, nrow(data))
  family$response$validate(y, response_data, response_columns)

  invisible(TRUE)
}


# Fixed residual standard deviations must be valid without relying on the
# numerical lower bound used for observation-level calculations.
assert_fixed_sigma = function(prior_table, predictors, family) {
  sigma_parameters = predictors$code_name[predictors$dpar == "sigma"]
  if (length(sigma_parameters) == 0)
    return(invisible(TRUE))
  if (get_dpar_spec(family, "sigma")$modeled)
    return(invisible(TRUE))

  fixed_sigma = prior_table$parameter %in% sigma_parameters & prior_table$kind == "constant"
  sigma_values = suppressWarnings(as.numeric(prior_table$value[fixed_sigma]))
  if (any(sigma_values <= 0)) {
    bad_parameters = prior_table$parameter[fixed_sigma][sigma_values <= 0]
    stop(
      "Fixed residual standard deviation parameter(s) must be positive: ",
      and_collapse(bad_parameters), "."
    )
  }

  invisible(TRUE)
}


# Change-point locations supplied to simulation must use the same ordering as
# the population-level parameters in the fitted model.
assert_ordered_population_cps = function(cps, args) {
  if (nrow(cps) < 2)
    return(invisible(NULL))

  locations = lapply(cps$name, function(name) args[[name]])
  n = max(lengths(locations))
  locations = do.call(cbind, lapply(locations, rep, length.out = n))
  if (any(locations[, -1, drop = FALSE] <= locations[, -ncol(locations), drop = FALSE]))
    stop("Population-level change points must remain strictly ordered.")

  invisible(NULL)
}


# Asserts whether `rel` is in a model
assert_rel = function(model) {
  has_rel = model %>%
    sapply(formula_to_char) %>%
    stringr::str_detect("rel\\(") %>%
    any()

  if (has_rel)
    stop(legacy_mcp_message(
      "`rel()` model terms were removed in mcp v0.4. Relative parameter estimates can be computed by subtracting posterior draws; there is no replacement for setting relative priors."
    ), call. = FALSE)

  TRUE
}


assert_typescale = function(type, scale) {
  rlang::arg_match0(type, c("predict", "fitted", "residuals", "loglik"))
  rlang::arg_match0(scale, c("response", "linear"))
  if (scale == "linear" && type != "fitted")
    stop("scale = 'linear' is only allowed when type = 'fitted'.")

  TRUE
}


# Check dpar with helpful errors and substitute NULL
assert_dpar = function(dpar, fit, type) {
  if (is.null(dpar))
    dpar = "epred"

  predictors = get_fit_model_tables(fit)$predictors
  allowed_dpars = unique(c(
    paste0(predictors$dpar, tidyr::replace_na(as.character(predictors$order), "")),   # "mu" "ar1" "ar2" "sigma", etc.
    fit$family$dpars[fit$family$dpars %notin% c("ar", "ma")],  # any model parameters that have no regression terms (~0)
    "epred"
  ))
  checkmate::assert_string(dpar)
  dpar = rlang::arg_match0(dpar, allowed_dpars)

  if (type != "fitted" && dpar != "epred")
    stop("dpar = '", dpar, "' is only allowed when type = 'fitted'.")

  dpar
}


#' Check sampled population- and group-level change-point locations
#'
#' @keywords internal
#' @noRd
#' @param draws An `mcmc.list` after group levels have been recovered.
#' @param cps The change-point table from `get_segment_tables()`.
#' @param x The observed change-point predictor.
assert_ordered_cp_draws = function(draws, cps, x) {
  if (is.null(draws) || nrow(cps) == 0)
    return(invisible(NULL))

  x_range = range(x, na.rm = TRUE)
  check_locations = function(locations, scope, check_range = FALSE) {
    if (any(!is.finite(locations)))
      stop("Sampled ", scope, " change points must be finite.")
    if (check_range && any(locations < x_range[1] | locations > x_range[2]))
      stop("Sampled ", scope, " change points must lie within the observed range of the change-point predictor.")
    if (ncol(locations) > 1 && any(locations[, -1, drop = FALSE] <= locations[, -ncol(locations), drop = FALSE]))
      stop("Sampled ", scope, " change points must remain strictly ordered in every draw.")
  }

  for (chain in draws) {
    missing_cps = setdiff(cps$name, colnames(chain))
    if (length(missing_cps) > 0)
      stop_github("Sampled output is missing change point(s): ", and_collapse(missing_cps), ".")

    population = chain[, cps$name, drop = FALSE]
    check_locations(population, "population-level", check_range = TRUE)

    varying = which(cps$varying)
    if (length(varying) == 0 || nrow(cps) < 2)
      next

    first_name = cps$group_name[varying[1]]
    first_prefix = paste0(first_name, "[")
    group_cols = colnames(chain)[startsWith(colnames(chain), first_prefix)]
    if (length(group_cols) == 0)
      stop_github("Sampled output is missing group-level change point `", first_name, "`.")
    group_levels = substring(group_cols, nchar(first_prefix) + 1L, nchar(group_cols) - 1L)
    for (group_level in group_levels) {
      realized = population
      for (j in varying) {
        group_col = paste0(cps$group_name[j], "[", group_level, "]")
        if (!group_col %in% colnames(chain))
          stop_github("Sampled output is missing group-level change point `", group_col, "`.")
        realized[, j] = realized[, j] + chain[, group_col]
      }
      check_locations(realized, "group-level")
    }
  }

  invisible(NULL)
}
