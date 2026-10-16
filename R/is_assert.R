# ABOUT: These functions are used internally for defensive programming.
# -----------------

# Check if object is an mcp family
is_family = function(x) {
  if (inherits(x, "family") == FALSE)
    return(FALSE)

  checkmate::assert_string(x$family)
  checkmate::assert_string(x$link)

  TRUE
}

# Check if object is a valid mcp model formula list
is_mcpmodel = function(x) {
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


# Stop with an informative error directing users to GitHub issues
stop_github = function(...) {
  stop("This looks like an internal error in mcp. To fix this for you and other users, please raise an issue on https://github.com/lindeloev/mcp/issues with the minimum data/code that reproduces this error:\n", ...)
}


# Require offset() in every segment where a dpar applies if it is in any of them. Offsets usually
# describe the observation (e.g., exposure), so omitting one in a segment is likely a mistake.
# Dpars apply in all segments (joined when omitted); AR/MA only in segments that include them.
assert_offsets = function(components, family) {
  for (dpar in c(family$dpars, "ar", "ma")) {
    is_arma = dpar %in% c("ar", "ma")
    forms = lapply(components, function(x) if (is_arma) x[[dpar]]$form else x[[dpar]])
    applies = if (is_arma) !vapply(forms, is.null, logical(1)) else rep(TRUE, length(forms))
    has_offset = vapply(forms, function(form) !is.null(form) && !is.null(attr(stats::terms(form), "offset")), logical(1))
    missing_segments = which(applies & !has_offset)
    if (any(has_offset) && length(missing_segments) > 0)
      stop(
        "`offset()` must be included in all segments of `", dpar, "` but is missing in segment ",
        and_collapse(as.character(missing_segments)), ". If you have a use case for omitting it, ",
        "please request it on https://github.com/lindeloev/mcp/issues",
        call. = FALSE
      )
  }
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


# Assert that column names do not collide with reserved mcp output names
# - extra: Additional reserved names, e.g., the evaluated dpar column (`.sigma`).
assert_reserved_output_namespace = function(col_names, context = "data", extra = NULL) {
  reserved_output = c(
    ".mean", ".sd", ".lower", ".upper",
    ".epred", ".linpred", ".prediction", ".residual", ".loglik",
    ".draw", ".chain", ".iteration", ".row",
    extra
  )
  collisions = intersect(col_names, reserved_output)
  if (length(collisions) > 0) {
    target = if (context == "data") "Data column name(s)" else paste0("`", context, "` column name(s)")
    action = if (context == "data") "before fitting." else "before evaluation."
    stop(
      target, " collide with mcp's reserved output namespace: ",
      and_collapse(paste0("'", collisions, "'")),
      ". Rename them ", action
    )
  }
  invisible(TRUE)
}


# Assert that column names do not collide with reserved mcp offset names
assert_reserved_offset_namespace = function(col_names, context = "data") {
  collisions = grep("^offset(_[0-9]+)?$", col_names, value = TRUE)
  if (length(collisions) > 0) {
    target = if (context == "data") "Data column name(s)" else paste0("`", context, "` column name(s)")
    action = if (context == "data") "before fitting." else "before evaluation."
    stop(
      target, " collide with mcp's reserved offset namespace: ",
      and_collapse(paste0("'", collisions, "'")),
      ". Rename them ", action
    )
  }
  invisible(TRUE)
}


# Validate response auxiliaries independently of whether responses themselves
# are observed. This makes the same family-level invariants available to model
# construction and R-side simulation.
assert_response_data = function(family, segments, data, response_required = TRUE, aux_required = NULL) {
  aux_columns = get_family_aux_columns(family, segments)
  response_columns = c(list(y = segments$y[1]), as.list(aux_columns))
  response_data = get_family_response_data(family, segments, data)

  # For auxiliaries not required by the operation, supply valid placeholders
  if (!is.null(aux_required)) {
    for (aux_name in names(aux_columns)) {
      if (aux_name %notin% aux_required) {
        response_data[[aux_name]] = rep(1L, nrow(data))
      }
    }
  }

  y = if (response_required && segments$y[1] %in% names(data)) data[[segments$y[1]]] else rep(NA_real_, nrow(data))
  family$response$validate(y, response_data, response_columns)

  invisible(TRUE)
}


# Fixed constant dpars must be within their support (e.g., sigma > 0) without relying
# on the numerical floor used for observation-level calculations.
assert_fixed_dpars = function(prior_table, predictors, family) {
  bounds = get_constant_dpar_bounds(family)
  for (dpar in names(bounds)) {
    fixed = prior_table$parameter %in% predictors$name[predictors$dpar == dpar] & prior_table$kind == "constant"
    values = suppressWarnings(as.numeric(prior_table$value[fixed]))
    if (any(values <= bounds[[dpar]])) {
      bad_parameters = prior_table$parameter[fixed][values <= bounds[[dpar]]]
      stop(
        "Fixed `", dpar, "` parameter(s) must be greater than ", bounds[[dpar]], ": ",
        and_collapse(bad_parameters), "."
      )
    }
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
    stringr::str_detect("(^|[^A-Za-z0-9._])rel\\(") %>%
    any()

  if (has_rel)
    stop(legacy_mcp_message(
      "`rel()` model terms were removed in mcp v0.4. Relative parameter estimates can be computed by subtracting posterior draws; there is no replacement for setting relative priors."
    ), call. = FALSE)

  TRUE
}


# Width of a central interval. 0 and 1 are degenerate intervals, so they are not allowed.
# - false_ok: Whether `FALSE` (no interval) is allowed.
assert_width = function(width, name = "width", false_ok = FALSE) {
  if (false_ok && isFALSE(width))
    return(invisible(width))
  checkmate::assert_number(width, .var.name = name)
  if (width <= 0 || width >= 1)
    stop("`", name, "` must be strictly between 0 and 1", if (false_ok) ", or `FALSE`", ".", call. = FALSE)
  invisible(width)
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

