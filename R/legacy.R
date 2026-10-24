# Compatibility and Deprecation Detections for Legacy Code & Fits -------------
# Functions are prefixed `legacy_` so that call sites show what to remove in mcp 1.0.
# The S3 methods `$.mcpfit` and `[[.mcpfit` and the defunct `tidy_samples()` keep their names.

# Check if an mcpfit object was created with mcp < 0.4.0, which did not record `fit$version`
legacy_is_mcpfit = function(fit) {
  is.null(fit$version)
}


# Message for code that requires substantial migration from mcp v0.3.4
legacy_mcp_message = function(problem) {
  paste0(
    problem, "\n",
    "This appears to use mcp v0.3.4 code or an mcpfit saved with v0.3.4.\n",
    "To reproduce the old workflow while you migrate, downgrade to mcp v0.3.4:\n",
    "  remotes::install_github(\"lindeloev/mcp@v0.3.4\")\n",
    "mcp v0.3.4 does not include important v0.4.0 bug fixes."
  )
}


# Stop if an mcpfit cannot be used by mcp v0.4
check_mcpfit_version = function(fit) {
  checkmate::assert_class(fit, "mcpfit")
  if (legacy_is_mcpfit(fit))
    stop(legacy_mcp_message(
      "This `mcpfit` object lacks the internal structures introduced in mcp v0.4."
    ), call. = FALSE)
  TRUE
}


# Stop on the most common v0.3.4 parameter names
legacy_parameter_names = function(x, context) {
  old_names = x[stringr::str_detect(x, "(^|[^[:alnum:]_])int_[0-9]+($|[^[:alnum:]_])")]
  if (length(old_names) > 0)
    stop(legacy_mcp_message(paste0(
      "`", context, "` uses the v0.3.4 parameter name `int_i`; use `Intercept_i` in v0.4, which also changed syntax for some other variable names."
    )), call. = FALSE)

  old_names = x[stringr::str_detect(x, "(^|[^[:alnum:]_])cp_[0-9]+_sd($|[^[:alnum:]_])")]
  if (length(old_names) > 0)
    stop(
      "`", context, "` uses the v0.3.4 parameter name `cp_i_sd`; use `cp_i_<group>_sd` in v0.4, e.g., `cp_1_id_sd` for `(1 | id)`.",
      call. = FALSE
    )
}


# Translate deprecated left-hand sides in segment 2+ formulas to cp() syntax, with a warning
# E.g., `1 + (1 | id) ~ 0 + x` and `y ~ 1 + (1 | id) ~ 0 + x` --> `~ cp(1 + (1 | id)) + 0 + x`.
legacy_cp_lhs = function(form, i, response) {
  lhs = form[[2]]
  rhs = form[[3]]
  if (is_call_to(rhs, "~"))
    stop("Error in segment ", i, " (change point): empty change point term.", call. = FALSE)

  # `response ~ cp ~ rhs`
  if (is_call_to(lhs, "~")) {
    if (is_call_to(lhs[[2]], "~"))
      stop("Error in segment ", i, ": Got none or more than two ~ in a segment formula.", call. = FALSE)
    if (!identical(lhs[[2]], response))
      stop("Error in segment ", i, ": the response must be the same as in segment 1. Got `", deparse1(lhs[[2]]), "`.", call. = FALSE)
    lhs = lhs[[3]]
  }
  if (!is.null(split_cp(rhs)$cp))
    stop("Error in segment ", i, ": the change point is defined both in `cp()` and on the left-hand side. Use only `cp()`.", call. = FALSE)

  # If using the old lhs-cp, transform to rhs-cp during soft deprecation.
  new_rhs = if (identical(lhs, 1)) rhs else sum_exprs(c(call("cp", lhs), unpack_additive(rhs)))
  new_form = stats::as.formula(call("~", new_rhs), env = environment(form))
  warning(
    "Left-hand sides in segment 2+ formulas were deprecated in mcp v0.4.0. Use `cp()` for change points: ",
    "replace `", formula_to_char(form), "` with `", formula_to_char(new_form), "` in segment ", i, ".",
    call. = FALSE
  )
  new_form
}


# Warn on deprecated `which_y` in `...` and return it as a `dpar` (NULL if absent). Other
# arguments in `...` are errors. In mcp 1.0, replace with `rlang::check_dots_empty()`.
# - dots: `list(...)` of the calling function.
# - func_name: Name of calling function for warning message context.
legacy_which_y = function(dots, func_name) {
  if ("which_y" %in% names(dots)) {
    msg = if (func_name == "plot") {
      "`which_y` was deprecated in mcp v0.4.0. Use `plot_dpar()` instead."
    } else {
      "`which_y` was deprecated in mcp v0.4.0. Use `dpar` instead."
    }
    warning(msg, call. = FALSE)
  }

  unknown = setdiff(rlang::names2(dots), "which_y")
  if (length(unknown) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(ifelse(nzchar(unknown), unknown, "<unnamed>")), call. = FALSE)

  # "ct" (central tendency) was the v0.3.4 default
  if (identical(dots$which_y, "ct")) "epred" else dots$which_y
}


# Warn that the deprecated `cores` argument in mcp() is ignored
legacy_cores = function(cores, env = rlang::caller_env(), user_env = rlang::caller_env(2)) {
  if (!lifecycle::is_present(cores))
    return(invisible(NULL))
  if (!identical(cores, "all"))
    checkmate::assert_int(cores, lower = 1)

  cores_details = paste0(
    "`cores` is ignored. Parallel processing is now controlled by the active ",
    "future plan. For example, call ",
    "`future::plan(future::multisession, workers = 3)` before `mcp()`, and ",
    "`future::plan(future::sequential)` when the workers are no longer needed."
  )
  if (identical(cores, "all") || (is.numeric(cores) && cores > 1))
    cores_details = paste0(
      "Setting `cores` above one no longer enables parallel processing. ",
      cores_details
    )

  lifecycle::deprecate_warn(
    when = "0.4.0",
    what = "mcp(cores)",
    details = cores_details,
    env = env,
    user_env = user_env
  )
}


# Resolve an argument renamed in mcp v0.4.0 and return the value to use.
# - new, old: Values of the new and the deprecated argument, e.g., `ndraws` and `nsamples`.
# - new_missing: `missing(new)` in the calling function.
# - what: Name of the calling function, e.g., "fitted.mcpfit".
# - new_name, old_name: Argument names, e.g., "ndraws" and "nsamples".
# - deprecate: `lifecycle::deprecate_soft` or `lifecycle::deprecate_warn`.
legacy_arg = function(new, old, new_missing, what, new_name, old_name,
                      deprecate = lifecycle::deprecate_soft, details = NULL,
                      env = rlang::caller_env(), user_env = rlang::caller_env(2)) {
  if (!lifecycle::is_present(old))
    return(new)
  if (!new_missing)
    stop("Use only one of `", new_name, "` and deprecated `", old_name, "`.", call. = FALSE)

  deprecate(
    "0.4.0",
    paste0(what, "(", old_name, ")"),
    paste0(what, "(", new_name, ")"),
    details = details,
    env = env,
    user_env = user_env
  )
  old
}


# Resolve deprecated `probs` (from v0.3.4) to `width`. Only a single central interval can be translated.
legacy_width = function(width, probs, width_missing, what,
                        env = rlang::caller_env(),
                        user_env = rlang::caller_env(2)) {
  if (!lifecycle::is_present(probs))
    return(width)

  probs = legacy_arg(
    width, probs, width_missing, what, "width", "probs", lifecycle::deprecate_warn,
    details = "The interval columns are now `.lower` and `.upper`. For other quantiles, summarise the draws from `summary = FALSE`.",
    env = env, user_env = user_env
  )
  legacy_probs_to_width(probs, "probs")
}


# Translate deprecated quantile arguments (`TRUE`, `FALSE`, or a symmetric pair) to a central interval `width`.
legacy_probs_to_width = function(probs, old_name) {
  if (isTRUE(probs)) {
    0.95
  } else if (isFALSE(probs)) {
    FALSE
  } else if (is.numeric(probs) && length(probs) == 2 && isTRUE(all.equal(sum(probs), 1))) {
    abs(diff(probs))
  } else {
    stop("Deprecated `", old_name, "` can only be translated to a central interval, e.g., `", old_name, " = c(0.025, 0.975)`. Use `width` instead, or summarise the draws from `summary = FALSE` for other quantiles.", call. = FALSE)
  }
}


# Resolve deprecated `q_fit` and `q_predict` (from v0.3.4) in plot() to `interval` and `width`.
# Returns `list(interval, width)`.
legacy_plot_interval = function(interval, width, q_fit, q_predict, new_missing, what,
                                env = rlang::caller_env(), user_env = rlang::caller_env(2)) {
  old = list(fitted = q_fit, predict = q_predict)
  old = Filter(lifecycle::is_present, old)
  if (length(old) == 0)
    return(list(interval = interval, width = width))
  if (!new_missing)
    stop("Use only `interval` and `width`, not deprecated `q_fit` or `q_predict`.", call. = FALSE)

  old_names = c(fitted = "q_fit", predict = "q_predict")[names(old)]
  for (old_name in old_names)
    lifecycle::deprecate_warn(
      "0.4.0", paste0(what, "(", old_name, ")"), I(paste0("`", what, "(interval, width)`")),
      details = "For example, `q_fit = TRUE, q_predict = c(0.1, 0.9)` is now `interval = c(\"fitted\", \"predict\"), width = 0.8` (both intervals share `width`).",
      env = env, user_env = user_env
    )

  widths = mapply(legacy_probs_to_width, old, old_names)
  widths = widths[widths != FALSE]
  if (length(unique(widths)) > 1)
    stop("`q_fit` and `q_predict` imply different interval widths, but `width` applies to both intervals. Use one `width`, or add the other interval as a custom ggplot2 layer.", call. = FALSE)
  if (length(widths) == 0)
    return(list(interval = "none", width = width))
  list(interval = names(widths), width = unname(widths[1]))
}


# Warn when a binomial prediction relies on the pre-v0.4 proportion default.
legacy_binomial_rate = function(fit, newdata, what) {
  if (fit$family$family != "binomial")
    return(invisible(NULL))
  if (is.null(newdata))
    newdata = fit$data
  trials = newdata[[mcp_columns(fit)$trials]]
  if (any(trials != 1, na.rm = TRUE)) {
    rlang::warn(
      paste0("`", what, "()` now returns counts for binomial models by default. ",
             "Before mcp 0.4 it returned proportions. Set `rate = TRUE` for ",
             "proportions."),
      .frequency = "once", .frequency_id = paste0("mcp_", what, "_rate")
    )
  }
}


# Signal deprecated and removed mcpfit fields from v0.3.4. LEGACY: remove `$.mcpfit` and `[[.mcpfit` in mcp 1.0.
#' @export
`$.mcpfit` = function(x, name) {
  if (name == "fit" && is.null(.subset2(x, "fit"))) {
    warning("`mcp_example()` now returns an `mcpfit` object directly instead of a list with `$fit`. Returning the object itself.", call. = FALSE)
    return(x)
  }
  if (name %in% c("mcmc_post", "mcmc_prior")) {
    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = I(paste0("fit$", name)),
      with = I("as_draws(fit) or coda::as.mcmc(fit)")
    )
  }
  legacy_mcpfit_field(name)
  .subset2(x, name)
}

#' @export
`[[.mcpfit` = function(x, i, exact = TRUE) {
  if (is.character(i) && i %in% c("mcmc_post", "mcmc_prior")) {
    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = I(paste0("fit$", i)),
      with = I("as_draws(fit) or coda::as.mcmc(fit)")
    )
  }
  if (is.character(i))
    legacy_mcpfit_field(i)
  .subset2(x, i)
}


# Stop with a pointer to the replacement when accessing mcpfit fields removed after v0.3.4
legacy_mcpfit_field = function(name) {
  removed = c(pars = "mcp_pars(fit) and mcp_columns(fit)", mcmc_loglik = "log_lik(fit)")
  if (name %in% names(removed))
    lifecycle::deprecate_stop("0.4.0", I(paste0("fit$", name)), I(removed[[name]]))
}


# Defunct internal helper for MCMC draw extraction, advertised as `mcp:::tidy_samples()` in v0.3.4. LEGACY: remove in mcp 1.0.
tidy_samples = function(...) {
  lifecycle::deprecate_stop("0.4.0", "tidy_samples()", I("as_draws_df() or tidybayes::spread_draws()"))
}


#' Defunct functions in mcp
#'
#' These functions were exported in mcp v0.3.4 and removed in v0.4.0. They now
#' throw an error pointing to their replacements.
#'
#' * `criterion()`: Use [loo()] or [waic()].
#' * `get_segment_table()`: Use `fit = mcp(..., sample = FALSE)` followed by [mcp_pars()] or [mcp_columns()].
#'
#' @param ... Ignored.
#' @name mcp-defunct
#' @keywords internal
NULL

# LEGACY: remove in mcp 1.0.
#' @rdname mcp-defunct
#' @export
criterion = function(...) {
  lifecycle::deprecate_stop("0.4.0", "criterion()", I("`loo(fit)` or `waic(fit)`"))
}

# LEGACY: remove in mcp 1.0.
#' @rdname mcp-defunct
#' @export
get_segment_table = function(...) {
  lifecycle::deprecate_stop("0.4.0", "get_segment_table()", I("`mcp_pars(fit)` or `mcp_columns(fit)` on `fit = mcp(..., sample = FALSE)`"))
}


# Renamed to is_mcpfit() for consistency with `{posterior}`'s is_draws(). LEGACY: remove in mcp 1.0.
#' @rdname is_mcpfit
#' @export
is.mcpfit = function(x) {
  lifecycle::deprecate_warn("0.4.0", "is.mcpfit()", "is_mcpfit()")
  is_mcpfit(x)
}
