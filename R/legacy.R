# Compatibility and Deprecation Detections for Legacy Code & Fits -------------

# Check if an mcpfit object was created with mcp < 0.4.0
is_legacy_mcpfit = function(fit) {
  tables = fit$.internal$model_tables
  !is.mcpfamily(fit$family) || is.null(tables) || is.null(tables$parameters) ||
    !"role" %in% names(tables$parameters)
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
  if (is_legacy_mcpfit(fit))
    stop(legacy_mcp_message(
      "This `mcpfit` object lacks the internal structures introduced in mcp v0.4."
    ), call. = FALSE)
  TRUE
}


# Stop on the most common v0.3.4 parameter names
check_legacy_parameter_names = function(x, context) {
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
translate_cp_lhs = function(form, i, response) {
  lhs = form[[2]]
  rhs = form[[3]]
  if (is_tilde_call(rhs))
    stop("Error in segment ", i, " (change point): empty change point term.", call. = FALSE)

  # `response ~ cp ~ rhs`
  if (is_tilde_call(lhs)) {
    if (is_tilde_call(lhs[[2]]))
      stop("Error in segment ", i, ": Got none or more than two ~ in a segment formula.", call. = FALSE)
    if (!identical(lhs[[2]], response))
      stop("Error in segment ", i, ": the response must be the same as in segment 1. Got `", deparse1(lhs[[2]]), "`.", call. = FALSE)
    lhs = lhs[[3]]
  }
  if (!is.null(split_cp(rhs)$cp))
    stop("Error in segment ", i, ": the change point is defined both in `cp()` and on the left-hand side. Use only `cp()`.", call. = FALSE)

  # If using the old lhs-cp, transform to rhs-cp during soft deprecation.
  new_rhs = if (identical(lhs, 1)) rhs else Reduce(function(a, b) call("+", a, b), c(call("cp", lhs), unpack_additive(rhs)))
  new_form = stats::as.formula(call("~", new_rhs), env = environment(form))
  warning(
    "Left-hand sides in segment 2+ formulas were deprecated in mcp v0.4.0. Use `cp()` for change points: ",
    "replace `", formula_to_char(form), "` with `", formula_to_char(new_form), "` in segment ", i, ".",
    call. = FALSE
  )
  new_form
}


# Warn on deprecated which_y argument and return it as a `dpar` (NULL if absent)
# - args: List of dots arguments passed to a function.
# - func_name: Name of calling function for warning message context.
warn_which_y = function(args, func_name) {
  if ("which_y" %in% names(args)) {
    msg = if (func_name == "plot") {
      "`which_y` was deprecated in mcp v0.4.0. Use `plot_dpar()` instead."
    } else {
      "`which_y` was deprecated in mcp v0.4.0. Use `dpar` instead."
    }
    warning(msg, call. = FALSE)
  }

  # "ct" (central tendency) was the v0.3.4 default
  if (identical(args$which_y, "ct")) "epred" else args$which_y
}


# Resolve the deprecated `nsamples` argument
resolve_ndraws = function(ndraws, nsamples, ndraws_missing, what,
                         env = rlang::caller_env(),
                         user_env = rlang::caller_env(2)) {
  if (lifecycle::is_present(nsamples)) {
    lifecycle::deprecate_soft(
      "0.4.0",
      paste0(what, "(nsamples)"),
      paste0(what, "(ndraws)"),
      env = env,
      user_env = user_env
    )
    if (!ndraws_missing)
      stop("Use only one of `ndraws` and deprecated `nsamples`.")
    ndraws = nsamples
  }
  ndraws
}


# Resolve the deprecated `samples_format` argument
resolve_draws_format = function(draws_format, samples_format, draws_format_missing, what,
                                env = rlang::caller_env(),
                                user_env = rlang::caller_env(2)) {
  if (lifecycle::is_present(samples_format)) {
    lifecycle::deprecate_soft(
      "0.4.0",
      paste0(what, "(samples_format)"),
      paste0(what, "(draws_format)"),
      env = env,
      user_env = user_env
    )
    if (!draws_format_missing)
      stop("Use only one of `draws_format` and deprecated `samples_format`.")
    draws_format = samples_format
  }
  rlang::arg_match0(draws_format, c("tidy", "matrix"))
}


# Resolve the deprecated `varying` argument
resolve_group_effects = function(group_effects, varying, group_effects_missing, what,
                                 env = rlang::caller_env(),
                                 user_env = rlang::caller_env(2)) {
  if (lifecycle::is_present(varying)) {
    lifecycle::deprecate_warn(
      "0.4.0",
      paste0(what, "(varying)"),
      paste0(what, "(group_effects)"),
      env = env,
      user_env = user_env
    )
    if (!group_effects_missing)
      stop("Use only one of `group_effects` and deprecated `varying`.")
    group_effects = varying
  }
  group_effects
}


# Resolve deprecated `probs` (from v0.3.4) to `width`. Only a single central interval can be translated.
resolve_width = function(width, probs, width_missing, what,
                         env = rlang::caller_env(),
                         user_env = rlang::caller_env(2)) {
  if (lifecycle::is_present(probs)) {
    lifecycle::deprecate_warn(
      "0.4.0",
      paste0(what, "(probs)"),
      paste0(what, "(width)"),
      details = "The interval columns are now `.lower` and `.upper`. For other quantiles, summarise the draws from `summary = FALSE`.",
      env = env,
      user_env = user_env
    )
    if (!width_missing)
      stop("Use only one of `width` and deprecated `probs`.")
    if (isTRUE(probs)) {
      width = 0.95
    } else if (isFALSE(probs)) {
      width = FALSE
    } else if (is.numeric(probs) && length(probs) == 2 && isTRUE(all.equal(sum(probs), 1))) {
      width = abs(diff(probs))
    } else {
      stop("Deprecated `probs` can only be translated to a central interval, e.g., `probs = c(0.025, 0.975)`. Use `width` instead, or summarise the draws from `summary = FALSE` for other quantiles.", call. = FALSE)
    }
  }
  width
}


# Warn when a binomial prediction relies on the pre-v0.4 proportion default.
warn_binomial_rate = function(fit, newdata, what) {
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


# Signal deprecated and removed mcpfit fields from v0.3.4
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
  stop_removed_mcpfit_field(name)
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
    stop_removed_mcpfit_field(i)
  .subset2(x, i)
}


# Stop with a pointer to the replacement when accessing mcpfit fields removed after v0.3.4
stop_removed_mcpfit_field = function(name) {
  removed = c(pars = "mcp_pars(fit) and mcp_columns(fit)", mcmc_loglik = "log_lik(fit)")
  if (name %in% names(removed))
    lifecycle::deprecate_stop("0.4.0", I(paste0("fit$", name)), I(removed[[name]]))
}


# Defunct internal helper for MCMC draw extraction, advertised as `mcp:::tidy_samples()` in v0.3.4
tidy_samples = function(...) {
  lifecycle::deprecate_stop("0.4.0", "tidy_samples()", I("as_draws_df() or tidybayes::spread_draws()"))
}
