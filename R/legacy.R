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


# Warn on deprecated which_y argument
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
