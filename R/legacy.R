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
