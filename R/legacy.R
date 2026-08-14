# Compatibility and Deprecation Detections for Legacy Code & Fits -------------

#' Check if an mcpfit object is compatible with the current version of mcp
#'
#' Throws an informative deprecation error if an mcpfit object was saved using
#' mcp < 0.4.0, which lacks modern internal data structures and parameter names.
#'
#' @keywords internal
#' @noRd
#' @param fit An `mcpfit` object.
check_mcpfit_version = function(fit) {
  checkmate::assert_class(fit, "mcpfit")

  tables = fit$.internal$model_tables
  is_valid = !is.null(fit$.internal) &&
    is.mcpfamily(fit$family) &&
    !is.null(tables) &&
    !is.null(tables$pars) &&
    !is.null(tables$design_specs) &&
    "role" %in% names(tables$pars)

  if (!is_valid) {
    stop(
      "This `mcpfit` object was created with an older version of mcp (< 0.4.0).\n",
      "Internal data structures and parameter names changed in mcp v0.4.0.\n",
      "Please re-fit the model using mcp >= 0.4.0.",
      call. = FALSE
    )
  }

  TRUE
}


#' Warn on deprecated `which_y` argument
#'
#' @keywords internal
#' @noRd
#' @param args List of dots arguments passed to a function.
#' @param func_name Name of the calling function for context.
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
