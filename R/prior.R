# Prior assembly -----------------------------------------------------------

#' Get priors for all parameters in the model
#'
#' @keywords internal
#' @noRd
get_prior = function(segments, cps, predictors, group_effects, family, prior = list(), data) {
  checkmate::assert_true(is.mcpfamily(family), .var.name = "family")
  context = prior_context(data, segments)
  warn_legacy_prior_constants(prior, context)

  specs = dplyr::bind_rows(
    default_cp_specs(cps, context),
    default_predictor_specs(predictors, family),
    default_group_specs(group_effects, family)
  )
  specs = overlay_user_prior_specs(specs, prior, cps, context)

  all_names = specs$parameter
  table = compile_prior_specs(specs, all_names, context)
  resolved = stats::setNames(as.list(table$value), table$parameter)
  attr(resolved, "prior_table") = table
  attr(resolved, "prior_context") = context[c(
    "x_name", "y_name", "x_display", "y_display", "x_min", "x_max",
    "x_span", "n_cp", "n_segments", "segment_width"
  )]
  resolved
}


#' Summarise priors used by an mcp model
#'
#' Shows the effective, resolved priors on the familiar SD/scale
#' parameterization rather than JAGS precision. Use `verbose = TRUE` to also
#' see the symbolic rule, its description, source, and kind.
#'
#' @param fit An `mcpfit` object.
#' @param verbose Logical. Include rule, description, source, and kind.
#' @return A tibble with one row per model parameter, ordered and labeled
#'   the same way as `summary()`: change points first, then `mu`, then the
#'   other distributional parameters, then `ar`/`ma` components - each with
#'   `segment` and `dpar` columns.
#' @export
prior_summary = function(fit, verbose = FALSE) {
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_flag(verbose)
  table = fit$.internal$prior_table
  if (is.null(table))
    table = attr(fit$prior, "prior_table")
  if (is.null(table))
    stop("This `mcpfit` object does not contain a resolved prior table. Please re-fit the model using mcp >= 0.4.0.")

  pars = mcp_pars(fit)
  if (!is.null(pars)) {
    match_idx = match(table$parameter, pars$name)
    table$segment = pars$segment[match_idx]
    table$dpar = pars$dpar[match_idx]
    table = table[order(match_idx), ]
  }

  public = c("parameter", "segment", "dpar", "prior", "bounds")
  if (verbose)
    public = c(public, "rule", "description", "source", "kind")
  public = intersect(public, names(table))
  table[, public, drop = FALSE]
}
