# Prior assembly -----------------------------------------------------------

# Validate the plain named-list prior format introduced before sampler-
# agnostic prior objects. A future front end can dispatch classed prior objects
# separately and retain this function as the legacy JAGS-string path.
validate_prior_v1 = function(prior) {
  checkmate::assert_list(prior)
  if (length(prior) > 0 &&
      (is.null(names(prior)) || anyNA(names(prior)) || any(!nzchar(names(prior))))) {
    stop("`prior` must be a completely named list; every entry needs a nonempty parameter name.")
  }
  check_legacy_parameter_names(names(prior), "prior")

  which_duplicated = duplicated(names(prior))
  if (any(which_duplicated))
    stop("`prior` has duplicated entries for the same parameter: ", and_collapse(names(prior)[which_duplicated]))

  valid_value = vapply(prior, function(value) {
    is_numeric = is.numeric(value) && length(value) == 1 && is.finite(value)
    is_string = is.character(value) && length(value) == 1 &&
      !is.na(value) && nzchar(trimws(value))
    is_numeric || is_string
  }, logical(1))
  if (any(!valid_value)) {
    stop(
      "Each `prior` value must be one finite numeric scalar or one nonempty character string. Invalid value(s): ",
      and_collapse(names(prior)[!valid_value])
    )
  }

  invisible(prior)
}

# Get priors for all parameters in the model
get_prior = function(segments, cps, predictors, group_effects, family, prior = list(), data) {
  checkmate::assert_true(is.mcpfamily(family), .var.name = "family")
  context = prior_context(data, segments)
  warn_legacy_prior_constants(prior, context)

  specs = dplyr::bind_rows(
    default_cp_specs(cps, context),
    default_predictor_specs(predictors, family),
    default_group_specs(group_effects, family)
  )
  specs = overlay_user_prior_specs(specs, prior, cps, context, predictors, family)

  all_names = specs$parameter
  table = compile_prior_specs(specs, all_names, context)
  resolved = stats::setNames(as.list(table$value), table$parameter)
  attr(resolved, "prior_table") = table
  attr(resolved, "prior_context") = context[c(
    "x_name", "y_name", "x_display", "y_display", "x_min", "x_max",
    "x_span", "n_cp", "n_segments"
  )]
  resolved
}


#' Summarise priors used by an mcp model
#'
#' \lifecycle{experimental}
#'
#' Shows the effective, resolved prior distributions on the familiar SD/scale
#' parameterization rather than JAGS precision. Symbolic expressions in bounds
#' (e.g. `min(x)` or `max(x)`) may be retained in symbolic form in compact output,
#' while `verbose = TRUE` displays the underlying rule, description, source, and kind.
#'
#' @param fit An `mcpfit` object.
#' @param verbose Logical. Include rule, description, source, and kind.
#' @return A tibble with one row per model parameter, ordered and labeled
#'   the same way as `summary()`: change points first, then `mu`, then the
#'   other distributional parameters, then `ar`/`ma` components - each with
#'   `segment` and `dpar` columns.
#' @export
#' @examples
#' prior_summary(demo_fit)  # Show the effective priors and bounds
#' prior_summary(demo_fit, verbose = TRUE)  # Include their rules and sources
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
