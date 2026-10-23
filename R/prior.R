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
  legacy_parameter_names(names(prior), "prior")

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

# The prior table: resolved priors and metadata for all parameters in the model
get_prior_table = function(segments, cps, predictors, group_effects, family, prior = list(), data, design_specs = list(), reparam = TRUE) {
  checkmate::assert_class(family, "mcpfamily")
  context = prior_context(data, segments, design_specs)
  context$cp_names = cps$name
  context$cp_pars = stats::na.omit(c(cps$name, cps$group_name))
  legacy_prior_constants(prior, context)
  predictors$prior_type = get_prior_type(predictors)

  # Disallow cp_0 and cp_K (at least for now)
  prior_text = unlist(prior[vapply(prior, is.character, logical(1))])
  cp_limits = paste0("cp_", c(0, context$n_segments))
  if (any(stringr::str_detect(prior_text, paste0("\\b(", paste(cp_limits, collapse = "|"), ")\\b"))))
    stop(
      "Priors cannot refer to `", cp_limits[1], "` or `", cp_limits[2], "`, which are not parameters. ",
      "Use `min(", context$x_display, ")` and `max(", context$x_display, ")` for the ends of the observed range.",
      call. = FALSE
    )

  # Begin building actual prior
  specs = dplyr::bind_rows(
    default_cp_specs(cps, context),
    default_predictor_specs(predictors, family, context),
    default_group_specs(group_effects, family)
  )

  # Parse default and specified priors once. A constant response (e.g., a placeholder before
  # simulating data) has no scale. Use 1 instead, as rstanarm does for non-Gaussian families.
  if (!isTRUE(context$y_sd > 0))
    specs$code = gsub("sd\\(\\.y( - offset(_[0-9]+)?)?\\)", "1", specs$code)
  specs$parsed = lapply(instantiate_prior_template(specs$code, context), parse_prior)
  specs$code = NULL
  prior = lapply(prior, parse_prior)

  compile = function(prior) {
    specs = overlay_user_prior_specs(specs, prior, cps, context, predictors, family)
    specs = add_intercept_references(specs, predictors, context)
    compile_prior_specs(specs, specs$parameter, context)
  }
  prior_table = compile(prior)

  # Specified priors that resolve to the default, e.g., from `mcp(..., prior = fit$prior)`, 
  # are treated as defaults (not specified). This ensures default optimizations remain active.
  # Compare at 15 digits because R's parser can be 1 ulp off for 17-digit numbers on some platforms (macOS arm64).
  if (length(prior) > 0) {
    defaults = compile(list())
    value = function(prior_table) format_prior_display(prior_table$value[match(names(prior), prior_table$parameter)], digits = 15)
    is_default = value(prior_table) == value(defaults)
    if (any(is_default))
      prior_table = compile(prior[!is_default])
  }

  # How each parameter is sampled
  prior_table = add_nodes(prior_table, cps, group_effects$name)
  add_reparams(prior_table, predictors, cps, context, reparam)
}


#' Summarise priors used by an mcp model
#'
#' `r lifecycle::badge("experimental")`
#'
#' Shows the specified priors in the syntax of `mcp(..., prior)`, resolved to numbers.
#' As in `prior`, scales are SD/scale rather than JAGS precision (see [mcp-priors]). `bounds` are the lower and upper bounds of
#' each prior, i.e., the arguments of `dunif()` or `T()`, or the support of the distribution.
#' Change points are also "ordered": their specified priors are restricted to ordered change
#' points. Use `plot_pars(fit, prior = TRUE)` to see the marginal priors after this restriction.
#' `verbose = TRUE` also displays the underlying rule, description, source, and kind.
#' Numbers are rounded to 5 significant digits for display. `fit$prior` holds the exact
#' values used in the model. Parameter names within a prior describe its dependence on
#' other parameters, e.g., the default prior on `Intercept_1` given `x_1`.
#'
#' `prior_summary()` shares its generic with \pkg{rstantools} (used by `brms` and `rstanarm`), so it
#' works for all of these models regardless of package load order.
#'
#' @param object An `mcpfit` object.
#' @param verbose Logical. Include rule, description, source, and kind.
#' @param ... Must be empty. Reserved for future use.
#' @return A tibble with one row per model parameter, ordered and labeled
#'   the same way as `summary()`: change points first, then `mu`, then the
#'   other distributional parameters, then AR/MA lags - each with
#'   `segment` and `dpar` columns.
#' @seealso [mcp-priors] for how to specify priors, and
#'   [Default priors](https://lindeloev.github.io/mcp/articles/default_priors.html) for the rules behind them.
#' @export
#' @examples
#' prior_summary(demo_fit)  # Show the resolved priors and their bounds
#' prior_summary(demo_fit, verbose = TRUE)  # Include their rules and sources
prior_summary = function(object, ...) {
  # If mcp masks rstantools::prior_summary(), forward e.g. brmsfit objects to it
  if (!inherits(object, "mcpfit") && isNamespaceLoaded("rstantools"))
    return(rstantools::prior_summary(object, ...))
  UseMethod("prior_summary")
}


#' @rdname prior_summary
#' @export
#' @exportS3Method rstantools::prior_summary
prior_summary.mcpfit = function(object, verbose = FALSE, ...) {
  rlang::check_dots_empty()
  fit = object  # Standard name in mcp
  checkmate::assert_flag(verbose)
  prior_table = fit$.internal$prior_table

  pars = mcp_pars(fit)
  match_idx = match(prior_table$parameter, pars$parameter)
  prior_table$segment = pars$segment[match_idx]
  prior_table$dpar = pars$dpar[match_idx]
  prior_table = prior_table[order(match_idx), ]

  public = c("parameter", "segment", "dpar", "prior", "bounds")
  if (verbose)
    public = c(public, "rule", "description", "source", "kind")
  prior_table[, public, drop = FALSE]
}
