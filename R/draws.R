# ABOUT: Extracting, transforming, and indexing posterior and prior draws.
# -----------------------------------------------------------------------

# Select posterior or explicitly requested prior draws.
mcmclist_draws = function(fit, prior = FALSE, error = TRUE) {
  check_mcpfit_version(fit)
  checkmate::assert_flag(prior)
  draws = .subset2(fit, if (prior) "mcmc_prior" else "mcmc_post")
  if (coda::is.mcmc.list(draws))
    return(draws)
  if (error)
    stop(if (prior) "No prior draws are available." else
      "No posterior draws are available. Select prior draws explicitly where supported, using `prior = TRUE`.",
      call. = FALSE)
  NULL
}


# Convert the selected draws to a posterior draws array.
posterior_draws = function(fit, prior = FALSE, error = TRUE) {
  draws = mcmclist_draws(fit, prior = prior, error = error)
  if (is.null(draws))
    return(NULL)
  posterior::as_draws_array(draws)
}


#' Extract MCMC Draws from `mcpfit` Objects
#'
#' Extract posterior or prior draws using \pkg{posterior}, \pkg{coda}, or the optional \pkg{tidybayes} package's S3 generics.
#'
#' @aliases as_draws as_draws.mcpfit as_draws_df.mcpfit as_draws_array.mcpfit as_draws_matrix.mcpfit as_draws_rvars.mcpfit as.mcmc.mcpfit tidy_draws.mcpfit
#' @param x An \code{\link{mcpfit}} object.
#' @param prior Logical. Extract prior draws (`TRUE`) instead of posterior draws
#'   (`FALSE`)? Errors if the requested draws are unavailable.
#' @param ... Passed to \pkg{posterior} or \pkg{tidybayes} format conversion functions.
#' @return A \pkg{posterior} `draws` object or a \pkg{coda} `mcmc.list` object.
#' @examples
#' # Default posterior draws, with one row per iteration and chain
#' draws = as_draws(demo_fit)  # Return a posterior::draws object
#' head(as_draws_df(demo_fit))  # Convert draws to a data frame
#'
#' # Other posterior formats are useful in different downstream packages
#' as_draws_matrix(demo_fit)[1:3, 1:3]  # Matrix of draws by parameter
#' as_draws_array(demo_fit)[1:2, , 1:2]  # Iteration-by-chain-by-parameter array
#' as_draws_rvars(demo_fit)[c("cp_1", "cp_2")]  # Random-variable representation
#'
#' # mcp also supports the coda and tidybayes conventions
#' head(coda::as.mcmc(demo_fit)[[1]])  # First chain as a coda mcmc object
#' if (requireNamespace("tidybayes", quietly = TRUE))
#'   head(tidybayes::tidy_draws(demo_fit))
#' @exportS3Method posterior::as_draws
as_draws.mcpfit = function(x, prior = FALSE, ...) {
  posterior_draws(x, prior = prior)
}

#' @exportS3Method posterior::as_draws_df
as_draws_df.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_df(posterior_draws(x, prior = prior), ...)
}

#' @exportS3Method posterior::as_draws_array
as_draws_array.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_array(posterior_draws(x, prior = prior), ...)
}

#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_matrix(posterior_draws(x, prior = prior), ...)
}

#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_rvars(posterior_draws(x, prior = prior), ...)
}

#' @rdname as_draws.mcpfit
#' @export
as_draws = posterior::as_draws

#' @rdname as_draws.mcpfit
#' @export
as_draws_df = posterior::as_draws_df

#' @rdname as_draws.mcpfit
#' @export
as_draws_array = posterior::as_draws_array

#' @rdname as_draws.mcpfit
#' @export
as_draws_matrix = posterior::as_draws_matrix

#' @rdname as_draws.mcpfit
#' @export
as_draws_rvars = posterior::as_draws_rvars

#' @exportS3Method coda::as.mcmc
as.mcmc.mcpfit = function(x, prior = FALSE, ...) {
  mcmclist_draws(x, prior = prior)
}

#' @exportS3Method tidybayes::tidy_draws
tidy_draws.mcpfit = function(model, ...) {
  posterior::as_draws_df(model, ...)
}


#' Index \code{mcpfit} objects
#'
#' Index variables, iterations, chains, and draws.
#'
#' @inheritParams fitted.mcpfit
#' @param x An `mcpfit` object or a posterior draws object.
#' @details These methods require posterior draws. To count prior draws, use
#'   e.g. `ndraws(as_draws(fit, prior = TRUE))`.
#' @return An integer count of iterations, chains, or draws.
#' @name draws-index-mcp
#' @examples
#' niterations(demo_fit)
#' nchains(as_draws(demo_fit, prior = TRUE))
NULL


#' @aliases niterations.mcpfit
#' @describeIn draws-index-mcp Number of iterations per chain of an `mcpfit` object.
#' @exportS3Method posterior::niterations
niterations.mcpfit = function(x, ...) {
  coda::niter(mcmclist_draws(x))
}

#' @aliases nchains.mcpfit
#' @describeIn draws-index-mcp Number of chains of an `mcpfit` object.
#' @exportS3Method posterior::nchains
nchains.mcpfit = function(x, ...) {
  coda::nchain(mcmclist_draws(x))
}

#' @rdname draws-index-mcp
#' @exportS3Method posterior::ndraws
ndraws.mcpfit = function(x, ...) {
  draws = mcmclist_draws(x)
  sum(vapply(draws, nrow, integer(1)))
}

#' @rdname draws-index-mcp
#' @export
ndraws = posterior::ndraws

#' @rdname draws-index-mcp
#' @export
nchains = posterior::nchains

#' @rdname draws-index-mcp
#' @export
niterations = posterior::niterations


# Resolve the deprecated `nsamples` argument
resolve_ndraws = function(ndraws, nsamples, ndraws_missing, what,
                         samples = lifecycle::deprecated(),
                         env = rlang::caller_env(),
                         user_env = rlang::caller_env(2)) {
  if (lifecycle::is_present(samples)) {
    lifecycle::deprecate_soft(
      "0.4.0",
      paste0(what, "(samples)"),
      env = env,
      user_env = user_env
    )
  }
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


#' Get tidy draws with or without group-level effects
#'
#' Extract posterior or prior draws formatted as tidy data frames
#'
#' Returns in a format useful for `fit$simulate()` with population-level parameters in wide format
#' and group-level deviations in long format (the number of rows is multiplied
#' by the number of selected group levels).
#'
#' @aliases mcp_draws
#' @keywords internal
#' @noRd
#' @inheritParams mcmclist_draws
#' @inheritParams pp_eval
#' @param population
#'   * `TRUE` All population-level model parameters.
#'   * `FALSE` No population-level effects. Same as `c()`.
#'   * Character vector: Only include specified population-level parameters.
#' @param varying Group-level effects. One of:
#'   * `TRUE` All group-level deviations.
#'   * `FALSE` No group-level deviations (`c()`).
#'   * `"cp"` or `"predictor"`: All group-level deviations belonging to that part of
#'     the model.
#'   * Character vector: Only include specified group-level parameters.
#' @param absolute
#'   * `TRUE` Returns the absolute location of all group-specific change points.
#'   * `FALSE` Return the group-level deviations.
#'   * Character vector: Apply the absolute transform only to these group-level parameters.
#'
#' @return `tibble` of posterior draws in `tidybayes` format.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
mcp_draws = function(
  fit,
  population = TRUE,
  varying = TRUE,
  absolute = FALSE,
  prior = FALSE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated()
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "mcp_draws")

  # General argument checks
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_multi_class(population, c("logical", "character"))
  checkmate::assert_multi_class(varying, c("logical", "character"), null.ok = TRUE)
  checkmate::assert_multi_class(absolute, c("logical", "character"), null.ok = TRUE)
  checkmate::assert_flag(prior)
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)

  if (all(population == FALSE) && all(varying == FALSE))
    stop("At least one TRUE or one parameter must be provided through either the `varying` or the `population` arguments.")


  # ----- IDENTIFY PARAMETERS -----
  # Group-level parameters.
  group_info = unpack_group_effects(fit, pars = varying)

  # Population-level parameters. Result is `pars_population`.
  if (all(population == FALSE)) {
    pars_population = c()  # Empty if no absolute group-level change points
  } else if (all(population == TRUE)) {
    pars_population = mcp_pars(fit, scope = "population")$name
  } else if (is.character(population)) {
    if (!all(population %in% mcp_pars(fit, scope = "population")$name))
      stop("Not all `population` selections are population-level parameters.")

    pars_population = population
  }

  # Absolute effects. Results are `absolute_cps` and `absolute`.
  if (all(absolute == TRUE)) {
    cp_effects = group_info$effects[group_info$effects$part == "cp", , drop = FALSE]
    absolute = cp_effects$name
    absolute_cps = cp_effects$population_name
  } else if (all(absolute == FALSE)) {
    absolute_cps = NULL
  } else if (is.character(absolute)) {
    # Check
    is_group_selection = absolute %in% group_info$pars
    if (any(!is_group_selection))
      stop("The following parameter names in `absolute` are not in `varying`: ", and_collapse(absolute[!is_group_selection]))
    absolute_effects = group_info$effects[
      group_info$effects$name %in% absolute, , drop = FALSE
    ]
    if (any(absolute_effects$part != "cp"))
      stop("`absolute` can select change-point group-level effects only.")

    absolute_cps = absolute_effects$population_name
  }

  # ----- GET THESE PARAMETERS AS TIDY DRAWS -----
  # Select draws before expanding across group levels.
  draws = tibble::as_tibble(posterior::as_draws_df(fit, prior = prior))
  if (!is.null(ndraws))
    draws = dplyr::sample_n(draws, ndraws)

  # Prepare for tidyr::pivot_longer_spec:
  # Describe the reshape using original factor levels
  groups = split(group_info$effects$name, group_info$effects$group_col)
  specs = lapply(names(groups), function(col) {
    spec = tidyr::expand_grid(.value = groups[[col]], !!col := unique(fit$data[[col]]))
    spec$.name = paste0(spec$.value, "[", spec[[col]], "]")
    spec
  })
  group_nodes = unlist(lapply(specs, function(spec) spec$.name))
  draws = dplyr::select(draws, dplyr::all_of(unique(c(
    ".chain", ".iteration", ".draw", pars_population, absolute_cps, group_nodes
  ))))

  # Pivot shared group effects together; successive pivots cross grouping columns.
  for (spec in specs)
    draws = tidyr::pivot_longer_spec(draws, spec)

  # Add population-level change points to deviations, then remove helper columns.
  if (length(absolute_cps) > 0) {
    draws[, absolute_cps] = draws[, absolute_cps] + draws[, absolute]
    draws = dplyr::select(draws, -dplyr::all_of(absolute))
  }

  # Unassigned group-level deviations are simulated as zero (the population-level mean).
  remaining_group_cols = dplyr::setdiff(mcp_pars(fit, scope = "group")$name, colnames(draws))
  draws[, remaining_group_cols] = 0

  # Return with chain etc. first
  draws %>%
    dplyr::relocate(".chain", ".iteration", ".draw")
}


# Internal helper to extract retained posterior imputations of missing responses.
# These are draws of the missing responses conditional on the observed data,
# rather than fresh outcomes from predict(). For AR/MA models, later
# observations can inform an imputation. Prediction uses these values only
# to complete the history for subsequent observations.
imputed_draws = function(object, ndraws = NULL) {
  mcmclist_draws(object)  # Validate the fit and availability of posterior draws.
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)
  imputed = object$.internal$imputed_response
  if (is.null(imputed))
    stop("No retained posterior imputations are available in this fit.", call. = FALSE)

  draws = tibble::as_tibble(posterior::as_draws_df(posterior::as_draws_array(imputed)))
  if (!is.null(ndraws))
    draws = dplyr::sample_n(draws, ndraws)
  rows = object$.internal$imputed_response_rows
  spec = tibble::tibble(
    .name = paste0(mcp_columns(object)$response, "[", rows, "]"),
    .value = ".imputed",
    data_row = rows
  )
  tidyr::pivot_longer_spec(draws, spec) %>%
    dplyr::relocate(".chain", ".iteration", ".draw")
}


# Deprecated internal helper for MCMC draw extraction
tidy_samples = function(...) {
  lifecycle::deprecate_soft(
    when = "0.4.0",
    what = "tidy_samples()",
    with = "as_draws_df() or tidybayes::spread_draws()"
  )
  mcp_draws(...)
}
