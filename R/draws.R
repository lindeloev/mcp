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
#' Base R's `as.matrix()`, `as.array()`, and `as.data.frame()` return plain draws, like in \pkg{rstanarm} and \pkg{brms},
#' so \pkg{bayesplot} functions work directly on the fit.
#'
#' @aliases as_draws as_draws.mcpfit as_draws_df.mcpfit as_draws_array.mcpfit as_draws_matrix.mcpfit as_draws_rvars.mcpfit as.mcmc.mcpfit tidy_draws.mcpfit
#' @param x An \code{\link{mcpfit}} object.
#' @param prior Logical. Extract prior draws (`TRUE`) instead of posterior draws
#'   (`FALSE`)? Errors if the requested draws are unavailable.
#' @param ... Passed to \pkg{posterior} or \pkg{tidybayes} format conversion functions.
#' @return A \pkg{posterior} `draws` object or a \pkg{coda} `mcmc.list` object.
#'   `as.matrix()` returns a draws-by-parameter matrix, `as.array()` an iteration-by-chain-by-parameter
#'   array, and `as.data.frame()` a draws-by-parameter data frame. The data frame has no chain information,
#'   so use `as_draws_df()` for chain-aware diagnostics such as `posterior::rhat()`.
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
#' # Base R formats, e.g., for bayesplot
#' as.matrix(demo_fit)[1:3, 1:3]  # Plain matrix of draws by parameter
#' bayesplot::mcmc_trace(demo_fit, pars = c("cp_1", "cp_2"))
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

# Base R coercions return plain draws, like rstanarm and brms. This also makes bayesplot::mcmc_*(fit) work.
#' @rdname as_draws.mcpfit
#' @export
as.matrix.mcpfit = function(x, prior = FALSE, ...) {
  unclass(as_draws_matrix.mcpfit(x, prior = prior, ...))  # Keeps attr "nchains" for posterior::rhat() etc.
}

#' @rdname as_draws.mcpfit
#' @export
as.array.mcpfit = function(x, prior = FALSE, ...) {
  unclass(as_draws_array.mcpfit(x, prior = prior, ...))
}

#' @rdname as_draws.mcpfit
#' @param row.names,optional Unused. For compatibility with [as.data.frame()].
#' @export
as.data.frame.mcpfit = function(x, row.names = NULL, optional = FALSE, prior = FALSE, ...) {
  as.data.frame(as.matrix.mcpfit(x, prior = prior, ...))
}

#' @exportS3Method tidybayes::tidy_draws
tidy_draws.mcpfit = function(model, ...) {
  # A plain tibble, like tidybayes' own method, so spread_draws() doesn't warn about dropping draws_df metadata
  tibble::as_tibble(posterior::as_draws_df(model, ...))
}


#' Index \code{mcpfit} objects
#'
#' Index variables, iterations, chains, and draws.
#'
#' @inheritParams fitted.mcpfit
#' @param x An `mcpfit` object or a posterior draws object.
#' @details These methods require posterior draws. To index prior draws, use
#'   e.g. `ndraws(as_draws(fit, prior = TRUE))`.
#' @return `variables()` returns the names of the draws' variables, e.g., `"cp_1"` and
#'   `"cp_1_id[John]"`. The other methods return an integer count of iterations, chains, or draws.
#' @name draws-index-mcp
#' @examples
#' variables(demo_fit)
#' niterations(demo_fit)
#' nchains(as_draws(demo_fit, prior = TRUE))
NULL


#' @aliases variables.mcpfit
#' @describeIn draws-index-mcp Variable names of an `mcpfit` object, as in `as_draws(fit)`.
#' @exportS3Method posterior::variables
variables.mcpfit = function(x, ...) {
  posterior::variables(posterior_draws(x))
}

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
variables = posterior::variables

#' @rdname draws-index-mcp
#' @export
ndraws = posterior::ndraws

#' @rdname draws-index-mcp
#' @export
nchains = posterior::nchains

#' @rdname draws-index-mcp
#' @export
niterations = posterior::niterations


#' Get draws with group-level deviations for selected group levels
#'
#' One row per draw with population-level parameters and group-level deviations as
#' `name[level]` columns, to be matched to data rows using `draws_to_pars()`. Group-level
#' deviations that are not selected are zero (the population-level mean).
#'
#' @inheritParams mcmclist_draws
#' @inheritParams fitted.mcpfit
#' @param group Group-level deviations. One of:
#'   * `TRUE` All group-level deviations.
#'   * `FALSE` No group-level deviations (`c()`).
#'   * `"cp"` or `"predictor"`: All group-level deviations belonging to that part of
#'     the model.
#'   * Character vector: Only include specified group-level parameters (`*` is a wildcard).
#' @param newdata Optional prediction data selecting group levels. Levels must
#'   be in the fitted data.
#' @return `tibble` of draws.
#' @noRd
mcp_draws = function(fit, group = TRUE, prior = FALSE, ndraws = NULL, newdata = NULL) {
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_flag(prior)
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)
  groups = unpack_group_effects(fit, pars = group)$groups
  if (!is.null(newdata))
    assert_data_cols(newdata, names(groups))

  draws = tibble::as_tibble(posterior::as_draws_df(fit, prior = prior))
  if (!is.null(ndraws))
    draws = dplyr::arrange(dplyr::sample_n(draws, ndraws), .data$.draw)

  # Columns of the selected group-level deviations at the selected levels
  group_nodes = unlist(lapply(names(groups), function(group_col) {
    group_levels = if (is.null(newdata)) unique(fit$data[[group_col]]) else unique(newdata[[group_col]])
    novel = setdiff(group_levels, fit$data[[group_col]])
    if (length(novel) > 0)
      stop(
        "`newdata` has level(s) of the grouping factor '", group_col, "' that are not in the fitted data: ",
        and_collapse(paste0("'", novel, "'")), ". Use `group_effects = FALSE` (similar to `re_formula = NA` in `posterior_*()`) ",
        "to predict without group-level effects, or `fit$simulate()` to simulate new groups.",
        call. = FALSE
      )
    paste0(rep(groups[[group_col]], each = length(group_levels)), "[", group_levels, "]")
  }))
  draws = dplyr::select(draws, dplyr::all_of(c(".chain", ".iteration", ".draw", mcp_pars(fit, scope = "population")$parameter, group_nodes)))
  draws[, setdiff(mcp_pars(fit, scope = "group")$parameter, unlist(groups))] = 0
  draws
}


# Parameter values of draws from `mcp_draws()`, as used by `evaluate_formula()`:
# one value per draw for population-level parameters and one value per row of `data` and draw
# (a rows x draws matrix) for group-level deviations, matched by the grouping columns in `data`.
# - groups: Group-level deviation names by grouping column, e.g., `list(id = c("Intercept_1_id", "x_1_id"))`.
draws_to_pars = function(draws, data, groups) {
  is_par = !grepl("[", names(draws), fixed = TRUE) & names(draws) %notin% c(".chain", ".iteration", ".draw")
  pars = as.list(draws[is_par])
  for (group_col in names(groups)) {
    group_levels = unique(data[[group_col]])
    group_indices = match(data[[group_col]], group_levels)
    for (name in groups[[group_col]])
      pars[[name]] = t(unname(as.matrix(draws[paste0(name, "[", group_levels, "]")])))[group_indices, , drop = FALSE]
  }
  pars
}


# Internal helper to extract retained posterior imputations of missing responses.
# These are draws of the missing responses conditional on the observed data,
# rather than fresh outcomes from predict(). For AR/MA models, later
# observations can inform an imputation. Prediction uses these values only
# to complete the history for subsequent observations.
imputed_draws = function(object, ndraws = NULL) {
  mcmclist_draws(object)  # Validate the fit and availability of posterior draws.
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)
  imputed_response = object$.internal$imputed_response
  if (is.null(imputed_response))
    stop("No retained posterior imputations are available in this fit.", call. = FALSE)

  draws = tibble::as_tibble(posterior::as_draws_df(posterior::as_draws_array(imputed_response)))
  if (!is.null(ndraws))
    draws = dplyr::sample_n(draws, ndraws)
  imputed_response_rows = object$.internal$imputed_response_rows
  pivot_spec = tibble::tibble(
    .name = paste0(mcp_columns(object)$response, "[", imputed_response_rows, "]"),
    .value = ".imputed",
    .row = imputed_response_rows
  )
  tidyr::pivot_longer_spec(draws, pivot_spec) |>
    dplyr::relocate(".chain", ".iteration", ".draw")
}
