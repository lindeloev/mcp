
#' Posterior Predictive Checks For Mcpfit Objects
#'
#' Plot posterior (default) or prior (`prior = TRUE`) predictive checks. This is a convenience wrapper
#' around the `bayesplot::ppc_*()` methods. See
#' [prior predictive checks](https://lindeloev.github.io/mcp/articles/priors.html#prior-predictive-checks) for an example.
#'
#' @rdname pp_check
#' @aliases pp_check pp_check.mcpfit
#' @inheritParams fitted.mcpfit
#' @param type String specifying the type of plot. See \code{\link[bayesplot:PPC-overview]{bayesplot::available_ppc()}}
#'   for available plot types, omitting the \code{"ppc_"} prefix (e.g., \code{"dens_overlay"},
#'   \code{"ribbon"}, \code{"intervals"}).
#' @param facet_by Name of a grouping column used by group-level effects.
#' @param varying Deprecated. Use `group_effects` instead.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @param ... Further arguments passed to `bayesplot::ppc_type(y, yrep, ...)`
#' @details Missing responses are omitted from the observed-data check. For
#'   GARMA models, each replicated response series is generated recursively
#'   rather than conditioning on the observed response history. LOO predictive
#'   checks use posterior draws and the original fitted data, so they do not
#'   support `prior = TRUE` or `newdata`.
#' @return A `ggplot2` object for single plots. Enriched by `patchwork` for faceted plots.
#' @seealso \code{\link{plot.mcpfit}} \code{\link{predict.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export pp_check
#' @exportS3Method bayesplot::pp_check
#' @seealso plot_pars plot_dpar pp_check
#' @examples
#' \donttest{
#' pp_check(demo_fit)
#' pp_check(demo_fit, type = "ecdf_overlay")
#' #pp_check(some_group_fit, type = "loo_intervals", facet_by = "id")
#' }
pp_check.mcpfit = function(
  object,
  type = "dens_overlay",
  facet_by = NULL,
  newdata = NULL,
  prior = FALSE,
  group_effects = TRUE,
  arma = TRUE,
  ndraws = 100,
  nsamples = lifecycle::deprecated(),
  varying = lifecycle::deprecated(),
  ...
) {
  ndraws = legacy_arg(ndraws, nsamples, missing(ndraws), "pp_check", "ndraws", "nsamples")
  group_effects = legacy_arg(group_effects, varying, missing(group_effects), "pp_check", "group_effects", "varying", lifecycle::deprecate_warn)

  # Internal mcp naming convention
  fit = object
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_string(facet_by, null.ok = TRUE)
  checkmate::assert_flag(prior)
  checkmate::assert_multi_class(group_effects, c("logical", "character"))
  checkmate::assert_flag(arma)
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)

  if (is.null(newdata)) {
    eval_data = fit$data
  } else {
    checkmate::assert_data_frame(newdata)
    eval_data = data.frame(newdata)
  }
  y_column = mcp_columns(fit)$response
  assert_data_cols(eval_data, y_column)
  if (!is.null(facet_by))
    assert_data_cols(eval_data, facet_by)

  y_all = as.numeric(eval_data[, y_column])  # strip simulated data of attributes
  observed_rows = which(!is.na(y_all))
  if (length(observed_rows) == 0)
    stop("`pp_check()` requires at least one observed response.")
  facet_data = if (is.null(facet_by)) NULL else eval_data[, facet_by]

  allowed_types = stringr::str_remove(bayesplot::available_ppc(), "ppc_")
  allowed_types = allowed_types[stringr::str_detect(allowed_types, "_grouped") == FALSE]  # Grouped done mcp-side (see below)
  if (type %notin% allowed_types)
    stop("`type` must be one of '", paste0(allowed_types, collapse = "', '"), "'")
  is_loo = stringr::str_detect(type, "loo")
  if (is_loo && prior)
    stop("LOO predictive checks require posterior draws; `prior = TRUE` is not supported.")
  if (is_loo && !is.null(newdata))
    stop("LOO predictive checks require the original fitted data; `newdata` is not supported.")
  if (is_loo && (!isTRUE(group_effects) || !arma))
    stop("LOO predictive checks evaluate the fitted model; `group_effects = TRUE` and `arma = TRUE` are required.")

  # Get tidy draws to preserve group and draw identifiers
  draws = pp_eval(
    fit,
    newdata = newdata,
    summary = FALSE,
    type = "predict",
    width = FALSE,
    rate = FALSE,
    prior = prior,
    dpar = NULL,
    group_effects = group_effects,
    arma = arma,
    # LOO weights are computed for every posterior draw. Keep those draws
    # intact; bayesplot's `samples` argument controls plot sampling where
    # supported.
    ndraws = if (is_loo) NULL else ndraws,
    draws_format = "tidy",
    .garma_replicate = TRUE
  )

  # LOO weights are computed once for all observed rows and subset per plot below.
  if (is_loo) {
    message("Computing `loo(fit, save_psis = TRUE)`...")
    psis_object = loo(fit, save_psis = TRUE)$psis_object
    if (nrow(psis_object$log_weights) != length(unique(draws$.draw)))
      stop_github("Could not align PPC draws with the LOO weights.")
  }
  get_psis_object = function(observations) {
    if (is_loo) subset_psis_object(psis_object, match(observations, observed_rows)) else NULL
  }

  # Return plot with or without facets
  if (is.null(facet_by)) {
    y = y_all[observed_rows]
    yrep = tidy_to_matrix(draws, type = ".prediction", data_rows = observed_rows)
    plot_return = get_ppc_plot(type, y, yrep, ndraws, get_psis_object(observed_rows), ...)
    return(plot_return)
  } else {
    facet_levels = unique(facet_data[observed_rows])
    all_plots = list()
    for (facet_level in facet_levels) {
      # Compute/extract y and yrep for this facet level
      observations_this = observed_rows[facet_data[observed_rows] == facet_level]
      y_this = y_all[observations_this]
      yrep_this = tidy_to_matrix(draws, type = ".prediction", data_rows = observations_this)

      # Add plot to list
      all_plots[[as.character(facet_level)]] = get_ppc_plot(
        type, y_this, yrep_this, ndraws, get_psis_object(observations_this), ...
      ) +
        ggplot2::ggtitle(facet_level)
    }

    # Return faceted plot using patchwork
    plot_return = patchwork::wrap_plots(all_plots) + patchwork::plot_layout(guides = "collect")
    return(plot_return)
  }
}


#' Call one bayesplot::ppc_*() function
#'
#' @inheritParams pp_check
#' @param y Response vector
#' @param yrep S X N matrix of predicted responses
#' @param psis_object A PSIS object for the columns of `yrep`. Required for loo-type plots.
#' @param ... Arguments passed to `bayesplot::ppc_type(y, yrep, ...)`
#' @return A `ggplot2` object returned by `bayesplot::ppc_*(y, yrep, ...)`.
#' @noRd
get_ppc_plot = function(type, y, yrep, ndraws, psis_object = NULL, ...) {
  func_name = paste0("ppc_", type)
  func_obj = utils::getFromNamespace(func_name, "bayesplot")
  plot_formals = names(formals(func_obj))
  checkmate::assert_subset(
    names(list(...)),
    setdiff(plot_formals, c("y", "yrep", "psis_object", "lw", "...")),
    .var.name = "names of arguments in `...`"
  )

  if (is.null(psis_object))
    return(suppressWarnings(rlang::exec(func_obj, y, yrep, ...)))

  plot_args = list(y = y, yrep = yrep)
  if ("psis_object" %in% plot_formals) {
    plot_args$psis_object = psis_object
  } else if ("lw" %in% plot_formals) {
    plot_args$lw = loo::weights.importance_sampling(psis_object)
  }
  if ("samples" %in% plot_formals && !is.null(ndraws) &&
      "samples" %notin% names(list(...)))
    plot_args$samples = ndraws

  suppressWarnings(rlang::exec(func_obj, !!!plot_args, ...))
}


# Subset a psis object by index
subset_psis_object = function(psis_object, observations) {
  psis_object$log_weights = psis_object$log_weights[
    , observations, drop = FALSE
  ]
  psis_object$diagnostics = lapply(
    psis_object$diagnostics,
    function(x) x[observations]
  )

  for (attribute in c("norm_const_log", "tail_len", "r_eff"))
    attr(psis_object, attribute) = attr(psis_object, attribute)[observations]
  attr(psis_object, "dims") = dim(psis_object$log_weights)
  psis_object
}
