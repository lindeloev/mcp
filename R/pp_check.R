
#' Posterior Predictive Checks For Mcpfit Objects
#'
#' Plot posterior (default) or prior (`prior = TRUE`) predictive checks. This is convenience wrapper
#' around the `bayesplot::ppc_*()` methods.
#'
#' @aliases pp_check pp_check.mcpfit
#' @param type String specifying the type of plot. See \code{\link[bayesplot:PPC-overview]{bayesplot::available_ppc()}}
#'   for available plot types, omitting the \code{"ppc_"} prefix (e.g., \code{"dens_overlay"},
#'   \code{"ribbon"}, \code{"intervals"}).
#' @param facet_by Name of a grouping column used by group-level effects.
#' @param ndraws Number of posterior draws. Note that you may want to use all
#'   draws for summary geoms, e.g., `pp_check(fit, type = "ribbon", ndraws = NULL)`.
#'   LOO checks always evaluate all posterior draws to preserve their PSIS
#'   weights; where supported, `ndraws` is passed to bayesplot to control the
#'   number of plotted draws.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @param ... Further arguments passed to `bayesplot::ppc_type(y, yrep, ...)`
#' @details Missing responses are omitted from the observed-data check. For
#'   GARMA models, each replicated response series is generated recursively
#'   rather than conditioning on the observed response history. LOO predictive
#'   checks use posterior draws and the original fitted data, so they do not
#'   support `prior = TRUE` or `newdata`.
#' @return A `ggplot2` object for single plots. Enriched by `patchwork` for faceted plots.
#' @seealso \code{\link{plot.mcpfit}} \code{\link{pp_eval}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @seealso plot_pars plot_dpar pp_check
#' @examples
#' \donttest{
#' pp_check(demo_fit)
#' pp_check(demo_fit, type = "ecdf_overlay")
#' #pp_check(some_group_fit, type = "loo_intervals", facet_by = "id")
#' }
pp_check = function(
  object,
  type = "dens_overlay",
  facet_by = NULL,
  newdata = NULL,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = 100,
  nsamples = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "pp_check")

  # Internal mcp naming convention
  fit = object
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_string(facet_by, null.ok = TRUE)
  checkmate::assert_flag(prior)
  checkmate::assert_multi_class(varying, c("logical", "character"))
  checkmate::assert_flag(arma)
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)

  # Check and recode inputs
  if (!is.null(facet_by))
    if (!is.character(facet_by) || length(facet_by) != 1)
      stop("`facet_by` must be a single character string.")

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
  group_data = if (is.null(facet_by)) NULL else eval_data[, facet_by]

  allowed_types = stringr::str_remove(bayesplot::available_ppc(), "ppc_")
  allowed_types = allowed_types[stringr::str_detect(allowed_types, "_grouped") == FALSE]  # Grouped done mcp-side (see below)
  if (type %notin% allowed_types)
    stop("`type` must be one of '", paste0(allowed_types, collapse = "', '"), "'")
  is_loo = stringr::str_detect(type, "loo")
  if (is_loo && prior)
    stop("LOO predictive checks require posterior draws; `prior = TRUE` is not supported.")
  if (is_loo && !is.null(newdata))
    stop("LOO predictive checks require the original fitted data; `newdata` is not supported.")

  # Get tidy draws to preserve group and draw identifiers
  draws = pp_eval(
    fit,
    newdata = newdata,
    summary = FALSE,
    type = "predict",
    probs = FALSE,
    rate = FALSE,
    prior = prior,
    dpar = NULL,
    varying = varying,
    arma = arma,
    # LOO weights are computed for every posterior draw. Keep those draws
    # intact; bayesplot's `samples` argument controls plot sampling where
    # supported.
    ndraws = if (is_loo) NULL else ndraws,
    draws_format = "tidy",
    .garma_replicate = TRUE
  )
  # Return plot with or without facets
  if (is.null(facet_by)) {
    y = y_all[observed_rows]
    yrep = tidy_to_matrix(draws, type = ".prediction", data_rows = observed_rows)
    plot_return = get_ppc_plot(
      fit, type, y, yrep, ndraws,
      observations = observed_rows,
      varying = varying,
      arma = arma,
      ...
    )
    return(plot_return)
  } else {
    groups = unique(group_data[observed_rows])
    all_plots = list()
    for (group in groups) {
      # Compute/extract y and yrep for this group
      observations_this = observed_rows[group_data[observed_rows] == group]
      y_this = y_all[observations_this]
      yrep_this = tidy_to_matrix(draws, type = ".prediction", data_rows = observations_this)

      # Add plot to list
      all_plots[[as.character(group)]] = get_ppc_plot(
        fit, type, y_this, yrep_this, ndraws,
        observations = observations_this,
        varying = varying,
        arma = arma,
        ...
      ) +
        ggplot2::ggtitle(group)
    }

    # Return faceted plot using patchwork
    plot_return = patchwork::wrap_plots(all_plots) + patchwork::plot_layout(guides = "collect")
    return(plot_return)
  }
}


#' pp_check for loo statistics
#'
#' @aliases get_loo_plot_call
#' @keywords internal
#' @noRd
#' @inheritParams pp_check
#' @param y Response vector
#' @param yrep S X N matrix of predicted responses
#' @param observations (required for loo-type plots) Original data-row indices.
#' @param ... Arguments passed to `bayesplot::ppc_type(y, yrep, ...)`
#' @return A `ggplot2` object returned by `bayesplot::ppc_*(y, yrep, ...)`.
#' @return A string
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_ppc_plot = function(fit, type, y, yrep, ndraws,
                        observations = NULL, varying = TRUE, arma = TRUE, ...) {
  is_loo = stringr::str_detect(type, "loo")

  func_name = paste0("ppc_", type)
  func_obj = utils::getFromNamespace(func_name, "bayesplot")
  plot_formals = names(formals(func_obj))
  checkmate::assert_subset(
    names(list(...)),
    setdiff(plot_formals, c("y", "yrep", "psis_object", "lw", "...")),
    .var.name = "names of arguments in `...`"
  )

  if (is_loo == FALSE) {
    return(suppressWarnings(rlang::exec(func_obj, y, yrep, ...)))
  } else if (is_loo == TRUE) {
    # Compute loo if missing
    fit = with_loo(
      fit,
      save_psis = TRUE,
      varying = varying,
      arma = arma,
      info = "Computing `fit$loo = loo(fit, save_psis = TRUE)`..."
    )

    # Align PSIS columns with the observations shown in this plot. Posterior
    # draws are deliberately left intact so the PSIS object remains valid.
    loo_settings = attr(fit$loo, "mcp_settings")
    loo_observed_rows = loo_settings$observed_rows
    observation_cols = match(observations, loo_observed_rows)
    if (nrow(fit$loo$psis_object$log_weights) != nrow(yrep) ||
        anyNA(observation_cols))
      stop_github("Could not align PPC observations with the LOO weights.")

    psis_object = fit$loo$psis_object
    psis_object = subset_psis_object(psis_object, observation_cols)
    lw = loo::weights.importance_sampling(psis_object)

    plot_args = list(y = y, yrep = yrep)
    if ("psis_object" %in% plot_formals) {
      plot_args$psis_object = psis_object
    } else if ("lw" %in% plot_formals) {
      plot_args$lw = lw
    }
    if ("samples" %in% plot_formals && !is.null(ndraws) &&
        "samples" %notin% names(list(...)))
      plot_args$samples = ndraws

    return(suppressWarnings(rlang::exec(func_obj, !!!plot_args, ...)))
  }
}


#' Subset a PSIS object by observations
#'
#' @keywords internal
#' @noRd
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
