
#' Posterior Predictive Checks For Mcpfit Objects
#'
#' Plot posterior (default) or prior (`prior = TRUE`) predictive checks. This is convenience wrapper
#' around the `bayesplot::ppc_*()` methods.
#'
#' @aliases pp_check pp_check.mcpfit
#' @inheritParams pp_eval
#' @param type One of `bayesplot::available_ppc("grouped", invert = TRUE) %>% stringr::str_remove("ppc_")`
#' @param facet_by Name of a column in data modeled as varying effect(s).
#' @param nsamples Number of draws. Note that you may want to use all data for summary geoms.
#'   e.g., `pp_check(fit, type = "ribbon", nsamples = NULL)`.
#' @param ... Further arguments passed to `bayesplot::ppc_type(y, yrep, ...)`
#' @return A `ggplot2` object for single plots. Enriched by `patchwork` for faceted plots.
#' @seealso \code{\link{plot.mcpfit}} \code{\link{pp_eval}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \donttest{
#' pp_check(demo_fit)
#' pp_check(demo_fit, type = "ecdf_overlay")
#' #pp_check(some_varying_fit, type = "loo_intervals", facet_by = "id")
#' }
pp_check = function(
  object,
  type = "dens_overlay",
  facet_by = NULL,
  newdata = NULL,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  nsamples = 100,
  ...
) {
  # Internal mcp naming convention
  fit = object
  assert_types(fit, "mcpfit")
  assert_types(facet_by, "null", "character", len = c(0, 1))
  assert_logical(prior)
  assert_types(varying, "logical", "character")
  assert_logical(arma)
  assert_integer(nsamples, lower = 1, len = 1)

  # Check and recode inputs
  if (!is.null(facet_by))
    if (!is.character(facet_by) || length(facet_by) != 1)
      stop("`facet_by` must be a single character string.")

  if (is.null(newdata)) {
    y = as.numeric(fit$data[, fit$pars$y])  # strip simulated data of attributes
    varying_data = fit$data[, facet_by]
  } else {
    assert_types(newdata, "data.frame", "tibble")
    y = as.numeric(newdata[, fit$pars$y])  # strip simulated data of attributes
    varying_data = newdata[, facet_by]
  }

  allowed_types = stringr::str_remove(bayesplot::available_ppc(), "ppc_")
  allowed_types = allowed_types[stringr::str_detect(allowed_types, "_grouped") == FALSE]  # Grouped done mcp-side (see below)
  if (type %notin% allowed_types)
    stop("`type` must be one of '", paste0(allowed_types, collapse = "', '"), "'")

  # Get as tidy samples to preserve info on groups and sampled draws
  samples = pp_eval(
    fit,
    newdata = newdata,
    summary = FALSE,
    type = "predict",
    probs = FALSE,
    rate = FALSE,
    prior = prior,
    which_y = "mu",
    varying = varying,
    arma = arma,
    nsamples = nsamples,
    samples_format = "tidy"
  )

  # Return plot with or without facets
  if (is.null(facet_by)) {
    yrep = tidy_to_matrix(samples, "predict")
    plot_return = get_ppc_plot(fit, type, y, yrep, nsamples)
    return(plot_return)
  } else {
    groups = unique(varying_data)
    all_plots = list()
    for (group in groups) {
      # Compute/extract y and yrep for this group
      y_this = y[varying_data == group]
      samples_this = dplyr::filter(samples, !!rlang::sym(facet_by) == group)
      yrep_this = tidy_to_matrix(samples_this, "predict")

      # Add plot to list
      all_plots[[group]] = get_ppc_plot(fit, type, y_this, yrep_this, nsamples, draws = samples_this$.draw) +
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
#' @param draws (required for loo-type plots) Indices of draws to use.
#' @param ... Arguments passed to `bayesplot::ppc_type(y, yrep, ...)`
#' @return A `ggplot2` object returned by `tidybayes::ppc_*(y, yrep, ...)`.
#' @return A string
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_ppc_plot = function(fit, type, y, yrep, nsamples, draws = NULL, ...) {
  is_loo = stringr::str_detect(type, "loo")

  func_name = paste0("ppc_", type)
  func_obj = utils::getFromNamespace(func_name, "bayesplot")

  if (is_loo == FALSE) {
    return(suppressWarnings(func_obj(y, yrep, ...)))
    #bayesplot_call = paste0("bayesplot::ppc_", type, "(y, yrep, ...)")
  } else if (is_loo == TRUE) {
    # Compute loo if missing
    fit = with_loo(fit, save_psis = TRUE, info = "Computing `fit$loo = loo(fit, save_psis = TRUE)`...")

    # Extract psis_object and lw
    psis_object = fit$loo$psis_object
    lw = psis_object$log_weights[unique(draws), ]
    psis_object$log_weights = lw
    attr(psis_object, "dims") = c(dim(yrep))

    # Build call (setting `samples` overwrites bayesplot defaults)
    return(suppressWarnings(func_obj(y, yrep, psis_object = psis_object, lw = lw, samples = nsamples, ...)))
  }
}
