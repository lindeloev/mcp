#' Plot mcpfit
#'
#' plot(fit, "combo") calls bayesplot::mcmc_combo(fit$samples). Use it directly
#' for finer control.
#'
#' @aliases plot plot.mcpfit
#' @param x An mcpfit object
#' @param type String. One of "overlay" (default) or "combo".
#' @param draws Positive integer. Number of posterior draws to use when type = "overlay".
#' @param pars Vector of parameter names to plot or "population" (default) or
#'   "all" or "varying". Only relevant for type = "combo".
#' @param facet_by String. Name of a varying group. Only relevant for
#'   type = "overlay"
#' @param ... Currently ignored.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A \code{ggplot2} object.
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_point facet_wrap
#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
#' @importFrom stats sd
#' @export
#' @examples
#' \dontrun{
#' plot(fit)  # defaults to the below
#' plot(fit, type="overlay", draws=50) + ggtitle("Great fit!")
#' plot(fit, "combo")
#' }

plot.mcpfit = function(x, type="overlay", draws=25, pars="population", facet_by = NULL, ...) {
  # Check arguments
  if (class(x) != "mcpfit")
    stop("Can only plot mcpfit objects. x was class: ", class(x))
  if (! type %in% c("overlay", "combo"))
    stop("Type has to be one of 'overlay' or 'combo'. Was: ", type)
  if (draws < 1)
    stop("Draws has to be a positive integer.")
  if (!is.null(pars)) {
    if (pars != "population" & !all(pars %in% c(x$pars$population, x$pars$varying)))
      stop("Not all these pars are in the model: '", paste0(pars, collapse="' and '"))
  }

  # TEMPORARY: Include test of whether this is a random/nested effect
  varying_groups = logical0_to_null(unique(stats::na.omit(x$.other$ST$cp_group_col)))
  if (!is.null(facet_by)) {
    if (!facet_by %in% varying_groups)
      stop("facet_by is not a data column used as varying grouping.")
  }

  # Plot function on top
  if (type == "overlay") {

    # We'll need these vars during data-processing-before-ggplot
    func_y = x$func_y
    eval_at = seq(min(x$data[, x$pars$x]),
                  max(x$data[, x$pars$x]),
                  length.out = 100)
    pars_population_regex = paste0(x$pars$population, collapse="|")

    # No faceting
    if (is.null(facet_by)) {
      Q = x$samples %>%
        tidybayes::spread_draws(!!rlang::sym(pars_population_regex), regex = TRUE) %>%
        # TO DO: use spread_draws(n = draws) when tidybayes 1.2 is out
        tidybayes::sample_draws(draws)

    } else {
      # Prepare for faceting
      # Read more about this weird syntax at https://github.com/mjskay/tidybayes/issues/38
      varying_by_facet = stats::na.omit(x$.other$ST$cp_group[stringr::str_detect(x$.other$ST$cp_group, paste0("_", facet_by))])
      varying_by_facet = paste0(varying_by_facet, collapse="|")

      Q = x$samples %>%
        tidybayes::spread_draws(!!rlang::sym(pars_population_regex),
                     (!!rlang::sym(varying_by_facet))[!!rlang::sym(facet_by)],
                     regex = TRUE) %>%
        # TO DO: use spread_draws(n = draws) when tidybayes 1.2 is out
        tidybayes::sample_draws(draws)
    }

    # First, let's get all the predictors in shape for func_y
    Q = Q %>%
      tidyr::expand_grid(!!x$pars$x := eval_at) %>%  # correct name of x-var

      # Add fitted draws (vectorized)
      dplyr::mutate(!!x$pars$y := purrr::invoke(func_y, ., type = "fitted")) %>%

      # Add line ID to separate lines. Mark a new line when "eval_at" repeats.
      dplyr::mutate(
        line = !!dplyr::sym(x$pars$x) == min(eval_at),
        line = cumsum(line)
      )

    if (is.null(facet_by)) {
      # Return plot without faceting
      return(ggplot(x$data, aes_string(x = x$pars$x, y = x$pars$y)) +
        geom_point() +
        geom_line(aes(group = line), data = Q, color = grDevices::rgb(0.5, 0.5, 0.5, 0.4)))
    } else {
      # Return plot with faceting
      return(ggplot(x$data, aes_string(x = x$pars$x, y = x$pars$y)) +
        geom_point() +
        geom_line(aes(group = line), data = Q, color = grDevices::rgb(0.5, 0.5, 0.5, 0.4)) +
        facet_wrap(paste0("~", facet_by)))
    }
  }

  if (type == "combo") {
    # Use population parameters by default
    if (pars == "population") {
      pars = x$pars$population
    }
    if (length(pars) > 1) {
      #pars = paste0(paste0("^", pars, "$"), collapse="|")
      return(bayesplot::mcmc_combo(x$samples, pars = pars))
    } else {
      return(bayesplot::mcmc_combo(x$samples, regex_pars = pars))
    }
  }
}



#' Summarise mcpfit
#'
#' A simple summary of mcpfit objects. It will be updated later.
#'
#' @aliases summary summary.mcpfit
#' @param object \code{mcpfit} object.
#' @param width The width of the Higest Density Interval.
#' @param ... Currently ignored.
#' @export
#' @examples
#' \dontrun{
#' summary(fit)
#' }

summary.mcpfit = function(object, width = 0.95, ...) {
  if (!is.null(object$samples)) {
    # TO DO: Temporarily suppress warnings about tidyr1.0::unnest() until tidybayes 1.2 is out.
    suppressWarnings(
      object$samples %>%
        tidybayes::tidy_draws() %>%
        tidyr::pivot_longer(-tidyselect::starts_with(".")) %>%
        dplyr::group_by(name) %>%
        tidybayes::mean_hdci(value, .width = width) %>%
        dplyr::rename(mean = value) %>%
        dplyr::select(-.point, -.width, -.interval)
    )
  }
  else {
    message("No samples. Nothing to summarise.")
  }
}

#' Print mcpfit
#'
#' @aliases print print.mcpfit
#' @param x \code{mcpfit} object.
#' @param ... Currently ignored.
#' @export
print.mcpfit = function(x, ...) {
  if (!is.null(x$samples)) {
    print(summary(x))
  }
  else {
    message("No samples. Nothing to summarise.")
  }
}
