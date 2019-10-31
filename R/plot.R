#' Plot mcpfit
#'
#' Plotting posterior fitted lines on top of data (\code{plot(fit)}) or many
#' types of plots of parameter estimates (typically \code{plot(fit, "combo")}).
#' See examples for typical use cases.
#'
#' @aliases plot plot.mcpfit
#' @param x An mcpfit object
#' @param type String or vector of strings. Calls \code{bayesplot::mcmc_*type*()}.
#'   Common calls are "combo", "trace", and "dens_overlay". Current options include
#'   'acf', 'acf_bar', 'areas', 'areas_ridges', 'combo', 'dens', 'dens_chains',
#'   'dens_overlay', 'hist', 'intervals', 'rank_hist', 'rank_overlay', 'trace',
#'   'trace_highlight', and 'violin".
#' @param draws Positive integer. Number of posterior draws to use when type = "segments".
#' @param pars Character vector. One of:
#'   \itemize{
#'     \item Vector of parameter names.
#'     \item "population" (default) plots all population parameters.
#'     \item "varying" plots all varying effects. To plot a particular varying
#'       effect, use \code{regex_pars = "^name"}.
#'   }
#' @param regex_pars Vector of regular expressions. This will typically just be
#'   the beginning of the parameter name(s), i.e., "^cp_" plots all change
#'   points, "^my_varying" plots all levels of a particular varying effect, and
#'   "^cp_|^my_varying" plots both.
#' @param facet_by String. Name of a varying group.
#'   \code{facet_by} only applies for \code{type = "segments"}
#' @param rate Boolean. For binomial models, plot on raw data (\code{rate = FALSE}) or
#'   response divided by number of trials (\code{rate = TRUE}). If FALSE, linear
#'   interpolation on trial number is used to infer trials at a particular x.
#'   \code{rate} only applies for \code{type = "segments"}
#' @param ncol Number of columns in plot. This is useful when you have many
#'   parameters and only one plot \code{type}.
#'   \code{ncol} only when \code{type != "segments"}
#' @param ... Currently ignored.
#' @details
#'   For \code{type = "segments"}, it uses \code{fit$func_y} with \code{draws}
#'   posterior samples. These represent the joint posterior distribution of
#'   parameter values.
#'
#'   For other \code{type}, it calls \code{bayesplot::mcmc_*type*()}. Use these
#'   directly on \code{fit$samples} if you want finer control of plotting, e.g.,
#'   \code{bayesplot::mcmc_dens(fit$samples)}. There are also a number of useful
#'   plots in the \pkg{coda} package, i.e., \code{coda::gelman.plot(fit$samples)}
#'   and \code{coda::crosscorr.plot(fit$samples)}
#'
#'   In any case, if you see a few erratic lines or parameter estimates, this is
#'   a sign that you may want to increase arguments 'adapt', 'update', and
#'   'iter' in \code{\link{mcp}}.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A \code{ggplot2} object.
#' @export
#' @examples
#' \dontrun{
#' # Plots segments (default)
#' plot(fit)
#' plot(fit, draws = 50, rate = FALSE)  # more draws. for binomial.
#' plot(fit, facet_by = "my_varying")  # varying effects
#'
#' # Plot parameter estimates
#' plot(fit, "combo")
#' plot(fit, "combo", pars = "varying", ncol = 3)  # plot all varying effects
#' plot(fit, "combo", regex_pars = "my_varying", ncol = 3)  # plot all levels of a particular varying
#'
#' # More options for parameter estimates
#' plot(fit, "combo", pars = c("var1", "var2", "var3"), regex_pars = "^my_varying")
#' plot(fit, c("areas", "intervals"))
#'
#' # Some plots only take pairs
#' plot(fit, "hex", pars = c("var1", "var2"))
#'
#'
#' # Customize one-column plots using regular ggplot2
#' plot(fit) + theme_bw(15) + ggtitle("Great plot!")
#'
#' # Customize two-column plots using the "patchwork" package.
#' plot(fit, c("trace", "dens_overlay")) * theme_bw(10)
#' }

plot.mcpfit = function(x, type="segments", pars="population", regex_pars = character(0), facet_by = NULL, rate = TRUE, draws=25, ncol = 1, ...) {
  # Check arguments
  if (class(x) != "mcpfit")
    stop("Can only plot mcpfit objects. x was class: ", class(x))

  if (!is.character(type))
    stop("`type` has to be string/character.")

  if (length(type) > 1 & "segments" %in% type)
    stop("type = 'segments' can only be plotted alone - not in in combination with others.")

  if (!check_integer(draws, "draws", lower = 1))
    stop("`draws` has to be a positive integer.")

  if (!is.logical(rate))
    stop("`rate` has to be TRUE or FALSE")

  # Include test of whether this is a random/nested effect
  varying_groups = logical0_to_null(unique(stats::na.omit(x$.other$ST$cp_group_col)))
  if (!is.null(facet_by)) {
    if (!facet_by %in% varying_groups)
      stop("facet_by is not a data column used as varying grouping.")
  }

  if (!is.character(pars) | !is.character(regex_pars))
    stop("pars and regex_pars has to be string/character.")

  if (any(c("population", "varying") %in% pars) & length(pars )> 1)
    stop("pars cannot be a vector that contains multiple elements AND 'population' or 'varying'.")


  if (any(c("hex", "scatter") %in% type) & (length(pars) != 2 | length(regex_pars) > 0))
    stop("`type` = 'hex' or 'scatter' takes exactly two parameters which must be provided via the `pars` argument")


  # Call underlying plot functions
  if ("segments" %in% type) {
    plot_segments(x, draws, facet_by, rate)
  } else {
    plot_bayesplot(x, type, pars, regex_pars, ncol)
  }
}



#' Plot posterior draws of segments on top of data
#'
#' Call if from \code{\link{plot.mcpfit}} and see more details there
#'
#' @aliases plot_segments
#' @inheritParams plot.mcpfit
#' @param fit An mcpfit object.
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_point facet_wrap
#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
#' @importFrom stats sd
#' @importFrom dplyr .data
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
plot_segments = function(fit, draws = 25, facet_by = NULL, rate = TRUE, ...) {
  X_RESOLUTION = 100  # number of points to evaluate at x

  #################
  # GET PLOT DATA #
  #################

  # We'll need these vars during data-processing-before-ggplot
  func_y = fit$func_y
  eval_at = seq(min(fit$data[, fit$pars$x]),
                max(fit$data[, fit$pars$x]),
                length.out = X_RESOLUTION)
  regex_pars_pop = paste0(fit$pars$population, collapse="|")
  cp_vars = paste0("cp_", seq_len(length(fit$segments) - 1))  # change point columns

  # No faceting
  if (is.null(facet_by)) {
    samples = fit$samples %>%
      tidybayes::spread_draws(!!rlang::sym(regex_pars_pop), regex = TRUE) %>%
      # TO DO: use spread_draws(n = draws) when tidybayes 1.2 is out
      tidybayes::sample_draws(draws)

  } else {
    # Prepare for faceting
    # Read more about this weird syntax at https://github.com/mjskay/tidybayes/issues/38
    varying_by_facet = stats::na.omit(fit$.other$ST$cp_group[stringr::str_detect(fit$.other$ST$cp_group, paste0("_", facet_by))])
    varying_by_facet = paste0(varying_by_facet, collapse="|")

    samples = fit$samples %>%
      tidybayes::spread_draws(!!rlang::sym(regex_pars_pop),
                              (!!rlang::sym(varying_by_facet))[!!rlang::sym(facet_by)],
                              regex = TRUE) %>%
      # TO DO: use spread_draws(n = draws) when tidybayes 1.2 is out
      tidybayes::sample_draws(draws)
  }

  # First, let's get all the predictors in shape for func_y
  if (fit$family$family != "binomial") {
    samples = samples %>%
      tidyr::expand_grid(!!fit$pars$x := eval_at)  # correct name of x-var
  } else if (fit$family$family == "binomial") {
    if (!is.null(facet_by) & rate == FALSE)
      stop("Plot with rate = FALSE not implemented for varying effects (yet).")

    # Interpolate trials for binomial at the values in "eval_at"
    interpolated_trials = round(stats::approx(fit$data[, fit$pars$trials], n = X_RESOLUTION)$y)
    samples = samples %>%
      tidyr::expand_grid(!!fit$pars$x := eval_at) %>%  # correct name of x-var
      dplyr::mutate(!!fit$pars$trials := rep(interpolated_trials, nrow(samples)))
  }

  # Predict y from model
  samples = samples %>%
    # Add fitted draws (vectorized)
    dplyr::mutate(!!fit$pars$y := purrr::invoke(func_y, ., type = "fitted", rate = rate)) %>%

    # Add line ID to separate lines. Mark a new line when "eval_at" repeats
    dplyr::mutate(
      line = !!dplyr::sym(fit$pars$x) == min(eval_at),
      line = cumsum(.data$line)
    )


  ###########
  # PLOT IT #
  ###########
  # If this is a binomial rate, divide by the number of trials
  if (fit$family$family == "binomial" & rate == TRUE) {
    fit$data[, fit$pars$y] = fit$data[, fit$pars$y] / fit$data[, fit$pars$trials]
  }

  if (is.null(facet_by)) {
    # Return plot without faceting
    gg = ggplot(fit$data, aes_string(x = fit$pars$x, y = fit$pars$y)) +
      geom_point() +
      geom_line(aes(group = .data$line), data = samples, color = grDevices::rgb(0.5, 0.5, 0.5, 0.4))
  } else {
    # Return plot with faceting
    gg = ggplot(fit$data, aes_string(x = fit$pars$x, y = fit$pars$y)) +
      geom_point() +
      geom_line(aes(group = .data$line), data = samples, color = grDevices::rgb(0.5, 0.5, 0.5, 0.4)) +
      facet_wrap(paste0("~", facet_by))
  }

  # Add better y-label for the rate plot
  if (rate == TRUE)
    gg = gg + ggplot2::labs(y = paste0("Probability of success for ", fit$pars$y))

  return(gg)
}



#' Plot mcpfit objects with bayesplot
#'
#' Call if from \code{\link{plot.mcpfit}} and see more details there
#'
#' @aliases plot_bayesplot
#' @inheritParams plot.mcpfit
#' @param fit An mcpfit object.
#' @import patchwork
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
plot_bayesplot = function(fit, type, pars = "population", regex_pars = character(0), ncol = 1) {
  # Handle special codes
  if ("population" %in% pars) {
    if (length(regex_pars) == 0) {
      pars = fit$pars$population
    } else {
      # This probably means that the user left pars as default.
      pars = character(0)
    }
  } else if ("varying" %in% pars) {
    # Regex search for varying effects
    regex_pars = paste0("^", fit$pars$varying)
    pars = character(0)
  }

  # Handles combo. Returns a customizable ggplot which "combo" does not.
  if (type == "combo")
    type = c("dens_overlay", "trace")

  # TO DO: a lot of eval(parse()) here. Is there a more built-in method?
  #types = c("dens", "dens_overlay", "trace", "areas")
  takes_facet = c("areas", "dens", "dens_overlay", "trace", "hist", "intervals", "rank_hist", "rank_overlay", "trace", "trace_highlight", "violin")
  for (this_type in type) {
    this_facet = ifelse(this_type %in% takes_facet, paste0(", facet_args = list(ncol = ", ncol, ")"), "")
    command = paste0("plot_", this_type, " = bayesplot::mcmc_", this_type, "(fit$samples, pars = pars, regex_pars = regex_pars", this_facet, ")")
    eval(parse(text = command))
  }

  # Select which to plot
  if (length(type) == 1) {
    return_plot = eval(parse(text = paste0("plot_", type)))
    return_plot = return_plot + ggplot2::theme(legend.position = "none")  # remove legend
    return(return_plot)
  } else {
    command = paste0(paste0("plot_", type), collapse = " + ")
    return_plot = eval(parse(text = command))
    return_plot = return_plot * ggplot2::theme(legend.position = "none")  # remove legend
    return(return_plot)
  }
}
