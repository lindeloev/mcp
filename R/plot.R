#' Plot mcpfit
#'
#' Plotting posterior fitted lines on top of data (\code{plot(fit)}) or many
#' types of plots of parameter estimates (typically \code{plot(fit, "combo")}).
#' See examples for typical use cases.
#'
#' @aliases plot plot.mcpfit
#' @param x An mcpfit object
#' @param type String or vector of strings. Calls \code{bayesplot::mcmc_type()}.
#'   Common calls are "combo", "trace", and "dens_overlay". Current options include
#'   'acf', 'acf_bar', 'areas', 'areas_ridges', 'combo', 'dens', 'dens_chains',
#'   'dens_overlay', 'hist', 'intervals', 'rank_hist', 'rank_overlay', 'trace',
#'   'trace_highlight', and 'violin".
#'
#' @param pars Character vector. One of:
#'   \itemize{
#'     \item Vector of parameter names.
#'     \item \strong{"population" (default):} plots all population parameters.
#'     \item \strong{"varying":} plots all varying effects. To plot a particular varying
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
#' @param lines Positive integer or \code{FALSE}. Number of lines (posterior
#'   draws) to use when \code{type = "segments"}. FALSE (or \code{lines = 0})
#'   plots no lines.
#' @param quantiles Whether to plot quantiles.
#'   \itemize {
#'     \item \strong{TRUE:} Add 2.5% and 97.5% quantiles. Corresponds to
#'       \code{quantiles = c(0.025, 0.975)}.
#'     \item \strong{FALSE (default):} No quantiles
#'     \item A vector of quantiles. For example, \code{quantiles = 0.5}
#'       plots the median and \code{quantiles = c(0.2, 0.8)} plots the 20% and 80%
#'       quantiles.
#'   }
#' @param ... Currently ignored.
#' @details
#'   For \code{type = "segments"}, it uses \code{fit$func_y} with \code{draws}
#'   posterior samples. These represent the joint posterior distribution of
#'   parameter values.
#'
#'   For other \code{type}, it calls \code{bayesplot::mcmc_type()}. Use these
#'   directly on \code{fit$mcmc_post} or \code{fit$mcmc_prior} if you want finer
#'   control of plotting, e.g., \code{bayesplot::mcmc_dens(fit$mcmc_post)}. There
#'   are also a number of useful plots in the \pkg{coda} package, i.e.,
#'   \code{coda::gelman.plot(fit$mcmc_post)} and \code{coda::crosscorr.plot(fit$mcmc_post)}
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
#' plot(fit, lines = 50, rate = FALSE)  # more lines, for binomial.
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

plot.mcpfit = function(x,
                       type = "segments",

                       # Arguments for bayesplot
                       pars = "population",
                       regex_pars = character(0),
                       ncol = 1,

                       # Arguments for "segments"
                       facet_by = NULL,
                       rate = TRUE,
                       lines = 25,
                       quantiles = FALSE,
                       ...) {

  # Check arguments
  if (class(x) != "mcpfit")
    stop("Can only plot mcpfit objects. x was class: ", class(x))

  if (!is.character(type))
    stop("`type` has to be string/character.")

  if (length(type) > 1 & "segments" %in% type)
    stop("type = 'segments' can only be plotted alone - not in in combination with others.")

  if (lines != FALSE)
    check_integer(lines, "lines", lower = 1)

  if (lines > 1000) {
    lines = 1000
    warning("Setting `lines` to 1000 (maximum).")
  }

  if (all(quantiles == TRUE))
    quantiles = c(0.025, 0.975)

  if (!is.logical(quantiles) & !is.numeric(quantiles))
    stop("`quantiles` has to be TRUE, FALSE, or a vector of numbers.")

  if (is.numeric(quantiles) & (any(quantiles > 1) | any(quantiles < 0)))
      stop ("all `quantiles` have to be between 0 (0%) and 1 (100%).")

  if (!is.logical(rate))
    stop("`rate` has to be TRUE or FALSE.")

  if (!coda::is.mcmc.list(x$mcmc_post) & !coda::is.mcmc.list(x$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

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
    plot_segments(x, facet_by, rate, lines, quantiles, ...)
  } else {
    plot_bayesplot(x, type, pars, regex_pars, ncol, ...)
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
plot_segments = function(fit,
                         facet_by = NULL,
                         rate = TRUE,
                         lines = 25,
                         quantiles = FALSE,
                         ...) {

  # Select posterior/prior samples
  samples = get_samples(fit)

  # General settings
  xvar = rlang::sym(fit$pars$x)
  yvar = rlang::sym(fit$pars$y)
  func_y = fit$func_y
  if (all(quantiles == FALSE) & is.numeric(lines)) {
    HDI_SAMPLES = lines
  } else {
    HDI_SAMPLES = 1000 # Maximum number of draws to use for computing HDI
    HDI_SAMPLES = min(HDI_SAMPLES, length(samples) * nrow(samples[[1]]))
  }


  #################
  # GET PLOT DATA #
  #################
  regex_pars_pop = paste0(fit$pars$population, collapse="|")

  # No faceting
  if (is.null(facet_by)) {
    samples = samples %>%
      tidybayes::spread_draws(!!rlang::sym(regex_pars_pop), regex = TRUE)

  } else {
    # Prepare for faceting
    # Read more about this weird syntax at https://github.com/mjskay/tidybayes/issues/38
    varying_by_facet = stats::na.omit(fit$.other$ST$cp_group[stringr::str_detect(fit$.other$ST$cp_group, paste0("_", facet_by, "$"))])
    varying_by_facet = paste0(varying_by_facet, collapse="|")

    samples = samples %>%
      tidybayes::spread_draws(!!rlang::sym(regex_pars_pop),
                              (!!rlang::sym(varying_by_facet))[!!rlang::sym(facet_by)],
                              regex = TRUE)
  }

  # Remove some samples
  samples = tidybayes::sample_draws(samples, n = HDI_SAMPLES)  # TO DO: use spread_draws(n = draws) when tidybayes 1.2 is out

  # Get x-coordinates to evaluate func_y (etc.) at
  eval_at = get_eval_at(fit, facet_by)

  # First, let's get all the predictors in shape for func_y
  if (fit$family$family != "binomial") {
    samples = samples %>%
      tidyr::expand_grid(!!xvar := eval_at)  # correct name of x-var
  } else if (fit$family$family == "binomial") {
    if (!is.null(facet_by) & rate == FALSE)
      stop("Plot with rate = FALSE not implemented for varying effects (yet).")

    # Interpolate trials for binomial at the values in "eval_at"
    interpolated_trials = round(stats::approx(fit$data[, fit$pars$trials], xout = eval_at)$y)
    samples = samples %>%
      tidyr::expand_grid(!!xvar := eval_at) %>%  # correct name of x-var
      dplyr::mutate(!!fit$pars$trials := rep(interpolated_trials, nrow(samples)))
  }

  # Predict y from model
  samples = samples %>%
    # Add fitted draws (vectorized)
    dplyr::mutate(!!yvar := purrr::invoke(func_y, ., type = "fitted", rate = rate))



  ###########
  # PLOT IT #
  ###########
  # If this is a binomial rate, divide by the number of trials
  if (fit$family$family == "binomial" & rate == TRUE) {
    fit$data[, fit$pars$y] = fit$data[, fit$pars$y] / fit$data[, fit$pars$trials]
  }

  # Initiate plot
  gg = ggplot(fit$data, aes_string(x = fit$pars$x, y = fit$pars$y)) +
    geom_point()

  # Add lines?
  if (lines != FALSE) {
    # Sample the desired number of lines
    data_lines = samples %>%
      tidybayes::sample_draws(lines) %>%
      dplyr::mutate(
        # Add line ID to separate lines in the plot.
        line = !!xvar == min(!!xvar),
        line = cumsum(.data$line)
      )

    # Plot it
    gg = gg + geom_line(aes(group = .data$line), data = data_lines, color = grDevices::rgb(0.5, 0.5, 0.5, 0.4))
  }

  # Add quantiles?
  if (is.numeric(quantiles)) {
    # For each quantile and each x, add quantile of !!yvar
    data_quantiles = samples %>%
      tidyr::expand_grid(quant = quantiles) %>%
      dplyr::group_by(!!xvar, .data$quant)

    # (... and for each group, if requested)
    if (!is.null(facet_by)) {
     data_quantiles = data_quantiles %>%
       dplyr::group_by(!!rlang::sym(facet_by), add = TRUE)
    }

    # Continue...
    data_quantiles = data_quantiles %>%
      dplyr::summarise(
        y = stats::quantile(!!yvar, probs = .data$quant[1])
      )

    # Add quantiles to plot
    gg = gg +
      geom_line(aes(y = .data$y, group = .data$quant), data = data_quantiles, lty = 2, lwd = 0.7, color = "red")
  }

  # Add faceting?
  if (!is.null(facet_by)) {
    gg = gg + facet_wrap(paste0("~", facet_by))
  }

  # Add better y-label for the rate plot
  if (fit$family$family == "bernoulli" | (fit$family$family == "binomial" & rate == TRUE))
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
plot_bayesplot = function(fit,
                          type,
                          pars = "population",
                          regex_pars = character(0),
                          ncol = 1) {

  # Get posterior/prior samples
  samples = get_samples(fit)

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
    command = paste0("plot_", this_type, " = bayesplot::mcmc_", this_type, "(samples, pars = pars, regex_pars = regex_pars", this_facet, ")")
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



#' Get a list of x-coordinates to evaluate fit$func_y at
#'
#' Solves two problems: if setting the number of points too high, the
#' function becomes slow. If setting it too low, the posterior at large intercept-
#' changes at chnage points look discrete, because they are evaluated at very
#' few x in that interval.
#'
#' This function makes a vector of x-values with large spacing in general,
#' but finer resolution at change points.
#'
#' @aliases get_eval_at
#' @inheritParams plot.mcpfit
get_eval_at = function(fit, facet_by) {
  # Set resolutions in general and for change points
  X_RESOLUTION_ALL = 100  # Number of points to evaluate at x
  X_RESOLUTION_CP = 600
  X_RESOLUTION_FACET = 300
  CP_INTERVAL = 0.9  # HDI interval width
  xmin = min(fit$data[, fit$pars$x])  # smallest observed X
  xmax = max(fit$data[, fit$pars$x])

  # Just give up for faceting and return a reasonable resolution
  if (!is.null(facet_by) | is.null(fit$mcmc_post)) {
    eval_at = seq(xmin, xmax, length.out = X_RESOLUTION_FACET)
    return(eval_at)
  }

  # Make the coarse resolution
  eval_at = seq(xmin, xmax, length.out = X_RESOLUTION_ALL)

  # Add the finer resolution for each change point
  cp_vars = paste0("cp_", seq_len(length(fit$segments) - 1))  # change point columns
  cp_hdis = fixef(fit, width = CP_INTERVAL)  # get the intervals
  cp_hdis = cp_hdis[cp_hdis$name %in% cp_vars, ]  # select change points
  for (i in seq_len(nrow(cp_hdis))) {
    x_proportion = (cp_hdis$X95[i] - cp_hdis$X5[i]) / (xmax - xmin)  # how big a section of x is this CP's HDI?
    length.out = ceiling(X_RESOLUTION_CP * x_proportion)  # number of x-points to add
    eval_at = c(eval_at, seq(from = cp_hdis$X5[i], to = cp_hdis$X95[i], length.out = length.out))
  }

  return(eval_at)
}
