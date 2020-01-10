#' Plot full fits
#'
#' Plot prior or posterior model draws on top of data. Use `plot_pars` to
#' plot individual parameter estimates.
#'
#' @aliases plot plot.mcpfit
#' @param x An \code{\link{mcpfit}} object
#' @param facet_by String. Name of a varying group.
#' @param lines Positive integer or `FALSE`. Number of lines (posterior
#'   draws). FALSE or `lines = 0` plots no lines. Note that lines always plot
#'   fitted values - not predicted. For prediction intervals, see the
#'   `q_predict` argument.
#' @param geom_data String. One of "point" (default), "line" (good for time-series),
#'   or FALSE (don not plot).
#' @param cp_dens TRUE/FALSE. Plot posterior densities of the change point(s)?
#'   Currently does not respect `facet_by`. This will be added in the future.
#' @param q_fit Whether to plot quantiles of the posterior (fitted value).
#'   * \strong{TRUE:} Add 2.5% and 97.5% quantiles. Corresponds to
#'       `q_fit = c(0.025, 0.975)`.
#'   * \strong{FALSE (default):} No quantiles
#'   * A vector of quantiles. For example, `quantiles = 0.5`
#'       plots the median and `quantiles = c(0.2, 0.8)` plots the 20% and 80%
#'       quantiles.
#' @param q_predict Same as `q_fit`, but for the prediction interval.
#' @param rate Boolean. For binomial models, plot on raw data (`rate = FALSE`) or
#'   response divided by number of trials (`rate = TRUE`). If FALSE, linear
#'   interpolation on trial number is used to infer trials at a particular x.
#' @param prior TRUE/FALSE. Plot using prior samples? Useful for `mcp(..., sample = "both")`
#' @param which_y What to plot on the y-axis. One of
#'
#'   * `"ct"`: The central tendency which is often the mean after applying the
#'     link function (default).
#'   * `"sigma"`: The variance
#'   * `"ar1"`, `"ar2"`, etc. depending on which order of the autoregressive
#'     effects you want to plot.
#'
#' @param ... Currently ignored.
#' @details
#'   `plot()` uses `fit$simulate()` on posterior samples. These represent the
#'   (joint) posterior distribution.
#' @return A \pkg{ggplot2} object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_point facet_wrap
#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
#' @importFrom dplyr .data
#' @export
#' @examples
#' # Typical usage. ex_fit is an mcpfit object.
#' plot(ex_fit)
#' plot(ex_fit, prior = TRUE)  # The prior
#'
#' \donttest{
#' plot(ex_fit, lines = 0, q_fit = TRUE)  # 95% HDI without lines
#' plot(ex_fit, q_predict = c(0.1, 0.9))  # 80% prediction interval
#' plot(ex_fit, which_y = "sigma", lines = 100)  # The variance parameter on y
#' }
#'
#' # Show a panel for each varying effect
#' # plot(fit, facet_by = "my_column")
#'
#' # Customize plots using regular ggplot2
#' library(ggplot2)
#' plot(ex_fit) + theme_bw(15) + ggtitle("Great plot!")
#'
plot.mcpfit = function(x,
                       facet_by = NULL,
                       lines = 25,
                       geom_data = "point",
                       cp_dens = TRUE,
                       q_fit = FALSE,
                       q_predict = FALSE,
                       rate = TRUE,
                       prior = FALSE,
                       which_y = "ct",
                       ...) {

  # Just for consistent naming in mcp
  fit = x

  # Check arguments
  if (class(fit) != "mcpfit")
    stop("Can only plot mcpfit objects. x was class: ", class(fit))

  if (!coda::is.mcmc.list(fit$mcmc_post) & !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

  if (lines != FALSE) {
    check_integer(lines, "lines", lower = 1)
  } else {
    lines = 0
  }

  if (lines > 1000) {
    lines = 1000
    warning("Setting `lines` to 1000 (maximum).")
  }

  if (!geom_data %in% c("point", "line", FALSE))
    stop("`geom_data` has to be one of 'point', 'line', or FALSE.")

  if (!is.logical(cp_dens))
    stop("`cp_dens` must be TRUE or FALSE.")

  if (all(q_fit == TRUE))
    q_fit = c(0.025, 0.975)

  if (all(q_predict == TRUE))
    q_predict = c(0.025, 0.975)

  if (!is.logical(q_fit) & !is.numeric(q_fit))
    stop("`q_fit` has to be TRUE, FALSE, or a vector of numbers.")

  if (is.numeric(q_fit) & (any(q_fit > 1) | any(q_fit < 0)))
    stop ("All `q_fit` have to be between 0 (0%) and 1 (100%).")

  if (!is.logical(q_predict) & !is.numeric(q_predict))
    stop("`q_predict` has to be TRUE, FALSE, or a vector of numbers.")

  if (is.numeric(q_predict) & (any(q_predict > 1) | any(q_predict < 0)))
    stop ("All `q_predict` have to be between 0 (0%) and 1 (100%).")

  if (!is.logical(rate))
    stop("`rate` has to be TRUE or FALSE.")

  if (!is.logical(prior))
    stop("`prior` must be either TRUE or FALSE.")

  # Is facet_by a random/nested effect?
  varying_groups = logical0_to_null(unique(stats::na.omit(fit$.other$ST$cp_group_col)))
  if (!is.null(facet_by)) {
    if (!facet_by %in% varying_groups)
      stop("`facet_by` is not a data column used as varying grouping.")
  }


  # R CMD Check wants a global definition of ".". The formal way of doing it is
  # if(getRversion() >= "2.15.1") utils::globalVariables(".")
  # but that makes the tests fail.
  . = "ugly fix to please R CMD check"

  # Select posterior/prior samples
  samples = get_samples(fit, prior = prior)

  # General settings
  xvar = rlang::sym(fit$pars$x)
  yvar = rlang::sym(fit$pars$y)
  simulate = fit$simulate  # To make a function call work later
  dens_threshold = 0.001  # Do not display change point densities below this threshold.
  dens_height = 0.2  # proportion of plot y-axis that the change point density makes up
  if (all(q_fit == FALSE) & all(q_predict == FALSE)) {
    HDI_SAMPLES = lines
  } else {
    HDI_SAMPLES = 1000 # Maximum number of draws to use for computing HDI
    HDI_SAMPLES = min(HDI_SAMPLES, length(samples) * nrow(samples[[1]]))
  }
  is_arma = length(fit$pars$arma) > 0
  if (is_arma & (all(q_fit != FALSE) | all(q_predict != FALSE)))
    message("Plotting ar() with quantiles can be slow. Raise an issue at GitHub (or thumb-up existing ones) if you need this.")
  if (!is.null(facet_by)) {
    n_facet_levels = length(unique(fit$data[, facet_by]))
  } else {
    n_facet_levels = 1
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

  # Get x-coordinates to evaluate simulate (etc.) at
  eval_at = get_eval_at(fit, facet_by)

  # First, let's get all the predictors in shape for simulate
  if (fit$family$family != "binomial") {
    samples = samples %>%
      tidyr::expand_grid(!!xvar := eval_at)  # correct name of x-var
  } else if (fit$family$family == "binomial") {
    if (!is.null(facet_by) & rate == FALSE)
      stop("Plot with rate = FALSE not implemented for varying effects (yet).")

    # Interpolate trials for binomial at the values in "eval_at"
    # to make sure that there are actually values to plot
    #interpolated_trials = suppressWarnings(stats::approx(x = fit$data[, fit$pars$x], y = fit$data[, fit$pars$trials], xout = eval_at)$y)
    interpolated_trials = stats::approx(x = fit$data[, fit$pars$x], y = fit$data[, fit$pars$trials], xout = eval_at)$y %>%
      suppressWarnings() %>%
      round()  # Only integers

    samples = samples %>%
      tidyr::expand_grid(!!xvar := eval_at) %>%  # correct name of x-var
      dplyr::mutate(!!fit$pars$trials := rep(interpolated_trials, nrow(samples)))
  }

  # For ARMA prediction, we need the raw data
  # We know that eval_at is the same length as nrow(data), so we can simply add corresponding data$y for each draw
  if (is_arma) {
    samples = dplyr::mutate(samples, ydata = rep(fit$data[, fit$pars$y], HDI_SAMPLES * n_facet_levels))
  }

  # Predict y from model by adding fitted/predicted draws (vectorized)
  if (lines > 0 | (any(q_fit != FALSE))) {
    samples = samples %>%
      dplyr::mutate(!!yvar := rlang::exec(simulate, !!!., type = "fitted", rate = rate, which_y = which_y, add_attr = FALSE))
  }
  if (any(q_predict != FALSE)) {
    samples = samples %>%
      dplyr::mutate(predicted_ = rlang::exec(simulate, !!!., type = "predict", rate = rate, add_attr = FALSE))
  }



  ###########
  # PLOT IT #
  ###########
  # If this is a binomial rate, divide by the number of trials
  if (fit$family$family == "binomial" & rate == TRUE) {
    fit$data[, fit$pars$y] = fit$data[, fit$pars$y] / fit$data[, fit$pars$trials]
  }

  # Initiate plot.
  gg = ggplot(fit$data, aes_string(x = fit$pars$x, y = fit$pars$y))
  if (which_y == "ct") {
    if (geom_data == "point")
      gg = gg + geom_point()
    if (geom_data == "line")
      gg = gg + geom_line()
  }

  # Add lines?
  if (lines > 0) {
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
  if ((any(q_fit != FALSE))) {
    samples_fit = dplyr::mutate(samples, y_quant = !!yvar)
    gg = gg + geom_quantiles(samples_fit, q_fit, xvar, facet_by, color = "red")
  }
  if (any(q_predict != FALSE)) {
    samples_predict = dplyr::mutate(samples, y_quant = .data$predicted_)
    gg = gg + geom_quantiles(samples_predict, q_predict, xvar, facet_by, color = "green4")
  }

  # Add change point densities?
  if (cp_dens == TRUE & length(fit$model) > 1) {
    # The scale of the actual plot (or something close enough)
    if (which_y == "ct" & geom_data != FALSE) {
      y_data_max = max(fit$data[, fit$pars$y])
      y_data_min = min(fit$data[, fit$pars$y])
    } else if (any(q_predict != FALSE)) {
      y_data_max = max(samples$predicted_)
      y_data_min = min(samples$predicted_)
    } else if (as.character(yvar) %in% names(samples)) {
      y_data_max = max(dplyr::pull(samples, as.character(yvar)))
      y_data_min = min(dplyr::pull(samples, as.character(yvar)))
    } else {
      stop("Failed to draw change point density for this plot. Please raise an error on GitHub.")
    }

    # Function to get density for each grouping in the dplyr pipes below.
    density_xy = function(x) {
      tmp = stats::density(x)
      df = data.frame(x_dens = tmp$x, y_dens = tmp$y) %>%
        dplyr::filter(.data$y_dens > dens_threshold) %>%
        return()
    }

    # Get and group samples to be used for computing density
    cp_regex = "^cp_[0-9]+$"
    cp_dens_xy = get_samples(fit, prior = prior) %>%  # Use all samples for this
      tidybayes::gather_draws(!!rlang::sym(cp_regex), regex = TRUE) %>%
      dplyr::group_by(.data$.chain, add = TRUE) %>%

      # Get density as x-y values by chain and change point number.
      dplyr::summarise(dens = list(density_xy(.data$.value))) %>%
      tidyr::unnest(cols = c(.data$dens)) %>%

      # Then scale to plot.
      dplyr::mutate(y_dens = y_data_min + .data$y_dens * (y_data_max - y_data_min) * dens_height / max(.data$y_dens))

    # Add cp density to plot
    gg = gg + ggplot2::geom_line(
      data = cp_dens_xy,
      mapping = aes(
        x = .data$x_dens,
        y = .data$y_dens,
        color = .data$.chain,
        group = interaction(.data$.variable, .data$.chain))) +
      ggplot2::theme(legend.position = "none")

  }

  # Add faceting?
  if (!is.null(facet_by)) {
    gg = gg + facet_wrap(paste0("~", facet_by))
  }

  # Add better y-labels
  if (fit$family$family == "bernoulli" | (fit$family$family == "binomial" & rate == TRUE))
    gg = gg + ggplot2::labs(y = paste0("Probability of success for ", fit$pars$y))
  if (which_y != "ct")
    gg = gg + ggplot2::labs(y = which_y)

  return(gg)
}


#' Return a geom_line representing the quantiles
#'
#' Called by `plot.mcpfit`.
#'
#' @aliases geom_quantiles
#' @keywords internal
#' @inheritParams plot.mcpfit
#' @param q Quantiles
#' @param ... Arguments passed to geom_line
#' @return A `ggplot2::geom_line` object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
geom_quantiles = function(samples, q, xvar, facet_by, ...) {
  # Trick to declare no facet = common group for all
  if (length(facet_by) == 0)
    facet_by = xvar

  # First: add quantiles column
  data_quantiles = samples %>%
    tidyr::expand_grid(quant = q) %>%

    # Now compute the quantile for each parameter, quantile, and (optionally) facet:
    dplyr::group_by(!!xvar, .data$quant) %>%
    dplyr::group_by(!!rlang::sym(facet_by), add = TRUE) %>%
    dplyr::summarise(
      y = stats::quantile(.data$y_quant, probs = .data$quant[1])
    )

  # Return geom
  geom = ggplot2::geom_line(
    mapping = aes(
      y = .data$y,
      group = .data$quant
    ),
    data = data_quantiles,
    lty = 2,
    lwd = 0.7,
    ...)
  return(geom)
}




#' Plot individual parameters
#'
#' Plot many types of plots of parameter estimates. See examples for typical use
#' cases.
#'
#' @aliases plot_pars
#' @param fit An \code{\link{mcpfit}} object.
#' @param pars Character vector. One of:
#'   * Vector of parameter names.
#'   * \emph{"population" (default):} plots all population parameters.
#'   * \emph{"varying":} plots all varying effects. To plot a particular varying
#'       effect, use `regex_pars = "^name"`.
#' @param regex_pars Vector of regular expressions. This will typically just be
#'   the beginning of the parameter name(s), i.e., "^cp_" plots all change
#'   points, "^my_varying" plots all levels of a particular varying effect, and
#'   "^cp_|^my_varying" plots both.
#' @param type String or vector of strings. Calls `bayesplot::mcmc_>>type<<()`.
#'   Common calls are "combo", "trace", and "dens_overlay". Current options include
#'   'acf', 'acf_bar', 'areas', 'areas_ridges', 'combo', 'dens', 'dens_chains',
#'   'dens_overlay', 'hist', 'intervals', 'rank_hist', 'rank_overlay', 'trace',
#'   'trace_highlight', and 'violin".
#' @param ncol Number of columns in plot. This is useful when you have many
#'   parameters and only one plot `type`.
#' @param prior TRUE/FALSE. Plot using prior samples? Useful for `mcp(..., sample = "both")`
#'
#'@details
#'   For other `type`, it calls `bayesplot::mcmc_type()`. Use these
#'   directly on `fit$mcmc_post` or `fit$mcmc_prior` if you want finer
#'   control of plotting, e.g., `bayesplot::mcmc_dens(fit$mcmc_post)`. There
#'   are also a number of useful plots in the \pkg{coda} package, i.e.,
#'   `coda::gelman.plot(fit$mcmc_post)` and `coda::crosscorr.plot(fit$mcmc_post)`
#'
#'   In any case, if you see a few erratic lines or parameter estimates, this is
#'   a sign that you may want to increase argument 'adapt' and 'iter' in \code{\link{mcp}}.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A \pkg{ggplot2} object.
#' @import patchwork
#' @export
#' @examples
#' # Typical usage. ex_fit is an mcpfit object.
#' plot_pars(ex_fit)
#'
#' # More options
#' plot_pars(ex_fit, regex_pars = "^cp_")  # Plot only change points
#' plot_pars(ex_fit, pars = c("int_3", "time_3"))  # Plot these parameters
#' plot_pars(ex_fit, type = c("trace", "violin"))  # Combine plots
#'
#' \dontrun{
#' # Some plots only take pairs. hex is good to assess identifiability
#' plot_pars(ex_fit, type = "hex", pars = c("cp_1", "time_2"))
#'
#' # Visualize the priors:
#' plot_pars(ex_fit, prior = TRUE)
#'
#' # Useful for varying effects:
#' # plot_pars(my_fit, pars = "varying", ncol = 3)  # plot all varying effects
#' # plot_pars(my_fit, regex_pars = "my_varying", ncol = 3)  # plot all levels of a particular varying
#'
#' # Customize multi-column ggplots using "*" instead of "+" (patchwork)
#' library(ggplot2)
#' plot_pars(ex_fit, type = c("trace", "dens_overlay")) * theme_bw(10)
#' }

plot_pars = function(fit,
                     pars = "population",
                     regex_pars = character(0),
                     type = "combo",
                     ncol = 1,
                     prior = FALSE) {

  # Check arguments
  if (class(fit) != "mcpfit")
    stop("Can only plot mcpfit objects. x was class: ", class(fit))

  if (!coda::is.mcmc.list(fit$mcmc_post) & !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

  if (!is.character(pars) | !is.character(regex_pars))
    stop("`pars` and `regex_pars` has to be string/character.")

  if (any(c("population", "varying") %in% pars) & length(pars )> 1)
    stop("`pars` cannot be a vector that contains multiple elements AND 'population' or 'varying'.")

  if (any(c("hex", "scatter") %in% type) & (length(pars) != 2 | length(regex_pars) > 0))
    stop("`type` = 'hex' or 'scatter' takes exactly two parameters which must be provided via the `pars` argument")

  if ("combo" %in% type & length(type) > 1)
    stop("'combo' type cannot be combined with other types. Replace 'combo' with the types you want combo\'ed")

  check_integer(ncol, name = "ncol", lower = 1)

  if (!is.logical(prior))
    stop("`prior` must be either TRUE or FALSE.")

  # Get posterior/prior samples
  samples = get_samples(fit, prior = prior)

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
  if ("combo" %in% type)
    type = c("dens_overlay", "trace")

  # TO DO: a lot of eval(parse()) here. Is there a more built-in method?
  #types = c("dens", "dens_overlay", "trace", "areas")
  bayesplot::available_mcmc()  # quick fix to make R CMD Check happy that bayesplot is imported
  takes_facet = c("areas", "dens", "dens_overlay", "trace", "hist", "intervals", "trace", "trace_highlight", "violin")
  for (this_type in type) {
    this_facet = ifelse(this_type %in% takes_facet, paste0(", facet_args = list(ncol = ", ncol, ")"), "")
    command = paste0("plot_", this_type, " = bayesplot::mcmc_", this_type, "(samples, pars = pars, regex_pars = regex_pars", this_facet, ")")
    eval(parse(text = command)) + ggplot2::theme(strip.placement = NULL)

  }

  # Select which to plot
  if (length(type) == 1) {
    return_plot = eval(parse(text = paste0("plot_", type)))
    return_plot = return_plot
  } else {
    # Use patchwork
    command = paste0(paste0("plot_", type), collapse = " + ")
    return_plot = eval(parse(text = command))
  }

  # Return
  return_plot & ggplot2::theme(
    strip.placement = NULL,  # fixes bug: https://github.com/thomasp85/patchwork/issues/132
    legend.position = "none"  # no legend on chains. Takes up too much space
  )
}



#' Get a list of x-coordinates to evaluate fit$simulate at
#'
#' Solves two problems: if setting the number of points too high, the
#' function becomes slow. If setting it too low, the posterior at large intercept-
#' changes at change points look discrete, because they are evaluated at very
#' few x in that interval.
#'
#' This function makes a vector of x-values with large spacing in general,
#' but finer resolution at change points.
#'
#' @aliases get_eval_at
#' @keywords internal
#' @inheritParams plot.mcpfit
#' @param fit An mcpfit object.
get_eval_at = function(fit, facet_by) {
  # If there are ARMA terms, evaluate at the data
  if (length(fit$pars$arma) > 0) {
    return(c(fit$data[, fit$pars$x]))
  }

  # Set resolutions in general and for change points
  X_RESOLUTION_ALL = 100  # Number of points to evaluate at x
  X_RESOLUTION_CP = 600
  X_RESOLUTION_FACET = 300
  CP_INTERVAL = 0.9  # HDI interval width
  xmin = min(fit$data[, fit$pars$x])  # smallest observed X
  xmax = max(fit$data[, fit$pars$x])

  # Just give up for faceting and prior-plots (usually very distributed change points)
  # and return a reasonable resolution
  if (!is.null(facet_by) | is.null(fit$mcmc_post)) {
    eval_at = seq(xmin, xmax, length.out = X_RESOLUTION_FACET)
    return(eval_at)
  }

  # Make the coarse resolution
  eval_at = seq(xmin, xmax, length.out = X_RESOLUTION_ALL)

  # Add the finer resolution for each change point
  cp_vars = paste0("cp_", seq_len(length(fit$model) - 1))  # change point columns
  cp_hdis = fixef(fit, width = CP_INTERVAL)  # get the intervals
  cp_hdis = cp_hdis[cp_hdis$name %in% cp_vars, ]  # select change points
  for (i in seq_len(nrow(cp_hdis))) {
    x_proportion = (cp_hdis$upper[i] - cp_hdis$lower[i]) / (xmax - xmin)  # how big a section of x is this CP's HDI?
    length.out = ceiling(X_RESOLUTION_CP * x_proportion)  # number of x-points to add
    add_this = seq(from = max(xmin, cp_hdis$lower[i]),  # hack to avoid bug for relative cps. TO DO
                   to = min(xmax, cp_hdis$upper[i]),  # hack to avoid bug for relative cps. TO DO
                   length.out = length.out)
    eval_at = c(eval_at, add_this)
  }

  return(sort(eval_at))
}
