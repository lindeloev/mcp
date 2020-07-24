#' Plot full fits
#'
#' Plot prior or posterior model draws on top of data. Use `plot_pars` to
#' plot individual parameter estimates.
#'
#' @aliases plot plot.mcpfit
#' @inheritParams pp_eval
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
                       arma = TRUE,
                       nsamples = 2000,
                       ...) {

  # Just for consistent naming in mcp
  fit = x

  # Check arguments
  # The following are checked in pp_eval: q_fit, q_predict, rate
  check_mcpfit(fit)

  if (!coda::is.mcmc.list(fit$mcmc_post) & !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

  if (lines != FALSE) {
    check_integer(lines, "lines", lower = 1)
  } else {
    lines = 0
  }

  if (!geom_data %in% c("point", "line", FALSE))
    stop("`geom_data` has to be one of 'point', 'line', or FALSE.")

  if (!is.logical(cp_dens))
    stop("`cp_dens` must be TRUE or FALSE.")

  if (is.logical(q_fit) && all(q_fit == TRUE))
    q_fit = c(0.025, 0.975)

  if (is.logical(q_predict) && all(q_predict == TRUE))
    q_predict = c(0.025, 0.975)

  if (!is.logical(q_fit) & !is.numeric(q_fit))
    stop("`q_fit` has to be TRUE, FALSE, or a vector of numbers.")

  if (is.numeric(q_fit) & (any(q_fit > 1) | any(q_fit < 0)))
    stop ("All `q_fit` have to be between 0 (0%) and 1 (100%).")

  if (!is.logical(q_predict) & !is.numeric(q_predict))
    stop("`q_predict` has to be TRUE, FALSE, or a vector of numbers.")

  if (is.numeric(q_predict) & (any(q_predict > 1) | any(q_predict < 0)))
    stop ("All `q_predict` have to be between 0 (0%) and 1 (100%).")

  if (!is.null(nsamples)) {
    check_integer(nsamples, "nsamples", lower = 1)
    if (lines != FALSE & nsamples < lines)
      stop("`lines` must be less than or equal to `nsamples`.")
  }

  # No need for more samples if they are only used to draw lines.
  if (all(q_fit == FALSE) & all(q_predict == FALSE))
    nsamples = lines

  # Is facet_by a random/nested effect?
  if (!is.null(facet_by)) {
    if (is.character(facet_by)) {
      varying_groups = logical0_to_null(unique(stats::na.omit(fit$.other$ST$cp_group_col)))
      if (!facet_by %in% varying_groups)
        stop("`facet_by` must be a data column and modeled as a varying effect.")
    } else {
      stop("`facet_by` must be a character string. Got ", class(facet_by), ".")
    }
  }

  # Useful vars
  xvar = rlang::sym(fit$pars$x)
  yvar = rlang::sym(fit$pars$y)
  is_arma = length(fit$pars$arma) > 0


  ############################
  # MAKE NEWDATA AND PREDICT #
  ############################
  newdata = tibble::tibble(!!xvar := get_eval_at(fit, facet_by))

  if (is.null(facet_by) == TRUE) {
    varying_pars = NULL
    if (is_arma)
      newdata$.ydata = fit$data[, fit$pars$y]
  } else {
    varying_pars = unpack_varying(fit, cols = facet_by)$pars
    if (is_arma) {
      # If ARMA, replace newdata with the original data that includes xvar, yvar, and the varying effect.
      newdata = dplyr::rename(fit$data, .ydata = as.character(yvar))
    } else {
      # Else, evaluate the same x for all varying levels
      newdata = tidyr::expand_grid(newdata, !!facet_by := unique(dplyr::pull(fit$data, facet_by)))
    }
  }

  if (fit$family$family == "binomial") {
    # Interpolate trials for binomial at all xvar to make sure that there are actually values to plot
    newdata[, fit$pars$trials] = stats::approx(x = fit$data[, fit$pars$x], y = fit$data[, fit$pars$trials], xout = dplyr::pull(newdata, xvar))$y %>%
      suppressWarnings() %>%
      round()  # Only integers
  }

  # Predict
  local_pp_eval = function(type) {
    pp_eval(
      object = fit,
      newdata = newdata,
      summary = FALSE,  # Get samples
      type = type,
      rate = rate,
      prior = prior,
      which_y = which_y,
      varying = varying_pars,
      arma = arma,
      nsamples = nsamples,
      samples_format = "tidy"
    ) %>%
      dplyr::rename(!!yvar := !!type)  # from "predict"/"fitted" to yvar (response name)
  }

  # Get data with fitted values. Optionally add predictions.
  samples_expanded = local_pp_eval("fitted")
  if (any(q_predict != FALSE))
    samples_expanded$.predicted = local_pp_eval("predict") %>% dplyr::pull(yvar)


  ###########
  # PLOT IT #
  ###########
  # If this is a binomial rate, divide by the number of trials
  if (fit$family$family == "binomial" & rate == TRUE) {
    fit$data[, fit$pars$y] = fit$data[, fit$pars$y] / fit$data[, fit$pars$trials]
  }

  # If this is time series, strip fit$data$y for the "ts" class to avoid ggplot2 warning about scale picking..
  # TO DO: hack.
  fit$data[, fit$pars$y] = as.numeric(fit$data[, fit$pars$y])

  # Initiate plot and show raw data (only applicable when which_y == "ct")
  gg = ggplot(fit$data, aes_string(x = fit$pars$x, y = fit$pars$y))
  if (which_y == "ct") {
    if (geom_data == "point") {
      if (is.null(fit$pars$weights)) {
        gg = gg + geom_point()
      } else {
        gg = gg + geom_point(aes(size = fit$data[, fit$pars$weights[1]])) +
          ggplot2::scale_size_area(max_size = 2 * 1.5/sqrt(1.5))  # See https://stackoverflow.com/questions/63023877/setting-absolute-point-size-for-geom-point-with-scale-size-area/63024297?noredirect=1#comment111454629_63024297
      }
    } else if (geom_data == "line") {
      gg = gg + geom_line()
    }
  }

  # Add lines?
  if (lines > 0) {
    data_lines = tidybayes::sample_draws(samples_expanded, lines)  # Only this number of lines
    gg = gg + geom_line(aes(group = .data$.draw), data = data_lines, color = grDevices::rgb(0.5, 0.5, 0.5, 0.4))
  }

  # Add quantiles?
  if ((any(q_fit != FALSE))) {
    gg = gg + geom_quantiles(samples_expanded, q_fit, xvar, yvar, facet_by, color = "red", lty = 2, lwd = 0.7)
  }
  if (any(q_predict != FALSE)) {
    yvar_predict = rlang::sym(".predicted")
    gg = gg + geom_quantiles(samples_expanded, q_predict, xvar, yvar_predict, facet_by, color = "green4", lty = 2, lwd = 0.7)
  }

  # Add change point densities?
  if (cp_dens == TRUE & length(fit$model) > 1) {

    # The scale of the actual plot (or something close enough)
    # This is faster than limits_y = ggplot2::ggplot_build(gg)$layout$panel_params[[1]]$y.range
    if (which_y == "ct" & geom_data != FALSE) {
      limits_y = c(min(fit$data[, fit$pars$y]),
                   max(fit$data[, fit$pars$y]))
    } else if (any(q_predict != FALSE)) {
      limits_y = c(min(samples_expanded$.predicted),
                   max(samples_expanded$.predicted))
    } else if (as.character(yvar) %in% names(samples_expanded)) {
      limits_y = c(min(dplyr::pull(samples_expanded, as.character(yvar))),
                   max(dplyr::pull(samples_expanded, as.character(yvar))))
    } else {
      stop("Failed to draw change point density for this plot. Please raise an error on GitHub.")
    }

    gg = gg + geom_cp_density(fit, facet_by, limits_y) +
      ggplot2::coord_cartesian(
        ylim = c(limits_y[1], NA),  # Remove density flat line from view
        xlim = c(min(fit$data[, fit$pars$x]), max(fit$data[, fit$pars$x]))  # Very broad varying change point posteriors can expand beyond observed range. TO DO
      )
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

  gg = gg + ggplot2::theme(legend.position = "none")
  return(gg)
}



#' Density geom for `plot.mcpfit()`
#'
#' @aliases geom_cp_density
#' @keywords internal
#' @param fit An `mcpfit` object
#' @param facet_by `NULL` or a a string, like `plot.mcpfit(..., facet_by = "id")`.
#' @param include Boolean. If `TRUE` and `!is.null(facet_by)`, only return
#'   densities for the change points "affected" by `facet_by`. If `FALSE` and `!is.null(facet_by)`,
#'   return all densities except those "affected" by `facet_by`. Has no effect if `is.null(facet_by)`
#' @return A `ggplot2::stat_density` geom representing the change point densities.
geom_cp_density = function(fit, facet_by, limits_y) {
  dens_scale = 0.2  # Proportion of plot height
  dens_cut = 0.05 + 0.007  # How much to move density down. 5% is ggplot default. Move a bit further.

  # Get varying and population change point parameter names
  if (!is.null(facet_by)) {
    cp_matches_facet = fit$.other$ST$cp_group_col == facet_by  # Varies by this column
    cp_not_facet = !cp_matches_facet | is.na(cp_matches_facet)
    varying = stats::na.omit(fit$.other$ST$cp_group[cp_matches_facet])  # The rest
    population = stats::na.omit(fit$.other$ST$cp_name[cp_not_facet][-1])  # [-1] to remove cp_0
  } else {
    varying = NULL
    population = fit$.other$ST$cp_name[-1]
  }

  # Get samples in long format
  samples = tidy_samples(fit, population = population, varying = varying, absolute = TRUE) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("cp_"), names_to = "cp_name", values_to = "value")

  # Make the geom!
  ggplot2::stat_density(aes(
    x = value,
    y = ..scaled.. * diff(limits_y) * dens_scale +  # Scale to proportion of view
      limits_y[1] -  # Put on x-axis
      diff(limits_y) * dens_cut,  # Move a bit further down to remove zero-density line from view.
    group = paste0(.chain, cp_name),  # Apply scaling for each chain X cp_i combo
    color = .chain
  ),
  data = samples,
  position = "identity",
  geom = "line",
  show.legend = FALSE
  )
}


#' Return a geom_line representing the quantiles
#'
#' Called by `plot.mcpfit`.
#'
#' @aliases geom_quantiles
#' @keywords internal
#' @inheritParams plot.mcpfit
#' @inheritParams get_quantiles
#' @param ... Arguments passed to geom_line
#' @return A `ggplot2::geom_line` object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
geom_quantiles = function(samples, quantiles, xvar, yvar, facet_by, ...) {
  data_quantiles = get_quantiles(samples, quantiles, xvar, yvar, facet_by)

  # Return geom
  geom = ggplot2::geom_line(
    mapping = aes(
      y = .data$y,
      group = .data$quantile
    ),
    data = data_quantiles,
    ...
  )
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
  check_mcpfit(fit)

  if (!coda::is.mcmc.list(fit$mcmc_post) & !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

  if (!is.character(pars) | !is.character(regex_pars))
    stop("`pars` and `regex_pars` has to be string/character.")

  if (any(c("population", "varying") %in% pars) & length(pars ) > 1)
    stop("`pars` cannot be a vector that contains multiple elements AND 'population' or 'varying'.")

  if (any(c("hex", "scatter") %in% type) & (length(pars) != 2 | length(regex_pars) > 0))
    stop("`type` = 'hex' or 'scatter' takes exactly two parameters which must be provided via the `pars` argument")

  if ("combo" %in% type & length(type) > 1)
    stop("'combo' type cannot be combined with other types. Replace 'combo' with the types you want combo\'ed")

  check_integer(ncol, name = "ncol", lower = 1)

  if (!is.logical(prior))
    stop("`prior` must be either TRUE or FALSE.")

  # Get posterior/prior samples
  samples = mcmclist_samples(fit, prior = prior)

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
#' @return A vector of x-values to evaluate at.
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


#' Posterior Predictive Checks For Mcpfit Objects
#'
#' Plot posterior (default) or prior (`prior = TRUE`) predictive checks. This is convenience wrapper
#' around the `bayesplot::ppc_*()` methods.
#'
#' @aliases pp_check pp_check.mcpfit
#' @inheritParams pp_eval
#' @param type One of `bayesplot::available_ppc("grouped", invert = TRUE) %>% stringr::str_remove("ppc_")`
#' @param facet_by Name of a column in data modeled as varying effect(s).
#' @param nsamples Number of draws. Note that for summary geoms you may want to use all data,
#'   e.g., `pp_check(fit, type = "ribbon", nsamples = NULL)`.
#' @param ... Further arguments passed to `bayesplot::ppc_type(y, yrep, ...)`
#' @return A `ggplot2` object for single plots. Enriched by `patchwork` for faceted plots.
#' @seealso \code{\link{plot.mcpfit}} \code{\link{pp_eval}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' pp_check(ex_fit)
#' \donttest{
#' pp_check(ex_fit, type = "ecdf_overlay")
#' pp_check(another_fit, type = "loo_intervals", facet_by = "id")
#' }
#'
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

  # Check and recode inputs
  if (!is.null(facet_by))
    if (!is.character(facet_by) || length(facet_by) != 1)
      stop("`facet_by` must be a single character string.")

  if (is.null(newdata)) {
    y = as.numeric(fit$data[, fit$pars$y])  # strip simulated data of attributes
    varying_data = fit$data[, facet_by]
  } else {
    y = as.numeric(newdata[, fit$pars$y])  # strip simulated data of attributes
    varying_data = newdata[, facet_by]
  }

  allowed_types = stringr::str_remove(bayesplot::available_ppc(), "ppc_")
  allowed_types = allowed_types[stringr::str_detect(allowed_types, "_grouped") == FALSE]  # Grouped done mcp-side (see below)
  if (!(type %in% allowed_types))
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
    which_y = "ct",
    varying = varying,
    arma = arma,
    nsamples = nsamples,
    samples_format = "tidy"
  )

  # Return plot with or without facets
  if (is.null(facet_by)) {
    yrep = tidy_to_matrix(samples, "predict")
    plot_out = get_ppc_plot(fit, type, y, yrep, nsamples, samples$.draw, ...)  # One plot: use all of y and yrep
    return(plot_out)
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
    plot_out = patchwork::wrap_plots(all_plots) + patchwork::plot_layout(guides = "collect")
    return(plot_out)
  }
}


#' pp_check for loo statistics
#'
#' @aliases get_loo_plot_call
#' @keywords internal
#' @inheritParams pp_check
#' @param y Response vector
#' @param yrep S X N matrix of predicted responses
#' @param draws (required for loo-type plots) Indices of draws to use.
#' @param ... Arguments passed to `bayesplot::ppc_type(y, yrep, ...)`
#' @return A `ggplot2` object returned by `tidybayes::ppc_*(y, yrep, ...)`.
#' @return A string
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_ppc_plot = function(fit, type, y, yrep, nsamples, draws = NULL, ...) {
  is_loo = stringr::str_detect(type, "loo")

  if (is_loo == FALSE) {
    bayesplot_call = paste0("bayesplot::ppc_", type, "(y, yrep, ...)")
  } else if (is_loo == TRUE) {
    # Compute loo if missing
    fit = with_loo(fit, save_psis = TRUE, info = "Computing `fit$loo = loo(fit, save_psis = TRUE)`...")

    # Extract psis_object and lw
    psis_object = fit$loo$psis_object
    lw = psis_object$log_weights[unique(draws), ]
    psis_object$log_weights = lw
    attr(psis_object, "dims") = c(dim(yrep))

    # Build call (setting `samples` overwrites bayesplot defaults)
    bayesplot_call = paste0("bayesplot::ppc_", type, "(y, yrep, psis_object = psis_object, lw = lw, samples = nsamples, ...)")
  }

  plot_out = suppressWarnings(eval(parse(text = bayesplot_call)))
  return(plot_out)
}
