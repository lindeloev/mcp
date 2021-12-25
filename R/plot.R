# ABOUT: These are functions directly related to plotting
# ----------------

#' Underlies `plot()` and `plot_dpar()`
#'
#' @aliases get_plot
#' @keywords internal
#' @inheritParams pp_eval
#' @param x An \code{\link{mcpfit}} object
#' @param q_fit Whether to plot quantiles of the posterior (fitted value).
#'   * `TRUE` Add 2.5% and 97.5% quantiles. Corresponds to
#'       `q_fit = c(0.025, 0.975)`.
#'   * `FALSE` No quantiles
#'   * A vector of quantiles. For example, `quantiles = 0.5`
#'       plots the median and `quantiles = c(0.2, 0.8)` plots the 20% and 80%
#'       quantiles.
#' @param q_predict Same as `q_fit`, but for the prediction interval.
#' @param facet_by Character vector. Names of categorical data columns to split to facets.
#'   Can be varying or RHS categoricals.
#' @param color_by Character vector. Names of categorical data columns to color by.
#'   See `facet_by` for more.
#' @param lines Positive integer or `FALSE`. The number of fitted lines (draws).
#'   If there are categorical predictors, it is the number of fitted lines for
#'   each combination of `color_by` and `facet_by`. FALSE or `lines = 0` plots
#'   no lines. Note that lines always plot fitted values - not predicted.
#'   For prediction intervals, see the `q_predict` argument.
#' @param geom_data String. One of "point", "line" (good for time-series),
#'   or FALSE (do not plot).
#' @param cp_dens TRUE/FALSE. Plot posterior densities of the change point(s)?
#'   Currently does not respect `facet_by`. This will be added in the future.
#' @param ... Currently ignored.
#' @return A \pkg{ggplot2} object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_plot = function(x,
                       q_fit = FALSE,
                       q_predict = FALSE,
                       facet_by = NULL,
                       color_by = NULL,
                       lines = 25,
                       geom_data = "point",
                       cp_dens = TRUE,
                       rate = TRUE,
                       prior = FALSE,
                       which_y = "mu",
                       arma = TRUE,
                       nsamples = NULL,
                       scale = "response",
                       ...) {

  # Just for consistent naming in mcp
  fit = x

  ########################
  # ASSERTS AND RECODING #
  ########################
  assert_types(fit, "mcpfit")

  if (lines != FALSE) {
    assert_integer(lines, lower = 1, len = 1)
  } else {
    lines = 0
  }

  assert_value(geom_data, allowed = c("point", "line", FALSE))
  assert_logical(cp_dens)

  # Quantiles
  assert_types(q_fit, "logical", "numeric")
  assert_types(q_predict, "logical", "numeric")
  if (all(q_fit == TRUE))
    q_fit = c(0.025, 0.975)
  if (all(q_predict == TRUE))
    q_predict = c(0.025, 0.975)
  if (is.numeric(q_fit))
    assert_numeric(q_fit, lower = 0, upper = 1)
  if (is.numeric(q_predict))
    assert_numeric(q_predict, lower = 0, upper = 1)

  if (!is.null(nsamples)) {
    assert_integer(nsamples, lower = 1, len = 1)
    if (lines != FALSE && nsamples < lines)
      stop("`lines` must be less than or equal to `nsamples`.")
  }
  if (all(q_fit == FALSE) && all(q_predict == FALSE))
    # No need for more samples if they are only used to draw lines.
    nsamples = lines

  # Is facet_by a random/nested effect?
  assert_types(facet_by, "null", "character")
  assert_types(color_by, "null", "character")

  if (is.character(facet_by)) {
    categorical_cols = names(get_categorical_levels(fit$data))

    #varying_groups = logical0_to_null(unique(stats::na.omit(fit$.internal$ST$cp_group_col)))
    if (facet_by %notin% categorical_cols)
      stop("`facet_by` must be one of '", paste0(categorical_cols, collapse = "', '"), "', i.e., a data column that is modeled as categorical or varying effect.")
  }

  if (!coda::is.mcmc.list(fit$mcmc_post) && !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

  if (scale == "linear" && rate == FALSE)
    message("Known bug: the data points are plotted incorrectly when scale = 'linear' and rate = FALSE.")

  assert_ellipsis(...)

  # Useful vars
  xvar = rlang::sym(fit$pars$x)
  yvar = rlang::sym(fit$pars$y)
  varying_pars = unpack_varying(fit, cols = facet_by)$pars

  ############################
  # MAKE NEWDATA AND PREDICT #
  ############################
  newdata = interpolate_newdata(fit, facet_by)

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
      samples_format = "tidy",
      scale = scale
    ) %>%
      dplyr::select(-dplyr::any_of(c(as.character(yvar)))) %>%  # Only a problem for ar() models
      dplyr::rename(!!yvar := !!type)  # from "predict"/"fitted" to yvar (response name)
  }

  # Get data with fitted values. Optionally add predictions
  samples_expanded = local_pp_eval("fitted") %>%
    add_group()


  if (any(q_predict != FALSE))
    samples_expanded$.predicted = local_pp_eval("predict") %>% dplyr::pull(yvar)


  ###############################
  # PREP RESPONSE DATA FOR PLOT #
  ###############################
  ydata = fit$data[, fit$pars$y]  # Convenient shortname

  # If this is a binomial rate, divide by the number of trials
  if (fit$family$family == "binomial" && rate == TRUE)
    ydata = ydata / fit$data[, fit$pars$trials]

  # Show data
  if (scale == "linear") {
    ydata = fit$family$linkfun(ydata)
    if (any(is.infinite(ydata)))
      message("Removing points with infinite values on the linear scale. You may get a few warnings.")
    ydata[is.infinite(ydata)] = NA
  }

  # Color info
  fit$data = add_group(fit$data)
  use_color = length(unique(fit$data$.group)) > 1

  # If this is time series, strip fit$data$y for the "ts" class to avoid ggplot2 warning about scale picking..
  # TO DO: hack.
  fit$data[, fit$pars$y] = as.numeric(ydata)


  ###########
  # PLOT IT #
  ###########
  # Initiate plot and show raw data (only applicable when which_y == "mu"
  if (use_color) {
    gg = ggplot2::ggplot(fit$data, ggplot2::aes_string(x = fit$pars$x, y = fit$pars$y, color = ".group"))
  } else {
    gg = ggplot2::ggplot(fit$data, ggplot2::aes_string(x = fit$pars$x, y = fit$pars$y, color = NULL))
  }
  if (which_y == "mu") {
    if (geom_data == "point") {
      if (is.null(fit$pars$weights)) {
        gg = gg + ggplot2::geom_point()
      } else {
        gg = gg + ggplot2::geom_point(ggplot2::aes(size = fit$data[, fit$pars$weights[1]])) +
          ggplot2::scale_size_area(max_size = 2 * 1.5/sqrt(1.5))  # See https://stackoverflow.com/questions/63023877/setting-absolute-point-size-for-geom-point-with-scale-size-area/63024297?noredirect=1#comment111454629_63024297
      }
    } else if (geom_data == "line") {
      gg = gg + ggplot2::geom_line()
    }
  }

  # Add lines?
  if (lines > 0) {
    data_lines = tidybayes::sample_draws(samples_expanded %>% dplyr::group_by(.data$.group), lines) %>% dplyr::ungroup()  # Only this number of lines
    gg = gg + ggplot2::geom_line(ggplot2::aes(group = interaction(.data$.draw, .data$.group), color = .data$.group), data = data_lines, alpha = 0.4)  # color = grDevices::rgb(0.5, 0.5, 0.5, 0.4)
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
  if (cp_dens == TRUE && length(fit$model) > 1) {

    # The scale of the actual plot (or something close enough)
    # This is faster than limits_y = ggplot2::ggplot_build(gg)$layout$panel_params[[1]]$y.range
    if (which_y == "mu" && geom_data != FALSE) {
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

    gg = gg + geom_cp_density(fit, facet_by, prior, limits_y) +
      ggplot2::coord_cartesian(
        ylim = c(limits_y[1], NA),  # Remove density flat line from view
        xlim = c(min(fit$data[, fit$pars$x]), max(fit$data[, fit$pars$x]))  # Very broad varying change point posteriors can expand beyond observed range. TO DO
      )
  }

  # Add faceting?
  if (!is.null(facet_by)) {
    gg = gg + ggplot2::facet_wrap(paste0("~", facet_by))
  }

  # Add better y-labels
  if (scale == "linear")
    gg = gg + ggplot2::labs(y = paste0(fit$family$link, "(", fit$pars$y, ")"))
  if (scale == "response" && (fit$family$family == "bernoulli" || (fit$family$family == "binomial" && rate == TRUE)))
    gg = gg + ggplot2::labs(y = paste0("P(", fit$pars$y, " = TRUE)"))
  if (which_y != "mu")
    gg = gg + ggplot2::labs(y = which_y)

  # No color if no categorical predictors
  if (use_color == FALSE) {
    gg = gg +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_color_manual(values = "#858585")
  } else {
    gg = gg +
      ggplot2::scale_color_viridis_d(end = 0.9) +   # Yellow is not distinct from the background
      ggplot2::theme(legend.title = ggplot2::element_blank())
  }

  # Return
  gg
}


#' Plot full fits
#'
#' Plot prior or posterior model draws on top of data.
#'
#' @aliases plot plot.mcpfit
#' @export
#' @inheritParams get_plot
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @seealso plot_pars plot_dpar pp_check
#' @details
#'   `plot()` uses `fit$simulate()` on posterior samples. These represent the
#'   (joint) posterior distribution.fit
#' @return A \pkg{ggplot2} object.
#' @examples
#' # Typical usage. demo_fit is an mcpfit object.
#' plot(demo_fit)
#' \donttest{
#' plot(demo_fit, prior = TRUE)  # The prior
#'
#' plot(demo_fit, lines = 0, q_fit = TRUE)  # 95% HDI without lines
#' plot(demo_fit, q_predict = c(0.1, 0.9))  # 80% prediction interval
#' plot(demo_fit, which_y = "sigma", lines = 100)  # The variance parameter on y
#'
#' # Show a panel for each varying effect
#' # plot(fit, facet_by = "my_column")
#'
#' # Customize plots using regular ggplot2
#' library(ggplot2)
#' plot(demo_fit) + theme_bw(15) + ggtitle("Great plot!")
#' }
plot.mcpfit = function(x,
                    q_fit = FALSE,
                    q_predict = FALSE,
                    facet_by = NULL,
                    color_by = NULL,
                    lines = 25,
                    geom_data = "point",
                    cp_dens = TRUE,
                    rate = TRUE,
                    prior = FALSE,
                    arma = TRUE,
                    nsamples = NULL,
                    ...) {

  args = list(...)
  if ("which_y" %in% names(args))
    warning("plot(fit, which_y = dpar) was deprecated since mcp v0.4. Use plot_dpar() instead.")

  get_plot(
    x,
    q_fit = q_fit,
    q_predict = q_predict,
    facet_by = facet_by,
    color_by = color_by,
    lines = lines,
    geom_data = geom_data,
    cp_dens = cp_dens,
    rate = rate,
    prior = prior,
    arma = arma,
    nsamples = nsamples,
    scale = "response",
    ...
  )
}


#' @aliases plot_dpar
#' @describeIn plot.mcpfit Plot distributional parameters
#' @export
plot_dpar = function(x,
                     dpar = "mu",
                     q_fit = FALSE,
                     facet_by = NULL,
                     color_by = NULL,
                     lines = 25,
                     cp_dens = TRUE,
                     prior = FALSE,
                     arma = TRUE,
                     nsamples = NULL,
                     scale = "response",
                     ...) {

  get_plot(
    x,
    dpar,
    q_fit = q_fit,
    q_predict = FALSE,
    facet_by = facet_by,
    color_by = color_by,
    lines = lines,
    geom_data = FALSE,
    cp_dens = cp_dens,
    rate = TRUE,
    prior = prior,
    arma = arma,
    nsamples = nsamples,
    scale = scale,
    ...
  )

}


#' Density geom for `plot.mcpfit()`
#'
#' Note that `geom_density(fill = ...)` always fill area to y = 0 but we want to fill
#' to the bottom of the plot. So we use `geom_polygon(fill = ...)` instead.
#'
#' @aliases geom_cp_density
#' @keywords internal
#' @noRd
#' @inheritParams plot.mcpfit
#' @param limits_y A vector of length 2 with c(lower, upper) limits on the plot.
#'   Used for scaling the densities to a proportion of the plot height.
#' @return A `ggplot2::geom_polygon()` representing the change point densities.
geom_cp_density = function(fit, facet_by, prior, limits_y) {
  dens_scale = 0.2  # Proportion of plot height
  dens_cut = 0.05  # How much to move density down. 5% is ggplot default. Move a bit further.

  # facet_by will expand by group in tidy_samples(). Categorical cols share
  # parameters across facets, so only expand for varying effects.
  cp_matches_facet = fit$.internal$ST$cp_group_col == facet_by  # Varies by this column
  cp_not_facet = cp_matches_facet == FALSE | is.na(cp_matches_facet)
  if (all(cp_not_facet))
    facet_by = NULL

  # Get varying and population change point parameter names
  if (!is.null(facet_by)) {
    varying = stats::na.omit(fit$.internal$ST$cp_group[cp_matches_facet])  # The rest
    population = stats::na.omit(fit$.internal$ST$cp_name[cp_not_facet][-1])  # [-1] to remove cp_0
  } else {
    varying = NULL
    population = fit$.internal$ST$cp_name[-1]
  }

  # Get samples in long format
  samples = tidy_samples(fit, population = population, varying = varying, absolute = TRUE, prior = prior) %>%
    tidyr::pivot_longer(cols = tidyselect::matches("^cp_[0-9]+$"), names_to = "cp_name", values_to = "value") %>%

    # Compute density per group
    #dplyr::mutate(!!facet_by := dplyr::if_else(facet_by %in% colnames(.), !!facet_by, NULL)) %>%  # Attempt at faceting non-varying
    dplyr::group_by(dplyr::across(c(.chain, cp_name, !!facet_by))) %>%
    dplyr::summarise(dens = list(stats::density(value, bw = "SJ", n = 2^10))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      densx = list(.data$dens$x),
      densy = list((.data$dens$y / max(.data$dens$y)) *  # Scale to 1 height
                     dens_scale * diff(limits_y) +  # Scale to desired proportion of plot
                     limits_y[1] -  # Put on x-axis
                     diff(limits_y) * dens_cut  # Move a bit further down to remove zero-density line from view.
      )
    ) %>%
    dplyr::select(-.data$dens) %>%
    tidyr::unnest(c(.data$densx, .data$densy))


  # Make the geom!
  ggplot2::geom_polygon(ggplot2::aes(
      x = .data$densx,
      y = .data$densy,
      group = interaction(.chain, cp_name),
      color = NULL
    ),
    data = samples,
    alpha = 1 / max(samples$.chain),  # Sum to opaque
    show.legend = FALSE
  )
}



#' Return a geom_line representing the quantiles
#'
#' Called by `plot.mcpfit`.
#'
#' @aliases geom_quantiles
#' @keywords internal
#' @noRd
#' @inheritParams plot.mcpfit
#' @inheritParams get_quantiles
#' @param ... Arguments passed to geom_line
#' @return A `ggplot2::geom_line` object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
geom_quantiles = function(samples, quantiles, xvar, yvar, facet_by, ...) {
  data_quantiles = get_quantiles(samples, quantiles, xvar, yvar, facet_by)

  # Return geom
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      y = .data$y,
      group = .data$quantile
    ),
    data = data_quantiles,
    ...
  )
}
