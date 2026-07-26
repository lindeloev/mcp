# ABOUT: These are functions directly related to plotting
# ----------------

# Add separate columns for curve identity and color mapping.
add_plot_groups = function(df, curve_by = names(get_categorical_levels(df)), color_by = NULL) {
  curve_by = intersect(curve_by, names(df))
  color_by = intersect(color_by, names(df))

  if (length(curve_by) == 0) {
    df$.group = factor(rep(1, nrow(df)))
  } else {
    df$.group = interaction(df[, curve_by, drop = FALSE], drop = TRUE, sep = ":")
  }

  if (length(color_by) == 0) {
    df$.color = factor(rep(1, nrow(df)))
  } else {
    df$.color = interaction(df[, color_by, drop = FALSE], drop = TRUE, sep = ":")
  }

  df
}


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
#' @param color_by `NULL` for no color grouping, or a character vector naming categorical or varying-effect data columns to color by.
#'   Multiple columns are combined as an interaction. Curves and quantiles remain separate for grouping columns not mapped to color.
#' @param lines Positive integer or `FALSE`. The number of fitted lines (draws).
#'   It is the number of joint posterior draws shown for every curve. FALSE or `lines = 0` plots no lines.
#'   Note that lines always plot fitted values - not predicted.
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
                       dpar = "epred",
                       arma = TRUE,
                       ndraws = NULL,
                       scale = "response",
                       nsamples = lifecycle::deprecated(),
                       ...) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "get_plot")

  # Just for consistent naming in mcp
  fit = x

  ########################
  # ASSERTS AND RECODING #
  ########################
  checkmate::assert_class(fit, "mcpfit")
  if (!is.mcpfamily(fit$family))
    fit$family = mcpfamily(fit$family)

  if (lines != FALSE) {
    checkmate::assert_int(lines, lower = 1)
  } else {
    lines = 0
  }

  checkmate::assert(
    checkmate::check_choice(geom_data, c("point", "line")),
    checkmate::check_false(geom_data),
    .var.name = "geom_data"
  )
  checkmate::assert_flag(cp_dens)
  assert_typescale(type = "fitted", scale = scale)
  dpar = assert_dpar(dpar, fit = fit, type = "fitted")
  # Quantiles
  checkmate::assert(
    checkmate::check_flag(q_fit),
    checkmate::check_numeric(q_fit, any.missing = FALSE),
    .var.name = "q_fit"
  )
  checkmate::assert(
    checkmate::check_flag(q_predict),
    checkmate::check_numeric(q_predict, any.missing = FALSE),
    .var.name = "q_predict"
  )
  if (all(q_fit == TRUE))
    q_fit = c(0.025, 0.975)
  if (all(q_predict == TRUE))
    q_predict = c(0.025, 0.975)
  if (is.numeric(q_fit))
    checkmate::assert_numeric(q_fit, lower = 0, upper = 1, any.missing = FALSE)
  if (is.numeric(q_predict))
    checkmate::assert_numeric(q_predict, lower = 0, upper = 1, any.missing = FALSE)

  if (!is.null(ndraws)) {
    checkmate::assert_int(ndraws, lower = 1)
    if (lines != FALSE && ndraws < lines)
      stop("`lines` must be less than or equal to `ndraws`.")
  }
  if (all(q_fit == FALSE) && all(q_predict == FALSE))
    # No need for more samples if they are only used to draw lines.
    ndraws = lines

  # Validate columns used for faceting and color.
  checkmate::assert_character(facet_by, null.ok = TRUE)
  checkmate::assert_character(color_by, null.ok = TRUE)
  facet_by = logical0_to_null(unique(facet_by))
  color_by = logical0_to_null(unique(color_by))

  categorical_cols = names(get_categorical_levels(fit$data))
  varying_cols = unique(stats::na.omit(fit$.internal$ST$cp_group_col))
  curve_by = unique(c(categorical_cols, varying_cols))

  validate_plot_groups = function(cols, arg) {
    invalid_cols = setdiff(cols, curve_by)
    if (length(invalid_cols) > 0) {
      valid_text = if (length(curve_by) == 0) {
        "There are no categorical or varying-effect columns in this model."
      } else {
        paste0("Valid columns are '", paste(curve_by, collapse = "', '"), "'.")
      }
      stop(
        "`", arg, "` must name categorical or varying-effect data columns. ",
        "Invalid: '", paste(invalid_cols, collapse = "', '"), "'. ",
        valid_text
      )
    }
  }
  validate_plot_groups(facet_by, "facet_by")
  validate_plot_groups(color_by, "color_by")

  if (!coda::is.mcmc.list(fit$mcmc_post) && !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

  if (scale == "linear" && rate == FALSE)
    message("Known bug: the data points are plotted incorrectly when scale = 'linear' and rate = FALSE.")

  rlang::check_dots_empty()

  # Useful vars
  xvar = rlang::sym(fit$pars$x)
  yvar = rlang::sym(fit$pars$y)
  by = unique(c(facet_by, color_by))
  varying_pars = unpack_varying(fit, cols = by)$pars

  ############################
  # MAKE NEWDATA AND PREDICT #
  ############################
  newdata = interpolate_newdata(fit, by = by)

  # Predict
  local_pp_eval = function(type, include_fitted = FALSE) {
    pp_eval(
      object = fit,
      newdata = newdata,
      summary = FALSE,  # Get samples
      type = type,
      rate = rate,
      prior = prior,
      dpar = dpar,
      varying = varying_pars,
      arma = arma,
      ndraws = ndraws,
      samples_format = "tidy",
      scale = scale,
      .include_fitted = include_fitted
    )
  }

  # Evaluate once so fitted values and predictions use exactly the same joint
  # draw IDs. Prediction draws already contain everything needed to compute
  # their corresponding fitted values.
  if (any(q_predict != FALSE)) {
    eval_draws = local_pp_eval("predict", include_fitted = TRUE)
    eval_draws = dplyr::rename(eval_draws, .predicted = "predict")
  } else {
    eval_draws = local_pp_eval("fitted")
  }

  eval_draws = eval_draws %>%
    # Only a problem for AR/MA models, where newdata contains the response.
    dplyr::select(-dplyr::any_of(as.character(yvar))) %>%
    dplyr::rename("{fit$pars$y}" := "fitted") %>%
    add_plot_groups(curve_by = curve_by, color_by = color_by)


  ###############################
  # PREP RESPONSE DATA FOR PLOT #
  ###############################
  ydata = fit$data[, fit$pars$y]  # Convenient shortname
  response_data = get_family_response_data(fit$family, fit$.internal$ST, fit$data)
  ydata = fit$family$response$observed(ydata, response_data, rate)

  # Show data
  if (scale == "linear") {
    ydata = fit$family$linkfun(ydata)
    if (any(is.infinite(ydata)))
      message("Removing points with infinite values on the linear scale. You may get a few warnings.")
    ydata[is.infinite(ydata)] = NA
  }

  # Color info
  fit$data = add_plot_groups(fit$data, curve_by = curve_by, color_by = color_by)
  use_color = !is.null(color_by) && length(unique(fit$data$.color)) > 1

  # Store the plotting-scale response and strip any "ts" class to avoid
  # ggplot2 warnings about scale selection.
  fit$data[, fit$pars$y] = as.numeric(ydata)


  ###########
  # PLOT IT #
  ###########
  # Initiate plot and show raw data (only applicable in plot.mcpfit())
  if (use_color) {
    gg = ggplot2::ggplot(fit$data, ggplot2::aes(x = .data[[fit$pars$x]], y = .data[[fit$pars$y]], group = .data$.group, color = .data$.color))
  } else {
    gg = ggplot2::ggplot(fit$data, ggplot2::aes(x = .data[[fit$pars$x]], y = .data[[fit$pars$y]], group = .data$.group))
  }
  if (dpar == "epred") {
    if (geom_data == "point") {
      point_size = fit$family$response$point_size
      point_size_col = get_family_aux_columns(fit$family, fit$.internal$ST)[point_size]
      if (length(point_size_col) == 0 || is.na(point_size_col)) {
        gg = gg + ggplot2::geom_point()
      } else {
        gg = gg + ggplot2::geom_point(ggplot2::aes(size = fit$data[, point_size_col])) +
          ggplot2::scale_size_area(max_size = 2 * 1.5/sqrt(1.5))  # See https://stackoverflow.com/questions/63023877/setting-absolute-point-size-for-geom-point-with-scale-size-area/63024297?noredirect=1#comment111454629_63024297
      }
    } else if (geom_data == "line") {
      gg = gg + ggplot2::geom_line()
    }
  }

  # Add lines?
  if (lines > 0) {
    # Select joint posterior draws once, then show the same draws for all curves.
    selected_draws = sample(unique(eval_draws$.draw), lines)
    data_lines = dplyr::filter(eval_draws, .data$.draw %in% selected_draws)
    line_mapping = if (use_color) {
      ggplot2::aes(group = interaction(.data$.draw, .data$.group), color = .data$.color)
    } else {
      ggplot2::aes(group = interaction(.data$.draw, .data$.group))
    }
    gg = gg + ggplot2::geom_line(line_mapping, data = data_lines, alpha = 0.4)
  }

  # Add quantiles?
  if ((any(q_fit != FALSE))) {
    gg = gg + geom_quantiles(eval_draws, q_fit, xvar, yvar, facet_by, use_color = use_color, lty = 2, lwd = 0.7)
  }
  if (any(q_predict != FALSE)) {
    yvar_predict = rlang::sym(".predicted")
    gg = gg + geom_quantiles(eval_draws, q_predict, xvar, yvar_predict, facet_by, use_color = use_color, lty = 3, lwd = 0.7)
  }

  # Add change point densities?
  if (cp_dens == TRUE && length(fit$model) > 1) {

    # The scale of the actual plot (or something close enough)
    # This is faster than limits_y = ggplot2::ggplot_build(gg)$layout$panel_params[[1]]$y.range
    if (dpar == "epred" && geom_data != FALSE) {
      limits_y = c(min(fit$data[, fit$pars$y]),
                   max(fit$data[, fit$pars$y]))
    } else if (any(q_predict != FALSE)) {
      limits_y = c(min(eval_draws$.predicted),
                   max(eval_draws$.predicted))
    } else if (as.character(yvar) %in% names(eval_draws)) {
      limits_y = c(min(dplyr::pull(eval_draws, as.character(yvar))),
                   max(dplyr::pull(eval_draws, as.character(yvar))))
    } else {
      stop("Failed to draw change point density for this plot. Please raise an error on GitHub.")
    }

    gg = gg + geom_cp_density(fit, facet_by, prior, limits_y) +
      ggplot2::coord_cartesian(
        ylim = c(limits_y[1], NA),  # Remove density flat line from view
        # Do not let broad varying change-point densities expand the observed x-range.
        xlim = c(min(fit$data[, fit$pars$x]), max(fit$data[, fit$pars$x]))
      )
  }

  # Add faceting?
  if (!is.null(facet_by)) {
    gg = gg + ggplot2::facet_wrap(ggplot2::vars(!!!rlang::syms(facet_by)))
  }

  # Add better y-labels
  if (scale == "linear")
    gg = gg + ggplot2::labs(y = paste0(fit$family$link, "(", fit$pars$y, ")"))
  if (scale == "response" && fit$family$response$probability(rate))
    gg = gg + ggplot2::labs(y = paste0("P(", fit$pars$y, " = TRUE)"))
  if (dpar != "epred")
    gg = gg + ggplot2::labs(y = dpar)

  # No color if no categorical predictors
  if (use_color == FALSE) {
    gg = gg +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_color_manual(values = "#858585")
  } else {
    gg = gg +
      ggplot2::scale_color_viridis_d(end = 0.9) +   # Yellow is not distinct from the background
      ggplot2::labs(color = paste(color_by, collapse = ":"))
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
#' plot(demo_fit, lines = 0, q_fit = TRUE)  # 95% central interval without lines
#' plot(demo_fit, q_predict = c(0.1, 0.9))  # 80% prediction interval
#' plot_dpar(demo_fit, dpar = "sigma", lines = 100)  # The variance parameter on y
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
                    ndraws = NULL,
                    nsamples = lifecycle::deprecated(),
                    ...) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "plot.mcpfit")

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
    dpar = NULL,
    arma = arma,
    ndraws = ndraws,
    scale = "response",
    ...
  )
}


#' @aliases plot_dpar
#' @describeIn plot.mcpfit Plot distributional parameters
#' @export
plot_dpar = function(x,
                     dpar = "epred",
                     q_fit = FALSE,
                     facet_by = NULL,
                     color_by = NULL,
                     lines = 25,
                     cp_dens = TRUE,
                     prior = FALSE,
                     arma = TRUE,
                     ndraws = NULL,
                     scale = "response",
                     nsamples = lifecycle::deprecated(),
                     ...) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "plot_dpar")
  get_plot(
    x,
    q_fit = q_fit,
    q_predict = FALSE,
    facet_by = facet_by,
    color_by = color_by,
    lines = lines,
    geom_data = FALSE,
    cp_dens = cp_dens,
    rate = TRUE,
    prior = prior,
    dpar = dpar,
    arma = arma,
    ndraws = ndraws,
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
  varying_cols = unique(stats::na.omit(fit$.internal$ST$cp_group_col))
  facet_by = intersect(facet_by, varying_cols)
  cp_matches_facet = fit$.internal$ST$cp_group_col %in% facet_by
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
    dplyr::group_by(dplyr::across(dplyr::all_of(c(".chain", "cp_name", facet_by)))) %>%
    dplyr::summarise(dens = list(stats::density(.data$value, bw = "SJ", n = 2^10))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      densx = list(.data$dens$x),
      densy = list((.data$dens$y / max(.data$dens$y)) *  # Scale to 1 height
                     dens_scale * diff(limits_y) +  # Scale to desired proportion of plot
                     limits_y[1] -  # Put on x-axis
                     diff(limits_y) * dens_cut  # Move a bit further down to remove zero-density line from view.
      )
    ) %>%
    dplyr::select(-"dens") %>%
    tidyr::unnest(c("densx", "densy"))


  # Make the geom!
  ggplot2::geom_polygon(ggplot2::aes(
      x = .data$densx,
      y = .data$densy,
      group = interaction(.data$.chain, .data$cp_name),
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
geom_quantiles = function(samples, quantiles, xvar, yvar, facet_by, use_color = TRUE, ...) {
  # Posterior summaries are defined only by data_row. Curve, color, facet, and
  # x metadata are rejoined afterward solely to control the plotted geometry.
  keep = unique(c(as.character(xvar), facet_by, intersect(c(".group", ".color"), colnames(samples))))
  data_quantiles = get_quantiles(samples, quantiles, type = as.character(yvar), keep = keep)

  quantile_mapping = if (use_color) {
    ggplot2::aes(y = .data[[as.character(yvar)]], group = interaction(.data$quantile, .data$.group), color = .data$.color)
  } else {
    ggplot2::aes(y = .data[[as.character(yvar)]], group = interaction(.data$quantile, .data$.group))
  }

  # Return geom
  ggplot2::geom_line(mapping = quantile_mapping, data = data_quantiles, ...)
}
