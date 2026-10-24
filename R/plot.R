# ABOUT: These are functions directly related to plotting
# ----------------

# Add `.group` (curve identity) and `.color` columns: interactions of `curve_by` and `color_by` columns
add_plot_groups = function(df, curve_by = names(get_categorical_levels(df)), color_by = NULL) {
  curve_by = intersect(curve_by, names(df))
  color_by = intersect(color_by, names(df))

  make_group = function(cols) {
    if (length(cols) == 0)
      return(rep.int(1L, nrow(df)))

    # Evaluated draws repeat the same plotting metadata for every posterior
    # draw. Build the interaction once per evaluation row, then expand it.
    if (".row" %in% names(df)) {
      first = !duplicated(df$.row)
      group = interaction(df[first, cols, drop = FALSE], drop = TRUE, sep = ":")
      return(group[match(df$.row, df$.row[first])])
    } else {
      interaction(df[, cols, drop = FALSE], drop = TRUE, sep = ":")
    }
  }

  df$.group = make_group(curve_by)
  df$.color = make_group(color_by)
  df
}


# Underlies plot() and plot_dpar(). Arguments are documented in ?plot.mcpfit.
get_plot = function(x,
                       interval = "none",
                       width = 0.95,
                       facet_by = NULL,
                       color_by = "auto",
                       lines = 25,
                       geom_data = "point",
                       cp_dens = TRUE,
                       rate = TRUE,
                       prior = FALSE,
                       dpar = "epred",
                       arma = TRUE,
                       ndraws = 500,
                       scale = "response",
                       at = NULL) {
  # Just for consistent naming in mcp
  fit = x

  ########################
  # ASSERTS AND RECODING #
  ########################
  check_mcpfit_version(fit)
  checkmate::assert_flag(rate)
  assert_arma(arma, replicate_ok = FALSE)

  # Counts scale with auxiliary response columns such as trials(N), so they must be explicit
  epred_aux = names(Filter(function(aux) "epred" %in% aux$operations, fit$family$response$auxiliary))
  missing_cols = setdiff(unlist(mcp_columns(fit)[epred_aux]), names(at))
  if (!rate && length(missing_cols) > 0) {
    stop(
      "`plot(..., rate = FALSE)` requires explicit values for ", and_collapse(missing_cols), ". ",
      "Supply them as `at = list(", missing_cols[1], " = ...)`.",
      call. = FALSE
    )
  }

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
  # Intervals
  interval = match.arg(interval, c("none", "fitted", "predict"), several.ok = TRUE)
  assert_width(width)
  show_fitted = "fitted" %in% interval
  show_predict = "predict" %in% interval

  # Validate columns used for faceting and color.
  checkmate::assert_character(facet_by, null.ok = TRUE)
  checkmate::assert_character(color_by, null.ok = TRUE)
  facet_by = logical0_to_null(unique(facet_by))
  color_by = logical0_to_null(unique(color_by))

  # Curves are identified by categorical predictors, series, and plotted grouping columns.
  # `color_by = "auto"` colors by the only categorical predictor, if there is one.
  model_tables = get_fit_model_tables(fit)
  data_columns = mcp_columns(fit)
  all_categorical_cols = names(get_categorical_levels(fit$data))
  group_cols = unique(stats::na.omit(mcp_pars(fit, scope = "group")$group_col))
  categorical_cols = intersect(all_categorical_cols, get_predictor_cols(fit))
  check_unmapped = !is.null(color_by)  # `color_by = NULL` allows unidentified curves
  if (identical(color_by, "auto"))
    color_by = if (length(categorical_cols) == 1) logical0_to_null(setdiff(categorical_cols, facet_by)) else NULL
  plot_by = unique(c(facet_by, color_by))
  curve_by = unique(c(categorical_cols, data_columns$series, intersect(group_cols, plot_by)))
  valid_group_cols = unique(c(categorical_cols, group_cols, data_columns$series))

  validate_plot_groups = function(cols, arg) {
    invalid_cols = setdiff(cols, valid_group_cols)
    if (length(invalid_cols) > 0) {
      valid_text = if (length(valid_group_cols) == 0) {
        "There are no categorical predictor or grouping columns in this model."
      } else {
        paste0("Valid columns are '", paste(valid_group_cols, collapse = "', '"), "'.")
      }
      stop(
        "`", arg, "` must name categorical predictor or grouping columns. ",
        "Invalid: '", paste(invalid_cols, collapse = "', '"), "'. ",
        valid_text
      )
    }
  }
  validate_plot_groups(facet_by, "facet_by")
  validate_plot_groups(color_by, "color_by")

  # Intervals are ambiguous for curves that are not mapped to color or facets
  unmapped = setdiff(categorical_cols, plot_by)
  if ((show_fitted || show_predict) && length(categorical_cols) > 1 &&
      check_unmapped && length(unmapped) > 0) {
    stop(
      "Unmapped categorical predictors: '", paste(unmapped, collapse = "', '"), "'. ",
      "Use `color_by` and/or `facet_by` to identify all fitted curves, or set `color_by = NULL`."
    )
  }

  if (!coda::is.mcmc.list(.subset2(fit, "mcmc_post")) && !coda::is.mcmc.list(.subset2(fit, "mcmc_prior")))
    stop("Cannot plot an mcpfit without prior or posterior draws.")

  available_draws = sum(vapply(mcmclist_draws(fit, prior = prior), nrow, integer(1)))  # Like niterations(fit), but also supporting prior = TRUE
  if (!is.null(ndraws))
    ndraws = min(ndraws, available_draws)
  lines = min(lines, available_draws)

  # Useful vars
  y_col = data_columns$response
  group_effects = unpack_group_effects(fit, cols = plot_by)$pars
  if (is.null(group_effects))
    group_effects = FALSE

  ############################
  # MAKE NEWDATA AND PREDICT #
  ############################
  needs_arma = uses_garma(fit, arma, "fitted", dpar)
  newdata = interpolate_newdata(fit, by = plot_by, at = at, arma = needs_arma, prior = prior)

  # Evaluate like fitted() and predict(). The response is in newdata for AR/MA models, so drop it.
  local_pp_eval = function(type, summary, ndraws) {
    pp_eval(
      fit, newdata, summary = summary, type = type, width = width, rate = rate, prior = prior,
      dpar = if (type == "fitted") dpar, group_effects = group_effects, arma = needs_arma, ndraws = ndraws, scale = scale
    ) |>
      dplyr::select(-dplyr::any_of(y_col))
  }

  # Fitted lines need only `lines` draws, selected jointly across all curves. Values are in the response column.
  eval_lines = if (lines > 0) {
    local_pp_eval("fitted", summary = FALSE, ndraws = lines) |>
      dplyr::rename("{y_col}" := dplyr::all_of(get_value_col("fitted", dpar, scale))) |>
      add_plot_groups(curve_by = curve_by, color_by = color_by)
  }

  # Central intervals, with the lower and upper ends in the `quantile` column
  get_interval = function(type, value_col) {
    local_pp_eval(type, summary = TRUE, ndraws = ndraws) |>
      dplyr::select(-".mean", -".sd") |>
      tidyr::pivot_longer(c(".lower", ".upper"), names_to = "quantile", values_to = value_col) |>
      add_plot_groups(curve_by = curve_by, color_by = color_by)
  }
  fitted_data = if (show_fitted) get_interval("fitted", y_col)
  predict_data = if (show_predict) get_interval("predict", ".prediction")
  plot_newdata = add_plot_groups(newdata, curve_by = curve_by, color_by = color_by)


  ###############################
  # PREP RESPONSE DATA FOR PLOT #
  ###############################
  ydata = fit$data[, data_columns$response]  # Convenient shortname
  response_data = get_family_response_data(fit$family, model_tables$segments, fit$data)

  # Show data
  if (scale == "linear") {
    ydata = fit$family$linkfun(fit$family$response$observed(
      ydata, response_data, rate = TRUE
    ))
    n_infinite = sum(is.infinite(ydata))
    if (n_infinite > 0)
      message("Not showing ", n_infinite, " data point(s) with infinite values on the linear scale.")
    ydata[is.infinite(ydata)] = NA
  } else {
    ydata = fit$family$response$observed(ydata, response_data, rate)
  }

  # Color info
  fit$data = add_plot_groups(fit$data, curve_by = curve_by, color_by = color_by)
  ncolors = length(unique(plot_newdata$.color))
  use_color = !is.null(color_by) && ncolors > 1

  # Store the plotting-scale response and strip any "ts" class to avoid
  # ggplot2 warnings about scale selection.
  fit$data[, data_columns$response] = as.numeric(ydata)


  ###########
  # PLOT IT #
  ###########
  # Initiate plot and show raw data (only applicable in plot.mcpfit())
  if (use_color) {
    gg = ggplot2::ggplot(fit$data, ggplot2::aes(x = .data[[data_columns$par_x]], y = .data[[data_columns$response]], group = .data$.group, color = .data$.color))
  } else {
    gg = ggplot2::ggplot(fit$data, ggplot2::aes(x = .data[[data_columns$par_x]], y = .data[[data_columns$response]], group = .data$.group))
  }
  if (dpar == "epred") {
    if (geom_data == "point") {
      observed_data = fit$data[!is.na(fit$data[[data_columns$response]]), , drop = FALSE]
      weights_col = data_columns$weights
      if (is.null(weights_col)) {
        gg = gg + ggplot2::geom_point(data = observed_data)
      } else {
        gg = gg + ggplot2::geom_point(ggplot2::aes(size = .data[[weights_col]]), data = observed_data) +
          ggplot2::scale_size_area(max_size = 2 * 1.5/sqrt(1.5))  # See https://stackoverflow.com/questions/63023877/setting-absolute-point-size-for-geom-point-with-scale-size-area/63024297?noredirect=1#comment111454629_63024297
      }
    } else if (geom_data == "line") {
      gg = gg + ggplot2::geom_line()
    }
  }

  # Add lines?
  if (lines > 0) {
    line_mapping = if (use_color) {
      ggplot2::aes(group = interaction(.data$.draw, .data$.group), color = .data$.color)
    } else {
      ggplot2::aes(group = interaction(.data$.draw, .data$.group))
    }
    gg = gg + ggplot2::geom_line(line_mapping, data = eval_lines, alpha = 0.2)
  }

  # Add intervals?
  if (show_fitted) {
    gg = gg + geom_quantiles(fitted_data, y_col, use_color, color = "red", lty = "dashed", lwd = 0.85)
  }
  if (show_predict) {
    gg = gg + geom_quantiles(predict_data, ".prediction", use_color, color = "#009E73", lty = "dashed", lwd = 0.85)
  }

  # Add change point densities?
  if (cp_dens == TRUE && length(fit$model) > 1) {

    # The scale of the actual plot (or something close enough)
    # This is faster than limits_y = ggplot2::ggplot_build(gg)$layout$panel_params[[1]]$y.range
    if (dpar == "epred" && geom_data != FALSE) {
      limits_y = range(fit$data[[data_columns$response]], finite = TRUE)
    } else if (show_predict) {
      limits_y = range(predict_data$.prediction)
    } else if (show_fitted) {
      limits_y = range(dplyr::pull(fitted_data, y_col))
    } else if (y_col %in% names(eval_lines)) {
      limits_y = range(dplyr::pull(eval_lines, y_col))
    } else {
      stop_github("Failed to draw change point density for this plot.")
    }

    gg = gg + geom_cp_density(fit, facet_by, prior, limits_y, use_color) +
      ggplot2::coord_cartesian(
        ylim = c(limits_y[1], NA),  # Remove density flat line from view
        # Do not let broad group-specific change-point densities expand the observed x-range.
        xlim = c(min(fit$data[, data_columns$par_x]), max(fit$data[, data_columns$par_x]))
      )
  }

  # Add faceting?
  if (!is.null(facet_by)) {
    gg = gg + ggplot2::facet_wrap(ggplot2::vars(!!!rlang::syms(facet_by)))
  }

  # Add better y-labels
  if (scale == "linear")
    gg = gg + ggplot2::labs(y = paste0(fit$family$link, "(", data_columns$response, ")"))
  if (scale == "response" && fit$family$response$probability(rate))
    gg = gg + ggplot2::labs(y = "Success probability")
  if (dpar != "epred")
    gg = gg + ggplot2::labs(y = dpar)

  # Caption with the values at which other continuous predictors are held
  at_cols = setdiff(
    names(newdata),
    c(data_columns$par_x, data_columns$response, data_columns$series, all_categorical_cols, group_cols)
  )
  at_cols = at_cols[vapply(newdata[at_cols], dplyr::n_distinct, integer(1)) == 1]  # E.g., AR/MA plots use observed values
  if (length(at_cols) > 0) {
    at_text = paste0(at_cols, " = ", format(signif(unlist(newdata[1, at_cols]), 4), trim = TRUE))
    gg = gg + ggplot2::labs(caption = paste("Continuous predictors held at", paste(at_text, collapse = ", ")))
  }

  # No color if no categorical predictors
  if (use_color == FALSE) {
    gg = gg +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_color_manual(values = "#858585")
  } else {
    plot_colors = c(
      "#0072B2", "#D55E00", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F5C710"
    )
    color_scale = if (ncolors <= length(plot_colors)) {
      ggplot2::scale_color_manual(values = plot_colors)
    } else {
      ggplot2::scale_color_viridis_d(begin = 0.05, end = 0.75)
    }
    gg = gg + color_scale + ggplot2::labs(color = paste(color_by, collapse = ":"))
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
#' @param x An \code{\link{mcpfit}} object
#' @param interval Which central intervals to plot as dashed lines:
#'   * `"none"` (default): No intervals.
#'   * `"fitted"`: Interval of the fitted values (expected response), as in `fitted()`. Red lines.
#'   * `"predict"`: Posterior predictive interval, as in `predict()`. Green lines.
#'   * `c("fitted", "predict")`: Both. `plot_dpar()` supports only `"none"` and `"fitted"`.
#'
#'   For other quantiles, see [the plotting article](https://lindeloev.github.io/mcp/articles/plotting.html#making-your-own-plots) for computing them from
#'   `fitted(..., summary = FALSE)` or `predict(..., summary = FALSE)` and adding them as a ggplot2 layer.
#' @param width Width of the central interval(s) given in `interval` (strictly between 0 and 1).
#'   For example, `width = 0.8` plots the 10% and 90% quantiles.
#' @param facet_by Character vector. Names of categorical data columns to split to facets.
#'   Can be grouping-factor or categorical predictor columns.
#' @param color_by A character vector naming categorical predictor or grouping columns to color by.
#'   Multiple columns are combined as an interaction. Curves and quantiles remain separate for grouping columns not mapped to color.
#'   `"auto"` colors by the categorical predictor if the model has exactly one and it is not in `facet_by`. `NULL` disables color.
#' @param at Named list setting additional continuous predictors to fixed values.
#'   They default to their observed means. Family-specific response auxiliaries
#'   can be supplied as explicit scalar design values. Passed to
#'   `interpolate_newdata()`.
#' @param lines Positive integer or `FALSE`. The number of fitted lines (draws).
#'   It is the number of joint posterior draws shown for every curve. FALSE or `lines = 0` plots no lines.
#'   Note that lines always plot fitted values - not predicted.
#'   For posterior predictive intervals, use `interval = "predict"`.
#' @param geom_data String. One of "point", "line" (good for time-series),
#'   or FALSE (do not plot).
#' @param cp_dens TRUE/FALSE. Plot posterior densities of the change point(s)?
#' @param rate Logical scalar. For binomial models, plot the observed or expected
#'   success proportion (`rate = TRUE`, the default here) or counts (`rate = FALSE`).
#'   This default is the opposite of `fitted()` and `predict()`, which return
#'   counts by default; plotting proportions keeps curves comparable across
#'   varying trial counts. Predictions and count-scale fitted values require a
#'   trials column in `newdata`/`at`. Other families do not support `rate = TRUE`,
#'   so the default only applies to binomial and Bernoulli models.
#' @inheritParams fitted.mcpfit
#' @param q_fit,q_predict Deprecated. Use `interval` and `width` instead.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @param ... Must be empty. Reserved for future use.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @seealso plot_pars plot_dpar pp_check
#' @details
#'   `plot()` evaluates the model on posterior draws like `fitted()` and `predict()`. These represent the
#'   (joint) posterior distribution. Interval summaries are based on `ndraws` draws, and change-point
#'   densities are simple densities on all available draws.
#' @return A \pkg{ggplot2} object.
#' @examples
#' # Typical usage. demo_fit is an mcpfit object.
#' plot(demo_fit)
#' \donttest{
#' plot(demo_fit, prior = TRUE)  # The prior
#'
#' plot(demo_fit, lines = 0, interval = "fitted")  # 95% central interval without lines
#' plot(demo_fit, interval = "predict", width = 0.8)  # 80% posterior predictive interval
#' plot_dpar(demo_fit, dpar = "sigma", lines = 100)  # Residual standard deviation on y
#'
#' # Show a panel for each group-level effect
#' # plot(fit, facet_by = "my_column")
#'
#' # Customize plots using regular ggplot2
#' library(ggplot2)
#' plot(demo_fit) + theme_bw(15) + ggtitle("Great plot!")
#' }
plot.mcpfit = function(x,
                    interval = "none",
                    width = 0.95,
                    facet_by = NULL,
                    color_by = "auto",
                    lines = 25,
                    geom_data = "point",
                    cp_dens = TRUE,
                    rate = TRUE,
                    prior = FALSE,
                    arma = TRUE,
                    ndraws = 500,
                    scale = "response",
                    at = NULL,
                    q_fit = lifecycle::deprecated(),
                    q_predict = lifecycle::deprecated(),
                    nsamples = lifecycle::deprecated(),
                    ...) {
  ndraws = legacy_arg(ndraws, nsamples, missing(ndraws), "plot.mcpfit", "ndraws", "nsamples")
  legacy = legacy_plot_interval(interval, width, q_fit, q_predict, missing(interval) && missing(width), "plot.mcpfit")
  which_y = legacy_which_y(list(...), "plot")
  if (missing(rate))
    rate = x$family$response$has_rate

  get_plot(
    x,
    interval = legacy$interval,
    width = legacy$width,
    facet_by = facet_by,
    color_by = color_by,
    at = at,
    lines = lines,
    geom_data = geom_data,
    cp_dens = cp_dens,
    rate = rate,
    prior = prior,
    dpar = if (is.null(which_y)) "epred" else which_y,
    arma = arma,
    ndraws = ndraws,
    scale = scale
  )
}


#' @aliases plot_dpar
#' @describeIn plot.mcpfit Plot distributional parameters
#' @export
plot_dpar = function(x,
                     dpar = "epred",
                     interval = "none",
                     width = 0.95,
                     facet_by = NULL,
                     color_by = "auto",
                     lines = 25,
                     cp_dens = TRUE,
                     prior = FALSE,
                     arma = TRUE,
                     ndraws = 500,
                     scale = "response",
                     at = NULL,
                     ...) {
  rlang::check_dots_empty()
  interval = match.arg(interval, c("none", "fitted"))

  get_plot(
    x,
    interval = interval,
    width = width,
    facet_by = facet_by,
    color_by = color_by,
    at = at,
    lines = lines,
    geom_data = FALSE,
    cp_dens = cp_dens,
    rate = x$family$response$has_rate,
    prior = prior,
    dpar = dpar,
    arma = arma,
    ndraws = ndraws,
    scale = scale
  )
}


#' Density geom for `plot.mcpfit()`
#'
#' Note that `geom_density(fill = ...)` always fill area to y = 0 but we want to fill
#' to the bottom of the plot. So we use `geom_polygon(fill = ...)` instead.
#'
#' @inheritParams plot.mcpfit
#' @param limits_y A vector of length 2 with c(lower, upper) limits on the plot.
#'   Used for scaling the densities to a proportion of the plot height.
#' @return A `ggplot2::geom_polygon()` representing the change point densities.
#' @noRd
geom_cp_density = function(fit, facet_by, prior, limits_y, use_color = FALSE) {
  dens_scale = 0.2  # Proportion of plot height
  dens_cut = 0.05  # How much to move density down. 5% is ggplot default. Move a bit further.

  # Faceting by the grouping column of group-level change points shows each group's change
  # points, cp_1 + cp_1_id[level]. Categorical cols share parameters across facets.
  cps = get_fit_model_tables(fit)$cps
  cp_matches_facet = cps$group_col %in% facet_by
  facet_by = unique(cps$group_col[cp_matches_facet])
  if (length(facet_by) == 0)
    facet_by = NULL

  # Change point draws in long format: one row per draw, change point, and facet level
  draws = posterior::subset_draws(posterior_draws(fit, prior = prior), variable = c(cps$name, cps$group_name[cp_matches_facet]))
  draws = posterior::as_draws_df(draws)
  levels = if (is.null(facet_by)) NA else unique(fit$data[[facet_by]])
  draws = dplyr::bind_rows(lapply(seq_along(levels), function(level_index) dplyr::bind_rows(lapply(seq_len(nrow(cps)), function(j) {
    value = draws[[cps$name[j]]]
    if (cp_matches_facet[j])
      value = value + draws[[paste0(cps$group_name[j], "[", levels[level_index], "]")]]
    tibble::tibble(.chain = draws$.chain, cp_name = cps$name[j], level_index = level_index, value = value)
  }))))
  if (!is.null(facet_by))
    draws[[facet_by]] = levels[draws$level_index]

  # Compute density per group. Tolerate constant change points like cp_2 = 80.
  draws = draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(".chain", "cp_name", facet_by)))) |>
    dplyr::summarise(dens = list(
      if (is.na(stats::sd(.data$value)) || stats::sd(.data$value) == 0) {
        list(x = rep(.data$value[1], 2), y = c(0, 1))
      } else {
        tryCatch(
          stats::density(.data$value, bw = "SJ", n = 2^10),
          error = function(e) stats::density(.data$value, bw = "nrd0", n = 2^10)
        )
      }
    )) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      densx = list(.data$dens$x),
      densy = list((.data$dens$y / max(.data$dens$y)) *  # Scale to 1 height
                     dens_scale * diff(limits_y) +  # Scale to desired proportion of plot
                     limits_y[1] -  # Put on x-axis
                     diff(limits_y) * dens_cut  # Move a bit further down to remove zero-density line from view.
      )
    ) |>
    dplyr::select(-"dens") |>
    tidyr::unnest(c("densx", "densy"))


  # Make the geom!
  ggplot2::geom_polygon(ggplot2::aes(
      x = .data$densx,
      y = .data$densy,
      group = interaction(.data$.chain, .data$cp_name)
    ),
    data = draws,
    color = if (use_color) "#666666" else "#3182BD",
    fill = if (use_color) "#A6A6A6" else "#6BAED6",
    linewidth = 0.25,
    alpha = 0.4 / max(draws$.chain),  # Combined opacity is approximately 35-40%
    show.legend = FALSE
  )
}



#' Return a geom_line representing the quantiles
#'
#' Called by `plot.mcpfit`.
#'
#' @inheritParams plot.mcpfit
#' @param color Fixed line color when `use_color = FALSE`.
#' @param ... Arguments passed to geom_line
#' @return A `ggplot2::geom_line` object.
#' @noRd
geom_quantiles = function(data_quantiles, y_col, use_color = TRUE, color = NULL, ...) {
  if (use_color) {
    mapping = ggplot2::aes(y = .data[[y_col]], group = interaction(.data$quantile, .data$.group), color = .data$.color)
    ggplot2::geom_line(mapping = mapping, data = data_quantiles, ...)
  } else {
    mapping = ggplot2::aes(y = .data[[y_col]], group = interaction(.data$quantile, .data$.group))
    ggplot2::geom_line(mapping = mapping, data = data_quantiles, color = color, ...)
  }
}
