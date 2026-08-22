# ABOUT: These are functions directly related to plotting
# ----------------

# Add separate columns for curve identity and color mapping.
add_plot_groups = function(df, curve_by = names(get_categorical_levels(df)), color_by = NULL) {
  curve_by = intersect(curve_by, names(df))
  color_by = intersect(color_by, names(df))

  make_group = function(cols) {
    if (length(cols) == 0)
      return(rep.int(1L, nrow(df)))

    # Evaluated draws repeat the same plotting metadata for every posterior
    # draw. Build the interaction once per evaluation row, then expand it.
    if ("data_row" %in% names(df)) {
      first = !duplicated(df$data_row)
      group = interaction(df[first, cols, drop = FALSE], drop = TRUE, sep = ":")
      return(group[match(df$data_row, df$data_row[first])])
    } else {
      interaction(df[, cols, drop = FALSE], drop = TRUE, sep = ":")
    }
  }

  df$.group = make_group(curve_by)
  df$.color = make_group(color_by)
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
#'   * A vector of quantiles. For example, `q_fit = 0.5`
#'       plots the median and `q_fit = c(0.2, 0.8)` plots the 20% and 80%
#'       quantiles.
#' @param q_predict Same as `q_fit`, but for the posterior predictive interval.
#' @param facet_by Character vector. Names of categorical data columns to split to facets.
#'   Can be grouping-factor or categorical predictor columns.
#' @param color_by A character vector naming categorical predictor or grouping columns to color by.
#'   If both `color_by` and `facet_by` are omitted, a sole categorical predictor is colored automatically. Set `color_by = NULL` explicitly to disable this.
#'   Multiple columns are combined as an interaction. Curves and quantiles remain separate for grouping columns not mapped to color.
#' @param at Named list setting additional continuous predictors to fixed values.
#'   They default to their observed means. Family-specific response auxiliaries
#'   can be supplied as explicit scalar design values. Passed to
#'   `interpolate_newdata()`.
#' @param .grouping Internal. Whether grouping arguments were omitted, mapped, or explicitly disabled.
#' @param lines Positive integer or `FALSE`. The number of fitted lines (draws).
#'   It is the number of joint posterior draws shown for every curve. FALSE or `lines = 0` plots no lines.
#'   Note that lines always plot fitted values - not predicted.
#'   For posterior predictive intervals, see the `q_predict` argument.
#' @param geom_data String. One of "point", "line" (good for time-series),
#'   or FALSE (do not plot).
#' @param cp_dens TRUE/FALSE. Plot posterior densities of the change point(s)?
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
                       ndraws = 1000,
                       scale = "response",
                       at = NULL,
                       .grouping = "auto",
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
  checkmate::assert_flag(rate)

  trials_col = mcp_columns(fit)$trials
  if (!rate && fit$family$family == "binomial" &&
      (is.null(at) || trials_col %notin% names(at))) {
    stop(
      "`plot(..., rate = FALSE)` requires an explicit binomial trial design. ",
      "Supply it as `at = list(", trials_col, " = ...)`.",
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
  show_q_fit = !isFALSE(q_fit)
  show_q_predict = !isFALSE(q_predict)

  # Validate columns used for faceting and color.
  checkmate::assert_character(facet_by, null.ok = TRUE)
  checkmate::assert_character(color_by, null.ok = TRUE)
  facet_by = logical0_to_null(unique(facet_by))
  color_by = logical0_to_null(unique(color_by))

  model_tables = get_fit_model_tables(fit)
  data_columns = mcp_columns(fit)
  all_categorical_cols = names(get_categorical_levels(fit$data))
  group_cols = unique(stats::na.omit(model_tables$group_effects$group_col))
  categorical_cols = setdiff(all_categorical_cols, c(group_cols, data_columns$series))
  plot_by = unique(c(facet_by, color_by))
  if (.grouping == "auto" && length(categorical_cols) == 1) {
    color_by = categorical_cols
    plot_by = unique(c(facet_by, color_by))
  }
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

  unmapped = setdiff(categorical_cols, plot_by)
  if ((show_q_fit || show_q_predict) && length(categorical_cols) > 1 &&
      .grouping != "none" && length(unmapped) > 0) {
    stop(
      "Unmapped categorical predictors: '", paste(unmapped, collapse = "', '"), "'. ",
      "Use `color_by` and/or `facet_by` to identify all fitted curves, or set `color_by = NULL` explicitly."
    )
  }

  if (!coda::is.mcmc.list(.subset2(fit, "mcmc_post")) && !coda::is.mcmc.list(.subset2(fit, "mcmc_prior")))
    stop("Cannot plot an mcpfit without prior or posterior draws.")

  available_draws = sum(vapply(mcmclist_draws(fit, prior = prior), nrow, integer(1)))  # Like niterations(fit), but also supporting prior = TRUE
  if (!is.null(ndraws))
    ndraws = min(ndraws, available_draws)
  lines = min(lines, available_draws)
  if (!is.null(ndraws)) {
    checkmate::assert_int(ndraws, lower = 1)
    if (lines != FALSE && ndraws < lines)
      stop("`lines` must be less than or equal to `ndraws`.")
  }
  if (!show_q_fit && !show_q_predict)
    # No need for more draws if they are only used to draw lines.
    ndraws = lines

  rlang::check_dots_empty()

  # Useful vars
  xvar = rlang::sym(data_columns$par_x)
  yvar = rlang::sym(data_columns$response)
  by = plot_by
  group_pars = unpack_group_effects(fit, cols = by)$pars

  ############################
  # MAKE NEWDATA AND PREDICT #
  ############################
  newdata = interpolate_newdata(fit, by = by, at = at)

  # Predict
  local_pp_eval = function(type, newdata, ndraws, include_fitted = FALSE) {
    pp_eval(
      object = fit,
      newdata = newdata,
      summary = FALSE,  # Get draws
      type = type,
      rate = rate,
      prior = prior,
      dpar = dpar,
      varying = group_pars,
      arma = arma,
      ndraws = ndraws,
      draws_format = "tidy",
      scale = scale,
      .include_fitted = include_fitted
    )
  }

  prepare_draws = function(draws) {
    draws = draws %>%
      # Only a problem for AR/MA models, where newdata contains the response.
      dplyr::select(-dplyr::any_of(as.character(yvar)))

    # ".epred" is only present for `type = "fitted"` (or "predict" with
    # `.include_fitted = TRUE`). Predict-only draws (no q_fit requested) have
    # no fitted values to rename.
    if (".epred" %in% names(draws))
      draws = dplyr::rename(draws, "{data_columns$response}" := ".epred")

    add_plot_groups(draws, curve_by = curve_by, color_by = color_by)
  }

  # Fitted lines need only `lines` draws, selected jointly across all curves.
  eval_lines = if (lines > 0) {
    prepare_draws(local_pp_eval("fitted", newdata, lines))
  } else {
    NULL
  }

  # Keep peak memory bounded for intervals by evaluating and immediately
  # summarising one curve at a time. Predictions also return their matching
  # fitted values, so q_fit and q_predict share the model evaluation.
  plot_newdata = add_plot_groups(newdata, curve_by = curve_by, color_by = color_by)
  interval_data = list()
  if (show_q_fit || show_q_predict) {
    curve_rows = split(seq_len(nrow(newdata)), plot_newdata$.group)
    interval_data = lapply(curve_rows, function(rows) {
      type = if (show_q_predict) "predict" else "fitted"
      draws = local_pp_eval(
        type, newdata[rows, , drop = FALSE], ndraws,
        show_q_fit && type == "predict"
      )
      if (type == "predict")
        draws = dplyr::rename(draws, .predicted = ".prediction")
      draws = prepare_draws(draws)
      keep = unique(c(as.character(xvar), facet_by, ".group", ".color"))

      list(
        fitted = if (show_q_fit) get_quantiles(draws, q_fit, as.character(yvar), keep) else NULL,
        predicted = if (show_q_predict) get_quantiles(draws, q_predict, ".predicted", keep) else NULL
      )
    })
  }
  q_fit_data = dplyr::bind_rows(lapply(interval_data, `[[`, "fitted"))
  q_predict_data = dplyr::bind_rows(lapply(interval_data, `[[`, "predicted"))


  ###############################
  # PREP RESPONSE DATA FOR PLOT #
  ###############################
  data_columns = mcp_columns(fit)
  ydata = fit$data[, data_columns$response]  # Convenient shortname
  response_data = get_family_response_data(fit$family, model_tables$segments, fit$data)
  ydata = fit$family$response$observed(ydata, response_data, rate)

  # Show data
  if (scale == "linear") {
    ydata = fit$family$linkfun(fit$family$response$observed(
      ydata, response_data, rate = TRUE
    ))
    if (any(is.infinite(ydata)))
      message("Removing points with infinite values on the linear scale. You may get a few warnings.")
    ydata[is.infinite(ydata)] = NA
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
      point_size = fit$family$response$point_size
      point_size_col = get_family_aux_columns(fit$family, model_tables$segments)[point_size]
      if (length(point_size_col) == 0 || is.na(point_size_col)) {
        gg = gg + ggplot2::geom_point(data = observed_data)
      } else {
        gg = gg + ggplot2::geom_point(ggplot2::aes(size = .data[[point_size_col]]), data = observed_data) +
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

  # Add quantiles?
  if (show_q_fit) {
    gg = gg + geom_quantiles(q_fit_data, xvar, yvar, use_color, color = "red", lty = "dashed", lwd = 0.85)
  }
  if (show_q_predict) {
    yvar_predict = rlang::sym(".predicted")
    gg = gg + geom_quantiles(q_predict_data, xvar, yvar_predict, use_color, color = "#009E73", lty = "dashed", lwd = 0.85)
  }

  # Add change point densities?
  if (cp_dens == TRUE && length(fit$model) > 1) {

    # The scale of the actual plot (or something close enough)
    # This is faster than limits_y = ggplot2::ggplot_build(gg)$layout$panel_params[[1]]$y.range
    if (dpar == "epred" && geom_data != FALSE) {
      limits_y = range(fit$data[[data_columns$response]], finite = TRUE)
    } else if (show_q_predict) {
      limits_y = range(q_predict_data$.predicted)
    } else if (show_q_fit) {
      limits_y = range(dplyr::pull(q_fit_data, as.character(yvar)))
    } else if (as.character(yvar) %in% names(eval_lines)) {
      limits_y = range(dplyr::pull(eval_lines, as.character(yvar)))
    } else {
      stop("Failed to draw change point density for this plot. Please raise an error on GitHub.")
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

  at_cols = setdiff(
    names(newdata),
    c(data_columns$par_x, data_columns$response, data_columns$series, all_categorical_cols, group_cols)
  )
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
#' @inheritParams get_plot
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @seealso plot_pars plot_dpar pp_check
#' @details
#'   `plot()` uses `fit$simulate()` on posterior draws. These represent the
#'   (joint) posterior distribution. Interval summaries use at most 1000 draws
#'   by default; use `ndraws = NULL` to use all draws. Change-point densities
#'   always use all available draws.
#' @return A \pkg{ggplot2} object.
#' @examples
#' # Typical usage. demo_fit is an mcpfit object.
#' plot(demo_fit)
#' \donttest{
#' plot(demo_fit, prior = TRUE)  # The prior
#'
#' plot(demo_fit, lines = 0, q_fit = TRUE)  # 95% central interval without lines
#' plot(demo_fit, q_predict = c(0.1, 0.9))  # 80% posterior predictive interval
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
                    ndraws = 1000,
                    at = NULL,
                    nsamples = lifecycle::deprecated(),
                    ...) {
  grouping = if (missing(color_by) && missing(facet_by)) "auto" else "mapped"
  if (!missing(color_by) && is.null(color_by)) grouping = "none"
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "plot.mcpfit")

  args = list(...)
  if ("which_y" %in% names(args)) {
    warn_which_y(args, "plot")
    return(plot_dpar(
      x,
      dpar = args$which_y,
      facet_by = facet_by,
      color_by = color_by,
      at = at,
      lines = lines,
      cp_dens = cp_dens,
      prior = prior,
      arma = arma,
      ndraws = ndraws
    ))
  }

  get_plot(
    x,
    q_fit = q_fit,
    q_predict = q_predict,
    facet_by = facet_by,
    color_by = color_by,
    at = at,
    lines = lines,
    geom_data = geom_data,
    cp_dens = cp_dens,
    rate = rate,
    prior = prior,
    dpar = NULL,
    arma = arma,
    ndraws = ndraws,
    scale = "response",
    .grouping = grouping,
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
                     ndraws = 1000,
                     scale = "response",
                     at = NULL,
                     nsamples = lifecycle::deprecated(),
                     ...) {
  grouping = if (missing(color_by) && missing(facet_by)) "auto" else "mapped"
  if (!missing(color_by) && is.null(color_by)) grouping = "none"
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "plot_dpar")
  get_plot(
    x,
    q_fit = q_fit,
    q_predict = FALSE,
    facet_by = facet_by,
    color_by = color_by,
    at = at,
    lines = lines,
    geom_data = FALSE,
    cp_dens = cp_dens,
    rate = TRUE,
    prior = prior,
    dpar = dpar,
    arma = arma,
    ndraws = ndraws,
    scale = scale,
    .grouping = grouping,
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
geom_cp_density = function(fit, facet_by, prior, limits_y, use_color = FALSE) {
  dens_scale = 0.2  # Proportion of plot height
  dens_cut = 0.05  # How much to move density down. 5% is ggplot default. Move a bit further.

  # facet_by will expand by group in mcp_draws(). Categorical cols share
  # parameters across facets, so only expand for group-level effects.
  model_tables = get_fit_model_tables(fit)
  cps = model_tables$cps
  group_cols = unique(stats::na.omit(model_tables$group_effects$group_col))
  facet_by = intersect(facet_by, group_cols)
  cp_matches_facet = cps$group_col %in% facet_by
  cp_not_facet = cp_matches_facet == FALSE | is.na(cp_matches_facet)
  if (all(cp_not_facet))
    facet_by = NULL

  # Get group-level and population-level change-point parameter names.
  if (!is.null(facet_by)) {
    varying = stats::na.omit(cps$group_name[cp_matches_facet])
    population = stats::na.omit(cps$name[cp_not_facet])
  } else {
    varying = NULL
    population = cps$name
  }

  # Get draws in long format
  draws = mcp_draws(fit, population = population, varying = varying, absolute = TRUE, prior = prior) %>%
    tidyr::pivot_longer(cols = tidyselect::matches("^cp_[0-9]+$"), names_to = "cp_name", values_to = "value") %>%

    # Compute density per group. Tolerate zero-variance CPs like cp_2 = 80.
    dplyr::group_by(dplyr::across(dplyr::all_of(c(".chain", "cp_name", facet_by)))) %>%
    dplyr::summarise(dens = list(
      if (is.na(stats::sd(.data$value)) || stats::sd(.data$value) == 0) {
        list(x = rep(.data$value[1], 2), y = c(0, 1))
      } else {
        tryCatch(
          stats::density(.data$value, bw = "SJ", n = 2^10),
          error = function(e) stats::density(.data$value, bw = "nrd0", n = 2^10)
        )
      }
    )) %>%
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
#' @aliases geom_quantiles
#' @keywords internal
#' @noRd
#' @inheritParams plot.mcpfit
#' @param color Fixed line color when `use_color = FALSE`.
#' @param ... Arguments passed to geom_line
#' @return A `ggplot2::geom_line` object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
geom_quantiles = function(data_quantiles, xvar, yvar, use_color = TRUE, color = NULL, ...) {
  quantile_mapping = if (use_color) {
    ggplot2::aes(y = .data[[as.character(yvar)]], group = interaction(.data$quantile, .data$.group), color = .data$.color)
  } else {
    ggplot2::aes(y = .data[[as.character(yvar)]], group = interaction(.data$quantile, .data$.group))
  }

  if (use_color) {
    return(ggplot2::geom_line(mapping = quantile_mapping, data = data_quantiles, ...))
  } else {
    ggplot2::geom_line(mapping = quantile_mapping, data = data_quantiles, color = color, ...)
  }
}
