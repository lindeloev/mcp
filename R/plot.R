# ABOUT: These are functions directly related to plotting
# ----------------

#' Plot full fits
#'
#' Plot prior or posterior model draws on top of data. Use `plot_pars` to
#' plot individual parameter estimates.
#'
#' @aliases plot plot.mcpfit
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
#'   or FALSE (don not plot).
#' @param cp_dens TRUE/FALSE. Plot posterior densities of the change point(s)?
#'   Currently does not respect `facet_by`. This will be added in the future.
#' @param ... Currently ignored.
#' @details
#'   `plot()` uses `fit$simulate()` on posterior samples. These represent the
#'   (joint) posterior distribution.fit
#' @return A \pkg{ggplot2} object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
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
  # IN PROGRESS HERE
  # all_group_args = c(facet_by, color_by)
  # rhs_cat_vars = get_categorical_levels(fit$data)
  # if (length(all_group_args) > 0) {
  # }

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


#' Density geom for `plot.mcpfit()`
#'
#' Note that `geom_density(fill = ...)` always fill area to y = 0 but we want to fill
#' to the bottom of the plot. So we use `geom_polygon(fill = ...)` instead.
#'
#' @aliases geom_cp_density
#' @keywords internal
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


#' Plot individual parameters
#'
#' Plot many types of plots of parameter estimates. See examples for typical use
#' cases.
#'
#' @aliases plot_pars
#' @param fit An \code{\link{mcpfit}} object.
#' @param pars Character vector. One of:
#'   * Vector of parameter names.
#'   * `"population"` plots all population parameters.
#'   * `"varying"` plots all varying effects. To plot a particular varying
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
#' @export
#' @examples
#' # Typical usage. demo_fit is an mcpfit object.
#' plot_pars(demo_fit)
#'
#' \dontrun{
#' # More options
#' plot_pars(demo_fit, regex_pars = "^cp_")  # Plot only change points
#' plot_pars(demo_fit, pars = c("Intercept_3", "time_3"))  # Plot these parameters
#' plot_pars(demo_fit, type = c("trace", "violin"))  # Combine plots
#' # Some plots only take pairs. hex is good to assess identifiability
#' plot_pars(demo_fit, type = "hex", pars = c("cp_1", "time_2"))
#'
#' # Visualize the priors:
#' plot_pars(demo_fit, prior = TRUE)
#'
#' # Useful for varying effects:
#' # plot_pars(my_fit, pars = "varying", ncol = 3)  # plot all varying effects
#' # plot_pars(my_fit, regex_pars = "my_varying", ncol = 3)  # plot all levels of a particular varying
#'
#' # Customize multi-column ggplots using "*" instead of "+" (patchwork)
#' library(ggplot2)
#' plot_pars(demo_fit, type = c("trace", "dens_overlay")) * theme_bw(10)
#' }
plot_pars = function(fit,
                     pars = "population",
                     regex_pars = character(0),
                     type = "combo",
                     ncol = 1,
                     prior = FALSE
                     ) {

  # Check arguments
  assert_types(fit, "mcpfit")

  if (!coda::is.mcmc.list(fit$mcmc_post) && !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

  if (!is.character(pars) || !is.character(regex_pars))
    stop("`pars` and `regex_pars` has to be string/character.")

  if (any(c("population", "varying") %in% pars) && length(pars ) > 1)
    stop("`pars` cannot be a vector that contains multiple elements AND 'population' or 'varying'.")

  if (any(c("hex", "scatter") %in% type) && (length(pars) != 2 || length(regex_pars) > 0))
    stop("`type` = 'hex' or 'scatter' takes exactly two parameters which must be provided via the `pars` argument")

  if ("combo" %in% type && length(type) > 1)
    stop("'combo' type cannot be combined with other types. Replace 'combo' with the types you want combo\'ed")

  assert_integer(ncol, lower = 1, len = 1)
  assert_logical(prior)
  bayesplot::available_mcmc()  # Quick fix to make R CMD Check happy that bayesplot is imported

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

  # Call the relevant bayesplot plot function for each type
  takes_facet = c("areas", "dens", "dens_overlay", "trace", "hist", "intervals", "trace", "trace_highlight", "violin")
  all_plots = list()
  for (this_type in type) {
    if (this_type %in% takes_facet) {
      facet_args = list(ncol = ncol)
    } else {
      facet_args = list()
    }

    func_name = paste0("mcmc_", this_type)
    func_obj = utils::getFromNamespace(func_name, "bayesplot")
    all_plots[[this_type]] = func_obj(samples, pars = pars, regex_pars = regex_pars, facet_args = facet_args)
  }

  # Then patch all_plots together and return
  patchwork::wrap_plots(all_plots) &
    ggplot2::theme(
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
#' @param fit An `mcpfit` object.
#' @return A vector of x-values to evaluate at.
get_eval_at = function(fit, facet_by = NULL, prior = FALSE) {
  N_BASIS = 100
  N_CP = 50
  X_RESOLUTION_FACET = 300  # Only varying

  xdata = fit$data[, fit$pars$x] %>% as.numeric()

  # If there are ARMA terms, evaluate at the data
  if (length(fit$pars$arma) > 0) {
    return(xdata)

    # Just give up for faceting and prior-plots (usually very distributed change points)
    # and return a reasonable resolution
  } else if (!is.null(facet_by) || is.null(fit$mcmc_post)) {
    eval_at = seq(min(xdata), max(xdata), length.out = X_RESOLUTION_FACET)
    return(eval_at)

    # Make regions of fine resolution within course resolution
  } else {
    # Get samples for these change points
    samples = mcmclist_samples(fit, prior = prior)
    cp_pars = get_cp_pars(fit$pars)
    call = paste0("tidybayes::spread_draws(samples, ", paste0(cp_pars, collapse = ", "), ")")
    samples = eval(parse(text = call))

    # Compute and return
    eval_at = sort(c(
      seq(min(xdata), max(xdata), length.out = N_BASIS),  # Default res for whole plot
      unlist(lapply(cp_pars, function(cp_par) unname(stats::quantile(samples[[cp_par]], probs = seq(0, 1, length.out = N_CP)))))  # Higher res at change points
    ))
    return(eval_at)
  }
}

# Add column ".group" which is the interaction of all categorical colnames
add_group = function(df) {
  categorical_colnames = df %>% get_categorical_levels() %>% names()
  if (length(categorical_colnames) == 0) {
    df$.group = as.factor(1)
    df
  } else {
    df %>%
      dplyr::mutate(.group = dplyr::select(., !!!categorical_colnames) %>% interaction())
  }
}

#' List of interpolated values at the values in "at".
#'
#' @aliases interpolate_continuous
#' @keywords internal
#' @param data fit$data
#' @param pars fit$pars
#' @param eval_at par_x values to interpolate continuous predictors at.
#' @return `data.frame` with one column for each continuous predictor.
#'   `NULL` if there are no continous predictors.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
interpolate_continuous = function(data, pars, eval_at) {
  assert_types(data, "data.frame", "tibble")
  assert_types(pars, "list")
  assert_numeric(eval_at)

  # Get numeric RHS data columns
  numeric_data = data[, sapply(data, is.numeric), drop = FALSE]
  numeric_data = numeric_data[, colnames(numeric_data) %notin% c(pars$x, pars$y, pars$weights, pars$varying), drop = FALSE]

  if (ncol(numeric_data) == 0)
    return(NULL)

  # Return interpolated
  numeric_data %>%
    lapply(function(col) stats::approx(x = dplyr::pull(data, pars$x), y = col, xout = eval_at)$y) %>%
    as.data.frame()
}


#' Returns a data.frame with all combos of predictors
#'
#' This function is used to synthesize predictors for all combos of RHS predictors.
#' It is used internally in `plot.mcpfit()` and may be useful if you want to
#' build your own custom plot.
#'
#' @aliases interpolate_newdata
#' @inheritParams plot.mcpfit
#' @param fit An `mcpfit` object.
#' @details
#' The `par_x` variable will be interpolated with higher resolution around the
#' change points where the values can change abruptly, but lower resolution in
#' between to speed up the computation.
#'
#' Categorical variables are combined factorially (all level combos) and all
#' continuous variables are interpolated at the x-values (see above paragraph) and
#' applied to all factor-combos.
#' @return `tibble` with
#'  * Cols for par_x
#'  * unique levels combos of factorial vars
#'  * interpolated continuous vars (interpolated within each factorial cell) (fills down/up if outside observed region)
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \dontrun{
#' # Get predictors for a fit
#' fit = mcp_example("multiple", sample = TRUE)$fit
#' newdata = interpolate_newdata(fit)
#'
#' # Fit summary
#' fitted(fit, newdata)
#'
#' # Predictions for each sample
#' prediction = predict(fit, newdata, summary = FALSE)
#' prediction[, c(".chain", ".iteration", ".draw", "x", "group", "z", "predict")]
#'
#' # Custom plot
#' library(ggplot2)
#' newdata = interpolate_newdata(fit)
#' plotdata = fitted(fit, newdata)
#' ggplot(plotdata, aes(x = x, y = fitted, color = group)) +
#'   geom_ribbon(aes(ymin = `Q2.5`, ymax = `Q97.5`, fill = group), alpha = 0.3) +
#'   geom_line(lwd = 2) +
#'   geom_point(aes(y = y), data = fit$data)
#' }
interpolate_newdata = function(fit, facet_by = NULL) {

  # Get unique predictors
  xvar = rlang::sym(fit$pars$x)
  eval_at = get_eval_at(fit, facet_by)
  categorical_interactions = get_categorical_levels(fit$data) %>% expand.grid()
  has_categorical = nrow(categorical_interactions) > 0 | length(colnames(categorical_interactions) %notin% facet_by) > 0
  has_continuous = interpolate_continuous(fit$data, fit$pars, eval_at[1]) %>% is.null() %>% `!`

  # Return with levels, if such exist
  if (!has_categorical & !has_continuous) {
    newdata = tibble::tibble(!!xvar := eval_at)
  } else  if (has_categorical & !has_continuous) {
    newdata = categorical_interactions %>%
      tidyr::expand_grid(!!xvar := eval_at)
  } else if (!has_categorical & has_continuous) {
    newdata = interpolate_continuous(fit$data, fit$pars, eval_at) %>%
      dplyr::mutate(!!xvar := eval_at)
  } else if (has_categorical & has_continuous) {
    # Interpolate continuous predictors within each factorial cell (row in categorical_interactions)
    # and up/down-fill if outside the observed region.
    df_list = list()
    for (i in seq_len(nrow(categorical_interactions))) {
      data_i = dplyr::left_join(categorical_interactions[i, , drop = FALSE], fit$data) %>% suppressMessages()
      interpolated_i = interpolate_continuous(data_i, fit$pars, eval_at) %>%
        tidyr::fill(dplyr::everything(), .direction = "downup")

      df_list[[i]] = categorical_interactions[i, , drop = FALSE] %>%
        tidyr::expand_grid(interpolated_i) %>%
        dplyr::mutate(!!xvar := eval_at)
    }

    newdata = dplyr::bind_rows(df_list)
  } else {
    stop_github("Not one of the possible combos of categorical and continuous.")
  }

  # Add response column for AR models
  if (is_arma(fit)) {
    if (nrow(newdata) != nrow(fit$data))
      stop_github("nrow(newdata) != nrow(fit$data) in interpolate_newdata for an AR model.")
    newdata[, fit$pars$y] = fit$data[, fit$pars$y]
  }

  as.data.frame(newdata)
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
