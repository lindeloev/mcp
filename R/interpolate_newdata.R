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
#' @aliases get_x_values
#' @keywords internal
#' @noRd
#' @inheritParams plot.mcpfit
#' @param fit An `mcpfit` object.
#' @param by Character vector of grouping columns to evaluate separately.
#' @return A vector of x-values to evaluate at.
get_x_values = function(fit, by = NULL, prior = FALSE) {
  N_BASIS = 100
  N_CP = 50
  X_RESOLUTION_GROUPED = 300

  xdata = fit$data[, fit$pars$x] %>% as.numeric()

  # If there are AR/MA terms, evaluate at the data
  if (length(fit$pars$arma) > 0) {
    x_values = xdata
  } else if (!is.null(by) || is.null(fit$mcmc_post)) {
    # Just give up for grouped and prior evaluations (usually very distributed change points)
    # and return a reasonable resolution
    x_values = seq(min(xdata), max(xdata), length.out = X_RESOLUTION_GROUPED)
  } else if (length(fit$pars$cp) == 0) {
    # No change points. Use default resolution for the whole plot
    x_values = seq(min(xdata), max(xdata), length.out = N_BASIS)
  } else {
    # Make fine resolution around change points in addition to course resolution
    # Get samples for these change points
    samples = mcmclist_samples(fit, prior = prior)
    cp_pars = fit$pars$cp
    call = paste0("tidybayes::spread_draws(samples, ", paste0(cp_pars, collapse = ", "), ")")
    samples = eval(str2lang(call))

    # Compute and return
    x_values = sort(c(
      seq(min(xdata), max(xdata), length.out = N_BASIS),  # Default resolution for the whole plot
      unlist(lapply(cp_pars, function(cp_par) unname(stats::quantile(samples[[cp_par]], probs = seq(0, 1, length.out = N_CP)))))  # Higher res at change points
    ))
  }
  return(x_values)
}

#' List of interpolated values at the values in "at".
#'
#' @aliases interpolate_continuous
#' @keywords internal
#' @noRd
#' @param data fit$data
#' @param pars fit$pars
#' @param x_values par_x values to interpolate continuous predictors at.
#' @param varying_cols Varying-effect columns to exclude from continuous
#'   interpolation.
#' @return `data.frame` with one column for each continuous predictor.
#'   `NULL` if there are no continous predictors.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
interpolate_continuous = function(data, pars, x_values, varying_cols = NULL) {
  checkmate::assert_data_frame(data)
  checkmate::assert_list(pars)
  checkmate::assert_numeric(x_values, any.missing = FALSE)

  # Get numeric predictor columns
  numeric_data = data[, sapply(data, is.numeric), drop = FALSE]
  numeric_data = numeric_data[, colnames(numeric_data) %notin% c(pars$x, pars$y, pars$weights, varying_cols), drop = FALSE]

  if (ncol(numeric_data) == 0)
    return(NULL)

  # Return interpolated
  numeric_data %>%
    lapply(function(col) stats::approx(x = dplyr::pull(data, pars$x), y = col, xout = x_values)$y) %>%
    as.data.frame()
}


#' Returns a data.frame with all combos of predictors
#'
#' This function synthesizes predictors for all combinations of predictor values.
#' It is used internally in `plot.mcpfit()` and may be useful if you want to
#' build your own custom plot.
#'
#' @aliases interpolate_newdata
#' @param fit An `mcpfit` object.
#' @param by Character vector of categorical or varying-effect columns to evaluate separately.
#'   Categorical model predictors are always included.
#' @param x_values Numeric vector of x-values to interpolate at.
#' @details
#' The `par_x` variable will be interpolated with higher resolution around the
#' change points where the values can change abruptly, but lower resolution in
#' between to speed up the computation.
#'
#' Categorical variables and requested varying-effect groups are combined factorially (all level combinations).
#' Continuous variables are interpolated at the x-values and applied to every curve.
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
#' fit = mcp_example("multiple")
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
interpolate_newdata = function(fit, by = NULL, x_values = get_x_values(fit, by)) {

  # Get unique predictors
  varying_cols = unique(stats::na.omit(get_fit_model_tables(fit)$group_effects$group_col))
  by = unique(c(names(get_categorical_levels(fit$data)), intersect(varying_cols, by)))
  # Numeric varying-group IDs are discrete even when they are not requested
  # for evaluation, so never interpolate them as continuous predictors.
  by_grid = lapply(fit$data[, by, drop = FALSE], unique) %>% expand.grid()
  has_groups = nrow(by_grid) > 0 | length(colnames(by_grid) %notin% by) > 0
  has_continuous = interpolate_continuous(fit$data, fit$pars, x_values[1], varying_cols) %>% is.null() %>% `!`

  # Return with levels, if such exist
  if (!has_groups & !has_continuous) {
    newdata = tibble::tibble("{fit$pars$x}" := x_values)
  } else  if (has_groups & !has_continuous) {
    newdata = by_grid %>%
      tidyr::expand_grid("{fit$pars$x}" := x_values)
  } else if (!has_groups & has_continuous) {
    newdata = interpolate_continuous(fit$data, fit$pars, x_values, varying_cols) %>%
      dplyr::mutate("{fit$pars$x}" := x_values)
  } else if (has_groups & has_continuous) {
    # Interpolate continuous predictors within each row of by_grid
    # and up/down-fill if outside the observed region.
    df_list = list()
    for (i in seq_len(nrow(by_grid))) {
      data_i = dplyr::left_join(by_grid[i, , drop = FALSE], fit$data) %>% suppressMessages()
      interpolated_i = interpolate_continuous(data_i, fit$pars, x_values, varying_cols) %>%
        tidyr::fill(dplyr::everything(), .direction = "downup")

      df_list[[i]] = by_grid[i, , drop = FALSE] %>%
        tidyr::expand_grid(interpolated_i) %>%
        dplyr::mutate("{fit$pars$x}" := x_values)
    }

    newdata = dplyr::bind_rows(df_list)
  } else {
    stop_github("Not one of the possible combos of categorical and continuous.")
  }

  # Add response column for AR/MA models
  if (is_arma(fit)) {
    if (nrow(newdata) != nrow(fit$data))
      stop_github("nrow(newdata) != nrow(fit$data) in interpolate_newdata for an AR/MA model.")
    newdata[, fit$pars$y] = fit$data[, fit$pars$y]
  }

  as.data.frame(newdata)
}
