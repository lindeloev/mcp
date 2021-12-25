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
#' @return A vector of x-values to evaluate at.
get_x_values = function(fit, facet_by = NULL, prior = FALSE) {
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
    x_values = seq(min(xdata), max(xdata), length.out = X_RESOLUTION_FACET)
    return(x_values)

    # Make regions of fine resolution within course resolution
  } else {
    # Get samples for these change points
    samples = mcmclist_samples(fit, prior = prior)
    cp_pars = get_cp_pars(fit$pars)
    call = paste0("tidybayes::spread_draws(samples, ", paste0(cp_pars, collapse = ", "), ")")
    samples = eval(str2lang(call))

    # Compute and return
    x_values = sort(c(
      seq(min(xdata), max(xdata), length.out = N_BASIS),  # Default res for whole plot
      unlist(lapply(cp_pars, function(cp_par) unname(stats::quantile(samples[[cp_par]], probs = seq(0, 1, length.out = N_CP)))))  # Higher res at change points
    ))
    return(x_values)
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
#' @noRd
#' @param data fit$data
#' @param pars fit$pars
#' @param x_values par_x values to interpolate continuous predictors at.
#' @return `data.frame` with one column for each continuous predictor.
#'   `NULL` if there are no continous predictors.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
interpolate_continuous = function(data, pars, x_values) {
  assert_types(data, "data.frame", "tibble")
  assert_types(pars, "list")
  assert_numeric(x_values)

  # Get numeric RHS data columns
  numeric_data = data[, sapply(data, is.numeric), drop = FALSE]
  numeric_data = numeric_data[, colnames(numeric_data) %notin% c(pars$x, pars$y, pars$weights, pars$varying), drop = FALSE]

  if (ncol(numeric_data) == 0)
    return(NULL)

  # Return interpolated
  numeric_data %>%
    lapply(function(col) stats::approx(x = dplyr::pull(data, pars$x), y = col, xout = x_values)$y) %>%
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
#' @param x_values Numeric vector of x-values to interpolate at.
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
interpolate_newdata = function(fit, facet_by = NULL, x_values = get_x_values(fit, facet_by)) {

  # Get unique predictors
  xvar = rlang::sym(fit$pars$x)
  categorical_interactions = get_categorical_levels(fit$data) %>% expand.grid()
  has_categorical = nrow(categorical_interactions) > 0 | length(colnames(categorical_interactions) %notin% facet_by) > 0
  has_continuous = interpolate_continuous(fit$data, fit$pars, x_values[1]) %>% is.null() %>% `!`

  # Return with levels, if such exist
  if (!has_categorical & !has_continuous) {
    newdata = tibble::tibble(!!xvar := x_values)
  } else  if (has_categorical & !has_continuous) {
    newdata = categorical_interactions %>%
      tidyr::expand_grid(!!xvar := x_values)
  } else if (!has_categorical & has_continuous) {
    newdata = interpolate_continuous(fit$data, fit$pars, x_values) %>%
      dplyr::mutate(!!xvar := x_values)
  } else if (has_categorical & has_continuous) {
    # Interpolate continuous predictors within each factorial cell (row in categorical_interactions)
    # and up/down-fill if outside the observed region.
    df_list = list()
    for (i in seq_len(nrow(categorical_interactions))) {
      data_i = dplyr::left_join(categorical_interactions[i, , drop = FALSE], fit$data) %>% suppressMessages()
      interpolated_i = interpolate_continuous(data_i, fit$pars, x_values) %>%
        tidyr::fill(dplyr::everything(), .direction = "downup")

      df_list[[i]] = categorical_interactions[i, , drop = FALSE] %>%
        tidyr::expand_grid(interpolated_i) %>%
        dplyr::mutate(!!xvar := x_values)
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
