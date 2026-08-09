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
  } else if (!is.null(by) || is.null(.subset2(fit, "mcmc_post"))) {
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

#' Get fixed values for continuous predictors
#'
#' @aliases get_continuous_at
#' @keywords internal
#' @noRd
#' @param data fit$data
#' @param pars fit$pars
#' @param at Named list overriding the default means.
#' @param group_cols Grouping columns to exclude from continuous
#'   predictors.
#' @return `data.frame` with one column for each continuous predictor.
#'   `NULL` if there are no continuous predictors.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_continuous_at = function(data, pars, at = NULL, group_cols = NULL) {
  checkmate::assert_data_frame(data)
  checkmate::assert_list(pars)

  # Get numeric predictor columns
  numeric_data = data[, sapply(data, is.numeric), drop = FALSE]
  numeric_data = numeric_data[, colnames(numeric_data) %notin% c(pars$x, pars$y, pars$weights, group_cols), drop = FALSE]

  checkmate::assert_list(at, types = "numeric", any.missing = FALSE, names = "unique", null.ok = TRUE)
  if (length(at) > 0 && is.null(names(at)))
    stop("`at` must be a named list.")
  invalid_at = setdiff(names(at), names(numeric_data))
  if (length(invalid_at) > 0)
    stop("`at` must name continuous predictors other than `par_x`. Invalid: '", paste(invalid_at, collapse = "', '"), "'.")
  if (any(lengths(at) != 1))
    stop("Every value in `at` must be a single number.")

  if (ncol(numeric_data) == 0)
    return(NULL)

  values = lapply(numeric_data, mean, na.rm = TRUE)
  values[names(at)] = at
  as.data.frame(values)
}


#' Returns a data.frame with all combos of predictors
#'
#' This function synthesizes predictors for all combinations of predictor values.
#' It is used internally in `plot.mcpfit()` and may be useful if you want to
#' build your own custom plot.
#'
#' @aliases interpolate_newdata
#' @param fit An `mcpfit` object.
#' @param by Character vector of categorical or group-level columns to evaluate separately.
#'   Categorical model predictors are always included.
#' @param x_values Numeric vector of x-values to evaluate at.
#' @param at Named list setting additional continuous predictors to fixed values.
#'   They default to their observed means. For example, `at = list(age = 40)`.
#' @details
#' The `par_x` variable will be interpolated with higher resolution around the
#' change points where the values can change abruptly, but lower resolution in
#' between to speed up the computation.
#'
#' Categorical variables and requested grouping factors are combined factorially (all level combinations).
#' Additional continuous predictors are held at their observed means, or at values supplied through `at`.
#' @return `tibble` with
#'  * Cols for par_x
#'  * unique levels combos of factorial vars
#'  * fixed values for additional continuous predictors
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
interpolate_newdata = function(fit, by = NULL, x_values = get_x_values(fit, by), at = NULL) {
  # Evaluate the default before `by` is normalised below. Otherwise an absent
  # grouping argument becomes character(0), which selects the denser grouped
  # grid in get_x_values().
  x_values = force(x_values)

  # Get unique predictors
  group_cols = unique(stats::na.omit(get_fit_model_tables(fit)$group_effects$group_col))
  categorical_cols = setdiff(names(get_categorical_levels(fit$data)), group_cols)
  by = unique(c(categorical_cols, intersect(group_cols, by)))
  # Numeric group IDs are discrete even when they are not requested
  # for evaluation, so never interpolate them as continuous predictors.
  by_grid = if (length(by) == 0) {
    data.frame(.row = 1)[, FALSE, drop = FALSE]
  } else {
    lapply(fit$data[, by, drop = FALSE], unique) %>% expand.grid()
  }
  continuous_at = get_continuous_at(fit$data, fit$pars, at, group_cols)
  newdata = by_grid %>% tidyr::expand_grid("{fit$pars$x}" := x_values)
  if (!is.null(continuous_at))
    newdata = tidyr::expand_grid(newdata, continuous_at)

  # Add response column for AR/MA models
  if (is_arma(fit)) {
    if (nrow(newdata) != nrow(fit$data))
      stop_github("nrow(newdata) != nrow(fit$data) in interpolate_newdata for an AR/MA model.")
    newdata[, fit$pars$y] = fit$data[, fit$pars$y]
  }

  as.data.frame(newdata)
}
