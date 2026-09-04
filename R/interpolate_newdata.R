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
get_x_values = function(fit, by = NULL, prior = FALSE, arma = NULL) {
  N_BASIS = 100
  N_CP = 50
  X_RESOLUTION_GROUPED = 300

  if (is.null(arma))
    arma = has_arma_terms(fit)

  data_columns = mcp_columns(fit)
  xdata = fit$data[, data_columns$par_x] %>% as.numeric()

  # If there are AR/MA terms and arma == TRUE, evaluate at the data
  if (isTRUE(arma)) {
    x_values = xdata
  } else if (!is.null(by) || is.null(.subset2(fit, "mcmc_post"))) {
    # Just give up for grouped and prior evaluations (usually very distributed change points)
    # and return a reasonable resolution
    x_values = seq(min(xdata), max(xdata), length.out = X_RESOLUTION_GROUPED)
  } else if (nrow(get_fit_model_tables(fit)$cps) == 0) {
    # No change points. Use default resolution for the whole plot
    x_values = seq(min(xdata), max(xdata), length.out = N_BASIS)
  } else {
    # Make fine resolution around change points in addition to coarse resolution
    draws = mcmclist_draws(fit, prior = prior)
    cp_pars = get_fit_model_tables(fit)$cps$name
    draws_mat = as.matrix(draws)[, cp_pars, drop = FALSE]

    # Compute and return
    x_values = sort(c(
      seq(min(xdata), max(xdata), length.out = N_BASIS),  # Default resolution for the whole plot
      unlist(lapply(cp_pars, function(cp_par) unname(stats::quantile(draws_mat[, cp_par], probs = seq(0, 1, length.out = N_CP)))))  # Higher res at change points
    ))
  }
  x_values
}

#' Get fixed values for continuous predictors
#'
#' @aliases get_continuous_at
#' @keywords internal
#' @noRd
#' @param data fit$data
#' @param data_columns Fitting-data column metadata.
#' @param at Named list overriding the default means.
#' @param group_cols Grouping columns to exclude from continuous
#'   predictors.
#' @return `data.frame` with one column for each continuous predictor.
#'   `NULL` if there are no continuous predictors.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_continuous_at = function(data, data_columns, at = NULL, group_cols = NULL) {
  checkmate::assert_data_frame(data)
  checkmate::assert_list(data_columns)

  # Get numeric predictor columns
  auxiliary_cols = unname(unlist(data_columns[setdiff(names(data_columns), c("par_x", "response", "series"))]))
  numeric_data = data[, sapply(data, is.numeric), drop = FALSE]
  numeric_data = numeric_data[, colnames(numeric_data) %notin% c(data_columns$par_x, data_columns$response, group_cols, auxiliary_cols), drop = FALSE]

  checkmate::assert_list(at, types = "numeric", any.missing = FALSE, names = "unique", null.ok = TRUE)
  if (length(at) > 0 && is.null(names(at)))
    stop("`at` must be a named list.")
  invalid_at = setdiff(names(at), c(names(numeric_data), auxiliary_cols))
  if (length(invalid_at) > 0)
    stop("`at` must name continuous predictors or response auxiliaries. Invalid: '", paste(invalid_at, collapse = "', '"), "'.")
  if (any(lengths(at) != 1))
    stop("Every value in `at` must be a single number.")

  values = lapply(numeric_data, mean, na.rm = TRUE)
  values[names(at)] = at
  if (length(values) == 0)
    return(NULL)
  as.data.frame(values)
}


#' Returns a data.frame with all combos of predictors
#'
#' \lifecycle{experimental}
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
#'   They default to their observed means. Family response auxiliaries can also
#'   be supplied as explicit scalar design values; e.g., `at = list(N = 20)`.
#' @param arma Logical. If `TRUE`, preserve the observed response history for
#'   conditional AR/MA evaluation. Defaults to `TRUE` when the fit includes
#'   `ar()` or `ma()` terms. Set to `FALSE` to interpolate unconditional trends.
#' @details
#' The `par_x` variable will be interpolated with higher resolution around the
#' change points where the values can change abruptly, but lower resolution in
#' between to speed up the computation.
#'
#' Categorical variables and requested grouping factors are combined factorially (all level combinations).
#' Additional continuous predictors are held at their observed means, or at values supplied through `at`.
#' Family-specific response auxiliaries are not interpolated. Supply binomial
#' trial counts as a scalar in `at` or use `newdata` for a varying design.
#' 
#' Likelihood weights are needed only when evaluating `log_lik()` and
#' therefore need not be supplied here.
#' @return `data.frame` with
#'  * Cols for par_x
#'  * unique levels combos of factorial vars
#'  * fixed values for additional continuous predictors
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \donttest{
#' # Get predictors for a fit
#' newdata = interpolate_newdata(demo_fit)
#'
#' # Fit summary
#' head(fitted(demo_fit, newdata))
#'
#' # Predictions for each draw
#' prediction = predict(demo_fit, newdata, summary = FALSE)
#' head(prediction)
#'
#' # Custom plot
#' library(ggplot2)
#' plotdata = fitted(demo_fit, newdata)
#' ggplot(plotdata, aes(x = time, y = fitted)) +
#'   geom_ribbon(aes(ymin = `Q2.5`, ymax = `Q97.5`), alpha = 0.3) +
#'   geom_line(lwd = 2) +
#'   geom_point(aes(y = response), data = demo_fit$data)
#' }
interpolate_newdata = function(fit, by = NULL, x_values = NULL, at = NULL, arma = NULL) {
  if (is.null(arma))
    arma = has_arma_terms(fit)

  # Conditional AR/MA predictions depend sequentially on the observed response
  # history. Preserve the observed design rather than creating a factorial grid.
  if (isTRUE(arma) && is.null(x_values)) {
    newdata = fit$data
    if (!is.null(at)) {
      for (col in intersect(names(at), names(newdata)))
        newdata[[col]] = at[[col]]
    }
    return(as.data.frame(newdata))
  }

  if (is.null(x_values))
    x_values = get_x_values(fit, by = by, arma = arma)
  # Evaluate the default before `by` is normalised below. Otherwise an absent
  # grouping argument becomes character(0), which selects the denser grouped
  # grid in get_x_values().
  x_values = force(x_values)

  # Get unique predictors
  data_columns = mcp_columns(fit)
  group_cols = unique(stats::na.omit(get_fit_model_tables(fit)$group_effects$group_col))
  series_col = data_columns$series
  categorical_cols = setdiff(
    names(get_categorical_levels(fit$data)),
    c(group_cols, series_col)
  )
  by = unique(c(categorical_cols, intersect(setdiff(group_cols, series_col), by)))
  # Numeric group IDs are discrete even when they are not requested
  # for evaluation, so never interpolate them as continuous predictors.
  by_grid = if (length(by) == 0) {
    data.frame(.row = 1)[, FALSE, drop = FALSE]
  } else {
    lapply(fit$data[, by, drop = FALSE], unique) %>% expand.grid()
  }
  continuous_at = get_continuous_at(
    fit$data, data_columns, at, unique(c(group_cols, series_col))
  )
  newdata = by_grid %>% tidyr::expand_grid("{data_columns$par_x}" := x_values)
  if (!is.null(continuous_at))
    newdata = tidyr::expand_grid(newdata, continuous_at)

  auxiliary_cols = unname(unlist(data_columns[setdiff(names(data_columns), c("par_x", "response", "series"))]))
  auxiliary_cols = setdiff(
    unname(unlist(data_columns[setdiff(names(data_columns), c("par_x", "response", "series"))])),
    data_columns$par_x  # Don't remove par_x when it is also a rhs parameter, e.g.: y | trials(N) ~ N
  )
  if (any(auxiliary_cols %in% names(newdata))) {
    model_tables = get_fit_model_tables(fit)
    response_data = get_family_response_data(fit$family, model_tables$segments, newdata)
    response_columns = c(
      list(y = data_columns$response),
      as.list(get_family_aux_columns(fit$family, model_tables$segments))
    )
    fit$family$response$validate(rep(NA_real_, nrow(newdata)), response_data, response_columns)
  }

  # Add response column for AR/MA models
  if (isTRUE(arma)) {
    if (nrow(newdata) != nrow(fit$data))
      stop(
        "Conditional AR/MA evaluation (`arma = TRUE`) requires an observed response history ",
        "matching the fitted data. Set `arma = FALSE` to evaluate unconditional trends on new or interpolated data."
      )
    if (!is.null(series_col))
      newdata[, series_col] = fit$data[, series_col]
    newdata[, data_columns$response] = fit$data[, data_columns$response]
  }

  as.data.frame(newdata)
}
