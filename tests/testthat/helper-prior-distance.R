# Maximum allowed distance (in prior scales) between a simulated value and its
# default prior's location. Default coefficient priors are normal, so values beyond
# 3 prior SDs are noticeably shrunk.
PRIOR_DISTANCE_MAX = 3


#' Distance of simulated values from their prior, in prior scales
#'
#' @keywords internal
#' @param fit An `mcpfit`, typically unsampled.
#' @param values Named list of simulated parameter values.
#' @return Data frame with `parameter`, `value`, `location`, `scale`, and `distance`.
#'   Only parameters with a numeric `dt()` or `dnorm()` prior are included.
prior_distance = function(fit, values) {
  rows = lapply(intersect(names(values), names(fit$prior)), function(name) {
    prior = fit$prior[[name]]
    if (!is.character(prior) || length(values[[name]]) != 1)
      return(NULL)
    call = parse_prior_call(split_prior_truncation(prior)$distribution)
    if (is.null(call) || call$name %notin% c("dt", "dnorm"))
      return(NULL)
    location = suppressWarnings(as.numeric(call$args[1]))
    scale = suppressWarnings(as.numeric(call$args[2]))
    if (is.na(location) || is.na(scale))
      return(NULL)
    value = values[[name]]
    data.frame(parameter = name, value = value, location = location, scale = scale, distance = abs(value - location) / scale)
  })
  do.call(rbind, rows)
}


#' Expect simulated values to lie within `PRIOR_DISTANCE_MAX` prior scales
#'
#' @keywords internal
#' @inheritParams prior_distance
#' @param label Description used in the failure message.
#' @param ignore Parameter names exempt from the check.
expect_prior_distance = function(fit, values, label, ignore = character()) {
  distances = prior_distance(fit, values)
  if (is.null(distances))
    return(invisible(NULL))
  far = distances[distances$distance > PRIOR_DISTANCE_MAX & distances$parameter %notin% ignore, ]
  testthat::expect(
    nrow(far) == 0,
    paste0(label, ": simulated values far from their default priors: ",
           paste0(far$parameter, " (", round(far$distance, 1), " scales)", collapse = ", "))
  )
  invisible(distances)
}
