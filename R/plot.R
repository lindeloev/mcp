#' Plot mcpfit
#'
#' @aliases plot plot.mcpfit
#' @param x An mcpfit object
#' @param type String. One of "overlay" (default) or "combo".
#' @param draws Positive integer. Number of posterior draws to use when type = "overlay".
#' @param ... Currently ignored.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A \code{ggplot2} object.
#' @import dplyr ggplot2
#' @importFrom stats sd
#' @importFrom grDevices rgb
#' @importFrom purrr invoke
#' @export
#' @examples
#' \dontrun{
#' plot(fit)  # defaults to the below
#' plot(fit, type="overlay", draws=50) + ggtitle("Great fit!")
#' plot(fit, "combo")
#' }

plot.mcpfit = function(x, type="overlay", draws=25, ...) {
  # Check arguments
  if(class(x) != "mcpfit")
    stop("Can only plot mcpfit objects. x was class: ", class(x))
  if(! type %in% c("overlay", "combo"))
    stop("Type has to be one of 'overlay' or 'combo'. Was: ", type)
  if(draws < 1)
    stop("Draws has to be a positive integer.")


  # Plot function on top
  if(type == "overlay") {
    func_y = x$func_y

    eval_at = seq(min(x$data[, x$pars$x]),
                  max(x$data[, x$pars$x]),
                  length.out = 100)

    # First, let's get all the predictors in shape for func_y
    Q = x$samples %>%
      tidybayes::tidy_draws() %>%
      sample_n(draws) %>%
      select(-starts_with(".")) %>%  # Not arguments for func_y so delete
      tidyr::expand_grid(!!x$pars$x := eval_at) %>%  # correct name of x-var

      # Add fitted draws (vectorized)
      mutate(!!x$pars$y := purrr::invoke(func_y, ., type="fitted")) %>%

      # Add group
      mutate(line = rep(1:draws, each=length(eval_at)))

    # Plot it
    ggplot(x$data, aes_string(x = x$pars$x, y = x$pars$y)) +
      geom_line(aes(group = line), data = Q, color=rgb(0.5, 0.5, 0.5, 0.4)) +
      geom_point()
  }

  else if(type == "combo") {
    pars = paste0("^int_|^cp_|^", x$pars$x, "_|sigma")
    bayesplot::mcmc_combo(x$samples, regex_pars=pars)
  }
}


#' Summarise mcpfit
#'
#' A simple summary of mcpfit objects. It will be updated later.
#'
#' @aliases summary summary.mcpfit
#' @param object \code{mcpfit} object.
#' @param width The width of the Higest Density Interval.
#' @param ... Currently ignored.
#' @export
#' @examples
#' \dontrun{
#' summary(fit)
#' }

summary.mcpfit = function(object, width = 0.95, ...) {
  if(!is.null(object$samples)) {
    object$samples %>%
      tidybayes::tidy_draws() %>%
      tidyr::pivot_longer(-starts_with(".")) %>%
      group_by(name) %>%
      tidybayes::mean_hdci(value, .width = width) %>%
      rename(mean = value) %>%
      select(-.point, -.width, -.interval)
  }
  else {
    message("No samples. Nothing to summarise.")
  }
}

#' Print mcpfit
#'
#' @aliases print print.mcpfit
#' @param x \code{mcpfit} object.
#' @param ... Currently ignored.
#' @export
print.mcpfit = function(x, ...) {
  if(!is.null(x$samples)) {
    print(summary(x))
  }
  else {
    message("No samples. Nothing to summarise.")
  }
}
