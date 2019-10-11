#' Plot mcpfit
#'
#' @aliases plot
#' @param x An mcpfit object
#' @param type String. One of "overlay" (default) or "combo".
#' @param draws Integer. Number of posterior draws to use when type = "overlay".
#' @param ... Currently ignored.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A \code{ggplot2} object.
#' @import dplyr ggplot2
#' @importFrom stats sd
#' @importFrom grDevices rgb
#' @export
#' @examples
#' plot(fit)  # defaults to the below
#' plot(fit, type="overlay", draws=50) + ggtitle("Great fit!")
#' plot(fit, "combo")

plot.mcpfit = function(x, type="overlay", draws=25, ...) {
  # Plot function on top
  if(type == "overlay") {
    eval_at = seq(min(x$data[, x$pars$x]),
                  max(x$data[, x$pars$x]),
                  length.out = 100)

    # First, let's get all the predictors in shape for func_y
    Q = x$samples %>%
      tidybayes::tidy_draws() %>%
      sample_n(draws) %>%
      select(-starts_with(".")) %>%
      tidyr::expand_grid(!!x$pars$x := eval_at) %>%  # correct name of x-var
      mutate(type = "fitted") %>%

      # Now we make a nested table for each row and apply them as args to func_y
      rowwise() %>%
      do(args = as.data.frame(.)) %>%
      mutate(
        !!x$pars$y := purrr::pmap(args, x$func_y),  # correct name of y-var
      ) %>%

      # Now finish up: add group (for ggplot) and unnest it all
      ungroup() %>%
      mutate(group = rep(1:draws, each=length(eval_at))) %>%
      tidyr::unnest(c(args, !!x$pars$y))

    ggplot(x$data, aes_string(x = x$pars$x, y = x$pars$y)) +
      geom_line(aes(group = group), data = Q, color=rgb(0.5, 0.5, 0.5, 0.4)) +
      geom_point()
  }

  else if(type == "combo") {
    pars = paste0("^int_|^cp_|^", x$pars$x, "_|sigma")
    bayesplot::mcmc_combo(x$samples, regex_pars=pars)
  }
}


#' Plot mcpfit
#'
#' A simple summary of mcpfit objects. It will be updated later.
#'
#' @aliases summary
#' @param object \code{mcpfit} object.
#' @param width The width of the Higest Density Interval.
#' @param ... Currently ignored.
#' @export
#' @example summary(fit)

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
#' @aliases print
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
