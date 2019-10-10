#' Plot mcpfit
#'
#' @param fit An mcpfit object
#' @param type String. One of "overlay" (default) or "combo".
#' @param draws Integer. Number of posterior draws to use when type = "overlay".
#' @import tidyr dplyr ggplot2
#' @export
#' @examples
#' plot(fit)  # defaults to the below
#' plot(fit, type="overlay", draws=50)
#' plot(fit, "combo")

plot.mcpfit = function(fit, type="overlay", draws=25) {
  # Plot function on top
  if(type == "overlay") {
    eval_at = seq(min(fit$data[, fit$pars$x]),
                  max(fit$data[, fit$pars$x]),
                  length.out = 100)

    # First, let's get all the predictors in shape for func_y
    Q = fit$samples %>%
      tidybayes::tidy_draws() %>%
      sample_n(draws) %>%
      select(-starts_with(".")) %>%
      tidyr::expand_grid(!!fit$pars$x := eval_at) %>%  # correct name of x-var
      mutate(type = "fitted") %>%

      # Now we make a nested table for each row and apply them as args to func_y
      rowwise() %>%
      do(args = as.data.frame(.)) %>%
      mutate(
        !!fit$pars$y := purrr::pmap(args, fit$func_y),  # correct name of y-var
      ) %>%

      # Now finish up: add group (for ggplot) and unnest it all
      ungroup() %>%
      mutate(group = rep(1:draws, each=length(eval_at))) %>%
      unnest(c(args, !!fit$pars$y))

    ggplot(fit$data, aes_string(x = fit$pars$x, y = fit$pars$y)) +
      geom_line(aes(group = group), data = Q, color=rgb(0.5, 0.5, 0.5, 0.4)) +
      geom_point()
  }

  else if(type == "combo") {
    pars = paste0("^int_|^cp_|^", fit$pars$x, "_|sigma")
    bayesplot::mcmc_combo(fit$samples, regex_pars=pars)
  }
}


#' Plot mcpfit
#'
#' A simple summary of mcpfit objects. It will be updated later.
#'
#' @param fit \code{mcpfit} object.
#' @param width The width of the Higest Density Interval.
#' @export
#' @example summary(fit)

summary.mcpfit = function(fit, width = 0.95) {
  fit$samples %>%
    tidy_draws() %>%
    pivot_longer(-starts_with(".")) %>%
    group_by(name) %>%
    mean_hdci(value, .width = width) %>%
    rename(mean = value) %>%
    select(-.point, -.width, -.interval)
}

#' Print mcpfit
#'
#' @param fit \code{mcpfit} object.
#' @export
print.mcpfit = function(fit) {
  print(summary(fit))
}
