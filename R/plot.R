#' Plot mcpfit
#'
#' @param fit An mcpfit object
#' @param type String. One of "overlay" (default) or "combo".
#' @param draws Integer. Number of posterior draws to use when type = "overlay".
#' @keywords mcmcfit, plot
#' @import dplyr
#' @import tidybayes
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
      tidy_draws() %>%
      sample_n(draws) %>%
      select(-starts_with(".")) %>%
      expand_grid(!!fit$pars$x := eval_at) %>%  # correct name of x-var
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
