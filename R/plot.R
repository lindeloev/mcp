#' Plot mcpfit
#'
#' @param fit An mcpfit object
#' @keywords mcmcfit, plot
#' @import bayesplot
#' @import tidybayes
#' @export
#' @examples
#' plot(fit)  # defaults to plot(fit, type="overlay", draws=50)
#' plot(fit, "combo")


plot.mcpfit = function(fit, type="overlay", draws=50) {
  # Plot function on top
  if(type == "overlay") {
    Q = fit$samples %>%
      spread_draws(year[i], y_[i])

    Q %>%
      sample_draws(draws) %>%
      ggplot(aes_string(x = fit$pars$x, y = "y_", group=".draw")) +
      geom_line(color = rgb(0, 0, 0, alpha=0.2)) +
      geom_point(aes_string(y = fit$pars$y, group=1), data = fit$data)
  }

  else if(type == "combo") {
    pars = paste0("^int_|^cp_|^", fit$pars$x, "_|sigma")
    bayesplot::mcmc_combo(fit$samples, regex_pars=pars)
  }
}

