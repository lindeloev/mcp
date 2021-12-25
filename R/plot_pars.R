#' Plot individual parameters
#'
#' Plot many types of plots of parameter estimates. See examples for typical use
#' cases.
#'
#' @aliases plot_pars
#' @param fit An \code{\link{mcpfit}} object.
#' @param pars Character vector. One of:
#'   * Vector of parameter names.
#'   * `"population"` plots all population parameters.
#'   * `"varying"` plots all varying effects. To plot a particular varying
#'       effect, use `regex_pars = "^name"`.
#' @param regex_pars Vector of regular expressions. This will typically just be
#'   the beginning of the parameter name(s), i.e., "^cp_" plots all change
#'   points, "^my_varying" plots all levels of a particular varying effect, and
#'   "^cp_|^my_varying" plots both.
#' @param type String or vector of strings. Calls `bayesplot::mcmc_>>type<<()`.
#'   Common calls are "combo", "trace", and "dens_overlay". Current options include
#'   'acf', 'acf_bar', 'areas', 'areas_ridges', 'combo', 'dens', 'dens_chains',
#'   'dens_overlay', 'hist', 'intervals', 'rank_hist', 'rank_overlay', 'trace',
#'   'trace_highlight', and 'violin".
#' @param ncol Number of columns in plot. This is useful when you have many
#'   parameters and only one plot `type`.
#' @param prior TRUE/FALSE. Plot using prior samples? Useful for `mcp(..., sample = "both")`
#'
#' @details
#'   For other `type`, it calls `bayesplot::mcmc_type()`. Use these
#'   directly on `fit$mcmc_post` or `fit$mcmc_prior` if you want finer
#'   control of plotting, e.g., `bayesplot::mcmc_dens(fit$mcmc_post)`. There
#'   are also a number of useful plots in the \pkg{coda} package, i.e.,
#'   `coda::gelman.plot(fit$mcmc_post)` and `coda::crosscorr.plot(fit$mcmc_post)`
#'
#'   In any case, if you see a few erratic lines or parameter estimates, this is
#'   a sign that you may want to increase argument 'adapt' and 'iter' in \code{\link{mcp}}.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A \pkg{ggplot2} object.
#' @export
#' @examples
#' # Typical usage. demo_fit is an mcpfit object.
#' plot_pars(demo_fit)
#'
#' \dontrun{
#' # More options
#' plot_pars(demo_fit, regex_pars = "^cp_")  # Plot only change points
#' plot_pars(demo_fit, pars = c("Intercept_3", "time_3"))  # Plot these parameters
#' plot_pars(demo_fit, type = c("trace", "violin"))  # Combine plots
#' # Some plots only take pairs. hex is good to assess identifiability
#' plot_pars(demo_fit, type = "hex", pars = c("cp_1", "time_2"))
#'
#' # Visualize the priors:
#' plot_pars(demo_fit, prior = TRUE)
#'
#' # Useful for varying effects:
#' # plot_pars(my_fit, pars = "varying", ncol = 3)  # plot all varying effects
#' # plot_pars(my_fit, regex_pars = "my_varying", ncol = 3)  # plot all levels of a particular varying
#'
#' # Customize multi-column ggplots using "*" instead of "+" (patchwork)
#' library(ggplot2)
#' plot_pars(demo_fit, type = c("trace", "dens_overlay")) * theme_bw(10)
#' }
plot_pars = function(fit,
                     pars = "population",
                     regex_pars = character(0),
                     type = "combo",
                     ncol = 1,
                     prior = FALSE
) {

  # Check arguments
  assert_types(fit, "mcpfit")

  if (!coda::is.mcmc.list(fit$mcmc_post) && !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")

  if (!is.character(pars) || !is.character(regex_pars))
    stop("`pars` and `regex_pars` has to be string/character.")

  if (any(c("population", "varying") %in% pars) && length(pars ) > 1)
    stop("`pars` cannot be a vector that contains multiple elements AND 'population' or 'varying'.")

  if (any(c("hex", "scatter") %in% type) && (length(pars) != 2 || length(regex_pars) > 0))
    stop("`type` = 'hex' or 'scatter' takes exactly two parameters which must be provided via the `pars` argument")

  if ("combo" %in% type && length(type) > 1)
    stop("'combo' type cannot be combined with other types. Replace 'combo' with the types you want combo\'ed")

  assert_integer(ncol, lower = 1, len = 1)
  assert_logical(prior)
  bayesplot::available_mcmc()  # Quick fix to make R CMD Check happy that bayesplot is imported

  # Get posterior/prior samples
  samples = mcmclist_samples(fit, prior = prior)

  # Handle special codes
  if ("population" %in% pars) {
    if (length(regex_pars) == 0) {
      pars = fit$pars$population
    } else {
      # This probably means that the user left pars as default.
      pars = character(0)
    }
  } else if ("varying" %in% pars) {
    # Regex search for varying effects
    regex_pars = paste0("^", fit$pars$varying)
    pars = character(0)
  }

  # Handles combo. Returns a customizable ggplot which "combo" does not.
  if ("combo" %in% type)
    type = c("dens_overlay", "trace")

  # Call the relevant bayesplot plot function for each type
  takes_facet = c("areas", "dens", "dens_overlay", "trace", "hist", "intervals", "trace", "trace_highlight", "violin")
  all_plots = list()
  for (this_type in type) {
    if (this_type %in% takes_facet) {
      facet_args = list(ncol = ncol)
    } else {
      facet_args = list()
    }

    func_name = paste0("mcmc_", this_type)
    func_obj = utils::getFromNamespace(func_name, "bayesplot")
    all_plots[[this_type]] = func_obj(samples, pars = pars, regex_pars = regex_pars, facet_args = facet_args)
  }

  # Then patch all_plots together and return
  patchwork::wrap_plots(all_plots) &
    ggplot2::theme(
      strip.placement = NULL,  # fixes bug: https://github.com/thomasp85/patchwork/issues/132
      legend.position = "none"  # no legend on chains. Takes up too much space
    )
}
