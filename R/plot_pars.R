#' Plot individual parameters
#'
#' Plot many types of plots of parameter estimates. See examples for typical use
#' cases.
#'
#' @aliases plot_pars
#' @param fit An \code{\link{mcpfit}} object.
#' @param pars Character vector. One of:
#'   * Vector of parameter names.
#'   * `"population"` plots all population-level parameters.
#'   * `"group"` plots all group-level deviations. To plot a particular group-level
#'       effect, use `regex_pars = "^name"`.
#' @param regex_pars Vector of regular expressions. This will typically just be
#'   the beginning of the parameter name(s), i.e., "^cp_" plots all change
#'   points, "^my_group_effect" plots all levels of a particular group-level effect, and
#'   "^cp_|^my_group_effect" plots both.
#' @param type String or vector of strings. Calls `bayesplot::mcmc_>>type<<()`.
#'   Common calls are "combo", "trace", and "dens_overlay". Current options include
#'   'acf', 'acf_bar', 'areas', 'areas_ridges', 'combo', 'dens', 'dens_chains',
#'   'dens_overlay', 'hist', 'intervals', 'rank_hist', 'rank_overlay', 'trace',
#'   'trace_highlight', and 'violin".
#' @param ncol Number of columns in plot. This is useful when you have many
#'   parameters and only one plot `type`.
#' @param nvariables Positive integer. Maximum number of parameters plotted per
#'   page. The default of 5 follows `brms::plot.brmsfit()`.
#' @param ask Logical. In an interactive session, prompt before printing each
#'   page after the first. Only used when there are multiple pages.
#' @param prior TRUE/FALSE. Plot using prior draws? Useful for `mcp(..., sample = "both")`
#'
#' @details
#'   For other `type`, it calls `bayesplot::mcmc_type()`. Use these
#'   directly on `coda::as.mcmc(fit)` or `as_draws(fit)` if you want finer
#'   control of plotting, e.g., `bayesplot::mcmc_dens(coda::as.mcmc(fit))`. There
#'   are also a number of useful plots in the \pkg{coda} package, i.e.,
#'   `coda::gelman.plot(coda::as.mcmc(fit))` and `coda::crosscorr.plot(coda::as.mcmc(fit))`
#'
#'   In any case, if you see a few erratic lines or parameter estimates, this is
#'   a sign that you may want to increase argument 'warmup' and 'iter' in \code{\link{mcp}}.
#'
#'   Up to `nvariables` parameters are shown on each page. Multi-page plots are
#'   printed in sequence; in interactive use, `ask = TRUE` pauses between pages.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @return A \pkg{ggplot2} object when all selected parameters fit on one page.
#'   For multiple pages, an invisible list of \pkg{ggplot2} objects.
#' @export
#' @seealso plot_dpar pp_check
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
#' # Useful for group-level effects:
#' # plot_pars(my_fit, pars = "group", ncol = 3)  # plot all group-level deviations
#' # plot_pars(my_fit, regex_pars = "my_group_effect", ncol = 3)  # one group-level effect
#' # pages = plot_pars(my_fit, pars = "group", ask = FALSE)
#' # pages[[1]]  # Inspect or customize one page
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
                     prior = FALSE,
                     nvariables = 5,
                     ask = TRUE
) {

  # Check arguments
  checkmate::assert_class(fit, "mcpfit")
  if (!is.character(pars) || !is.character(regex_pars))
    stop("`pars` and `regex_pars` has to be string/character.")
  check_legacy_parameter_names(pars, "plot_pars(pars)")

  if (!coda::is.mcmc.list(.subset2(fit, "mcmc_post")) && !coda::is.mcmc.list(.subset2(fit, "mcmc_prior")))
    stop("Cannot plot an mcpfit without prior or posterior draws.")

  if ("varying" %in% pars) {
    warning("`pars = \"varying\"` is deprecated; use `pars = \"group\"`.", call. = FALSE)
    pars[pars == "varying"] = "group"
  }

  if (any(c("population", "group") %in% pars) && length(pars ) > 1)
    stop("`pars` cannot combine 'population' or 'group' with other elements.")

  if (any(c("hex", "scatter") %in% type) && (length(pars) != 2 || length(regex_pars) > 0))
    stop("`type` = 'hex' or 'scatter' takes exactly two parameters which must be provided via the `pars` argument")

  if ("combo" %in% type && length(type) > 1)
    stop("'combo' type cannot be combined with other types. Replace 'combo' with the types you want combo\'ed")

  checkmate::assert_int(ncol, lower = 1)
  checkmate::assert_int(nvariables, lower = 1)
  checkmate::assert_flag(ask)
  checkmate::assert_flag(prior)
  if (any(c("hex", "scatter") %in% type) && nvariables < 2)
    stop("`nvariables` must be at least 2 for `type = 'hex'` or `type = 'scatter'`.")
  bayesplot::available_mcmc()  # Quick fix to make R CMD Check happy that bayesplot is imported

  # Get posterior/prior draws
  draws = posterior_draws(fit, prior = prior)

  # Handle special codes
  if ("population" %in% pars) {
    if (length(regex_pars) == 0) {
      pars = mcp_pars(fit, scope = "population")$name
    } else {
      # This probably means that the user left pars as default.
      pars = character(0)
    }
  } else if ("group" %in% pars) {
    # Regex search for group-level deviations
    regex_pars = paste0("^", mcp_pars(fit, scope = "group")$name, "\\[")
    pars = character(0)
  }

  select_parameters = utils::getFromNamespace("select_parameters", "bayesplot")
  pars = select_parameters(
    explicit = pars,
    patterns = regex_pars,
    complete_pars = posterior::variables(draws)
  )

  # Handles combo. Returns a customizable ggplot which "combo" does not.
  if ("combo" %in% type)
    type = c("dens_overlay", "trace")

  # Call the relevant bayesplot plot function for each type and page
  takes_facet = c(
    "areas", "dens", "dens_overlay", "trace", "hist", "intervals",
    "trace_highlight", "violin"
  )
  page_pars = split(pars, ceiling(seq_along(pars) / nvariables))
  pages = lapply(page_pars, function(this_page_pars) {
    all_plots = stats::setNames(lapply(type, function(this_type) {
      func = utils::getFromNamespace(paste0("mcmc_", this_type), "bayesplot")
      facet_args = if (this_type %in% takes_facet) list(ncol = ncol) else list()
      func(
        draws,
        pars = this_page_pars,
        regex_pars = character(0),
        facet_args = facet_args
      )
    }), type)

    patchwork::wrap_plots(all_plots) &
      ggplot2::theme(
        strip.placement = NULL,  # fixes bug: https://github.com/thomasp85/patchwork/issues/132
        legend.position = "none"  # no legend on chains. Takes up too much space
      )
  })

  if (length(pages) == 1)
    return(pages[[1]])

  for (page in seq_along(pages)) {
    if (page > 1 && ask && interactive())
      invisible(readline(prompt = "Press <Enter> to continue"))
    print(pages[[page]])
  }

  invisible(pages)
}
