`%>%` = magrittr::`%>%`  # instead of importing dplyr

#' Test a list of segments and simulation values
#'
#' * Simulates data from model and values
#' * Fits model to data
#' * Checks that values are recovered
#'
#' @aliases test_fit
#' @keywords internal
#' @param model A list of (unnamed) formulas
#' @param simulated Parameter values to be used for simulation.
#' @param newdata Optional simulation design. Defaults to 400 evenly spaced
#'   observations on `x`.
#' @param hyperparameters Optional generative parameters that are not arguments
#'   to `fit$simulate()`, such as group-level SDs.
#' @param family A family or `mcpfamily` used for both simulation and fitting.
#' @param series Optional independent-series column passed to `mcp()`.
#' @param chains Number of MCMC chains used when fitting.
#' @param warmup Number of warmup iterations used when fitting.
#' @param iter Number of post-adaptation iterations used when fitting.
#' @param seed Random seed for fitting. Defaults to 42.
test_fit = function(model, simulated, newdata = NULL, hyperparameters = NULL,
                     family = gaussian(), chains, warmup, iter,
                     min_ess, seed = 42) {
  if (Sys.getenv("MCP_TEST_LEVEL") != "release") {
    testthat::skip("Time-consuming fit recovery tests are only run when MCP_TEST_LEVEL='release'.")
  }

  # Simulate
  if (is.null(newdata)) {
    newdata = data.frame(
      x = seq(1, 200, length.out = 400),  # Needs to be reasonably high to get a correct estimate
      y = 0
    )
  } else if (!"y" %in% colnames(newdata)) {
    newdata$y = 0
  }
  empty = mcp(
    model, data = newdata, family = family, sample = FALSE,
    par_x = "x"
  )
  set.seed(42)
  simulated_y = suppressMessages(do.call(empty$simulate, c(list(fit = empty, newdata = newdata), simulated)))
  if (!is.null(hyperparameters)) {
    simulation_values = c(
      as.list(attr(simulated_y, "simulated")),
      hyperparameters
    )
    class(simulation_values) = c("mcplist", "list")
    attr(simulated_y, "simulated") = simulation_values
  }
  newdata[[mcp_columns(empty)$response]] = simulated_y

  # Fit
  fit = mcp(
    model, newdata, family = family, par_x = "x",
    chains = chains, warmup = warmup, iter = iter, seed = seed,
    diagnostics = FALSE, quiet = TRUE
  )  # Ensure convergence
  assign("fit", fit, envir = .GlobalEnv)  # for easier debugging

  test_matches_simulated(fit, min_ess)
}


#' Apply `test_fit` to each element of this list
#'
#' @aliases apply_test_fit
#' @keywords internal
#' @param all_models A list of lists. Each sub-list is an unnamed list of
#'   formulas with one named entry called "simulated" with parameter values to
#'   be used for simulation. It can optionally include a `family` entry.
#' @param family Default family used when a model does not supply one.
apply_test_fit = function(desc, all_models, family = gaussian()) {
  for (this in all_models) {
    # Split into formulas and simulation values
    simulated = this[["simulated"]]
    newdata = this[["newdata"]]
    hyperparameters = this[["hyperparameters"]]
    model_family = this[["family"]]
    chains = this[["chains"]]
    warmup = this[["warmup"]]
    iter = this[["iter"]]
    min_ess = this[["min_ess"]]
    seed = if (!is.null(this[["seed"]])) this[["seed"]] else 42
    if (is.null(model_family))
      model_family = family
    if (any(vapply(list(chains, warmup, iter, min_ess), is.null, logical(1))))
      stop("Every fit-recovery model must specify `chains`, `warmup`, `iter`, and `min_ess`.", call. = FALSE)
    model = this[names(this) == ""]

    # Test!
    testthat::test_that(desc, {
      test_fit(
        model, simulated, newdata, hyperparameters, model_family,
        chains, warmup, iter, min_ess, seed = seed
      )
    })
  }
}


#' Test whether posteriors matches simulated values
#'
#' Tests recovery of simulation parameters.
#' Tests effecitve N.
#'
#' @aliases test_matches_simulated
#' @keywords internal
#' @param fit An `mcpfit` object.
test_matches_simulated = function(fit, min_ess) {
  summaries = rbind(
    get_summary(fit, width = 0.97, scope = "population"),
    ranef(fit, width = 0.97)
  ) %>%
    dplyr::filter(is.na(sim) == FALSE)

  # Parameters recovery check:
  # - For K <= 3 parameters: strict 5% interval-width tolerance and 0 allowed failures.
  # - For K > 3 parameters: 20% interval-width tolerance and at most 1 allowed failure to prevent FWER false alarms.
  k = nrow(summaries)
  buffer_pct = if (k <= 3) 0.05 else 0.20
  max_failures = if (k <= 3) 0 else 1

  width = summaries$upper - summaries$lower
  new_lower = summaries$lower - buffer_pct * width
  new_upper = summaries$upper + buffer_pct * width
  matches = summaries$match == "OK" | (summaries$sim >= new_lower & summaries$sim <= new_upper)
  correctly_estimated = (sum(!matches) <= max_failures)

  # At least some effective samples for structural parameters. Individual
  # Group-level deviations are still checked for recovery above, but their ESS is not
  # a useful regression signal because it is dominated by partial pooling.
  ess_summaries = summaries[!grepl("[", summaries$variable, fixed = TRUE), ]
  good_eff = all(ess_summaries$ess_bulk > min_ess & ess_summaries$ess_tail > min_ess)

  # Test
  if (correctly_estimated == FALSE | good_eff == FALSE)
    print(summaries)
  testthat::expect_true(correctly_estimated)
  testthat::expect_true(good_eff)
}


#' Test the structure and rendering of an mcp_example() plot
#'
#' @keywords internal
#' @param plot The plot produced by `mcp_example()`.
#' @param example The name passed to `mcp_example()`.
test_example_plot = function(plot, example) {
  patchwork_examples = c("ar", "sigma")

  testthat::expect_s3_class(plot, "ggplot")
  testthat::expect_identical(inherits(plot, "patchwork"), example %in% patchwork_examples)
  testthat::expect_no_error(ggplot2::ggplotGrob(plot))

  if (example %in% c("group_mu", "group_cp"))
    testthat::expect_s3_class(plot$facet, "FacetWrap")
  if (example == "group_mu")
    testthat::expect_identical(plot$labels$colour, "condition")
  if (example == "multiple")
    testthat::expect_identical(plot$labels$colour, "group")
}


#' Snapshot an mcp_example() plot using an explicit theme
#'
#' @keywords internal
#' @param plot The plot produced by `mcp_example()`.
#' @param example The name passed to `mcp_example()`.
snapshot_example_plot = function(plot, example) {
  if (inherits(plot, "patchwork"))
    plot = plot & ggplot2::theme_gray()
  else
    plot = plot + ggplot2::theme_gray()

  vdiffr::expect_doppelganger(paste("mcp example", example), plot)
}


#' Test one mcp_example() fit and plot
#'
#' @keywords internal
test_mcp_example = function(example, snapshot = FALSE) {
  suppressMessages(capture.output({ fit = mcp_example(example, plot = TRUE) }))
  example_plot = ggplot2::last_plot()

  test_matches_simulated(fit, min_ess = 30)
  test_example_plot(example_plot, example)
  if (snapshot)
    snapshot_example_plot(example_plot, example)
}
