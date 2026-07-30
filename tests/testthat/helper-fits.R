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
test_fit = function(model, simulated, newdata = NULL, hyperparameters = NULL) {
  # COMMENT THIS LINE TO ENABLE EXTENSIVE TESTS. THEN REVERT. DO NOT COMMIT.
  # It's a slow test so we do it rarely.
  testthat::skip("This time-consuming test is only run locally before release.")

  # Simulate
  if (is.null(newdata)) {
    newdata = data.frame(
      x = seq(1, 200, length.out = 400),  # Needs to be reasonably high to get a correct estimate
      y = rnorm(400)
    )
  }
  empty = mcp(model, data = newdata, sample = FALSE, par_x = "x")
  simulated_y = do.call(empty$simulate, c(list(fit = empty, newdata = newdata), simulated))
  if (!is.null(hyperparameters)) {
    simulation_values = c(
      as.list(attr(simulated_y, "simulated")),
      hyperparameters
    )
    class(simulation_values) = c("mcplist", "list")
    attr(simulated_y, "simulated") = simulation_values
  }
  newdata$y = simulated_y

  # Fit
  quiet_out = purrr::quietly(mcp)(model, newdata, par_x = "x", chains = 5, adapt = 10000, iter = 3000)  # Ensure convergence
  fit = quiet_out$result
  assign("fit", fit, envir = .GlobalEnv)  # for easier debugging

  test_matches_simulated(fit)
}


#' Apply `test_fit` to each element of this list
#'
#' @aliases apply_test_fit
#' @keywords internal
#' @param all_models A list of lists. Each sub-list is an unnamed list of
#'   formulas with one named entry called "simulated" with parameter values to
#'   be used for simulation.
apply_test_fit = function(desc, all_models) {
  for (this in all_models) {
    # Split into formulas and simulation values
    simulated = this[["simulated"]]
    newdata = this[["newdata"]]
    hyperparameters = this[["hyperparameters"]]
    model = this[names(this) == ""]

    # Test!
    testthat::test_that(desc, {
      test_fit(model, simulated, newdata, hyperparameters)
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
test_matches_simulated = function(fit) {
  summaries = rbind(
    fixef(fit, width = 0.97),
    ranef(fit, width = 0.97)
  ) %>%
    dplyr::filter(is.na(sim) == FALSE)

  # Parameters within lower/upper + 10%
  new_lower = summaries$lower - 0.1*(summaries$mean - summaries$lower)
  new_upper = summaries$upper - 0.1*(summaries$mean - summaries$upper)
  correctly_estimated = all(summaries$match == "OK" | (summaries$sim > new_lower & summaries$sim < new_upper))

  # At least some effective samples
  good_eff = all(summaries$ess_bulk > 50 & summaries$ess_tail > 50)

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
  patchwork_examples = c("ar", "variance")
  titles = list(
    ar = c("plot(fit)", 'plot_dpar(fit, "ar1")'),
    binomial = "plot(fit)",
    demo = "plot(fit)",
    group = 'plot(fit, facet_by = "participant", color_by = "condition")',
    intercepts = "plot(fit)",
    multiple = 'plot(fit, color_by = "group")',
    quadratic = "plot(fit)",
    variance = c('plot(fit, q_predict = TRUE)', 'plot_dpar(fit, "sigma")'),
    varying = 'plot(fit, facet_by = "id")'
  )

  testthat::expect_s3_class(plot, "ggplot")
  testthat::expect_identical(inherits(plot, "patchwork"), example %in% patchwork_examples)
  testthat::expect_no_error(ggplot2::ggplotGrob(plot))

  plot_titles = if (inherits(plot, "patchwork")) {
    vapply(seq_along(plot), function(i) plot[[i]]$labels$title, character(1))
  } else {
    plot$labels$title
  }
  testthat::expect_identical(unname(plot_titles), titles[[example]])

  if (example %in% c("group", "varying"))
    testthat::expect_s3_class(plot$facet, "FacetWrap")
  if (example == "group")
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
