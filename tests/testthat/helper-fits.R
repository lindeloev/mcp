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
test_fit = function(model, simulated) {
  testthat::skip("This time-consuming test is only run locally before release.")

  # Simulate
  newdata = data.frame(
    x = seq(1, 200, length.out = 400),  # Needs to be reasonably high to get a correct estimate
    y = rnorm(400)
  )
  empty = mcp(model, data = newdata, sample = FALSE, par_x = "x")
  newdata$y = do.call(empty$simulate, c(list(fit = empty, newdata = newdata), simulated))

  # Fit
  options(mc.cores = NULL)  # Respect `cores`
  quiet_out = purrr::quietly(mcp)(model, newdata, par_x = "x", chains = 5, cores = 5, adapt = 10000, iter = 3000)  # Ensure convergence
  fit <<- quiet_out$result  # assign to global namespace for easier debugging

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
    simulated = this[names(this) == "simulated"][[1]]
    model = this[names(this) == ""]

    # Test!
    testthat::test_that(desc, {test_fit(model, simulated)})
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
  good_eff = all(summaries$n.eff > 50)

  # Test
  if (correctly_estimated == FALSE | good_eff == FALSE)
    print(summaries)
  testthat::expect_true(correctly_estimated)
  testthat::expect_true(good_eff)
}
