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
  testthat::skip_if(is.null(options("test_mcp_fits")[[1]]),
                    "This time-consuming test is only run locally before release.")

  # Simulate
  empty = mcp(model, sample = FALSE, par_x = "x")
  x = seq(1, 200, length.out = 400)
  data = data.frame(
    x = x,  # Needs to be reasonably high to get a correct estimate
    y = do.call(empty$simulate, c(list(x = x), simulated))
  )

  # Fit
  options(mc.cores = NULL)  # Respect `cores`
  quiet_out = purrr::quietly(mcp)(model, data, par_x = "x", chains = 4, cores = 4, adapt = 4000, iter = 2000)
  fit <<- quiet_out$result

  # Results table
  results_table = purrr::quietly(fixef)(fit, width = 0.98)$result
  recovered = all(results_table$match == "OK")  # Parameter recovery
  effective = all(results_table$n.eff > 100)  # Effective samples

  # Show table if the tests failed. Cannot be after tests for some reason...
  if (recovered == FALSE | effective == FALSE) {
    print(results_table)
  }

  # Tests
  testthat::expect_true(recovered, model)
  testthat::expect_true(effective, model)
}


#' Apply `test_fit` to each element of this list
#'
#' @aliases apply_test_fit
#' @keywords internal
#' @param all_models A list of lists. Each sub-list is an unnamed list of
#'   formulas with one named entry called "simulated" with parameter values to
#'   be used for simulation.
apply_test_fit = function(all_models, code) {
  for (this in all_models) {
    # Split into formulas and simulation values
    simulated = this[names(this) == "simulated"][[1]]
    model = this[names(this) == ""]

    # Test!
    testthat::test_that(
      test_fit(model, simulated),
      code = code
    )
  }
}
