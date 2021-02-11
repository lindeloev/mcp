# ABOUT: These functions take data/samples and run them through
# formula code to make model predictions
# ------------------------------------------------------------

#' Make the levels of categorical predictors match the original data
#' @aliases relevel_newdata
#' @keywords internal
#' @inheritParams add_rhs_predictors
#' @return `newdata` with all categorical columns as factors that have identical levels to the original data.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
relevel_newdata = function(newdata, fit) {
  # Check that the necessary data is available
  rhs_vars = get_rhs_vars(fit$model)
  assert_data_cols(rhs_vars, newdata, fail_types = c("na", "nan", "infinite"))

  # Make sure to carry over the exact level-structure of the original model
  for (col_name in rhs_vars) {
    org_col = fit$data[, col_name]
    new_col = newdata[, col_name]
    if (is.character(org_col) | is.factor(org_col)) {
      new_col = factor(new_col, levels = levels(factor(org_col)))

      # Helpful error
      if (any(is.na(new_col))) {
        new_levels = newdata[is.na(new_col), col_name] %>% as.character() %>% unique()
        stop("Got novel values (", and_collapse(new_levels), ") for column ", col_name, ". Only values used during fitting are allowed.")
      }
    }

    newdata[, col_name] = new_col
  }

  return(newdata)
}


#' Add predictors to `newdata`
#' @aliases add_rhs_predictors
#' @keywords internal
#' @param newdata A data.frame that contains the RHS-predictors from `fit$data` and `fit$model`.
#' @param fit An `mcpfit` object.
#' @return `newdata` with additional (dummy-coded) columns that make up the design matrix required for the model.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
add_rhs_predictors = function(newdata, fit) {
  assert_types(fit, "mcpfit")
  assert_types(newdata, "data.frame", "tibble")

  # Make categorical predictors match the original data
  newdata = newdata %>%
    as.data.frame() %>%
    relevel_newdata(fit)

  # Get rhs_matrix
  rhs_table_tmp = get_rhs_table(fit$model, newdata, fit$family, fit$pars$x, check_rank = FALSE)
  rhs_matrix = get_rhs_matrix(rhs_table_tmp)

  # Check that the variable structure matches
  if (nrow(fit$.internal$rhs_table) != nrow(rhs_table_tmp) || any(fit$.internal$rhs_table$code_name != rhs_table_tmp$code_name))
    stop_github("rhs_table_tmp does not match fit$.internal$rhs_table.")

  # All permutations of rows in newdata and parameters
  predictors = as.data.frame(rhs_matrix) %>%
    magrittr::set_colnames(paste0(".pred_", colnames(rhs_matrix))) %>%
    dplyr::bind_cols(newdata)

  return(predictors)
}


#' Returns parameters needed for simulation
#' @aliases get_sim_pars
#' @keywords internal
#' @return Character vector
get_sim_pars = function(rhs_table, pars) {
  sim_pars = c(
    pars$reg[stringr::str_starts(pars$reg, "cp_")],  # cp_*
    rhs_table$code_name  # mu, sigma, ar, etc.
  )

  # Default to varying = 0
  if (is.null(pars$varying) == FALSE)
    sim_pars = c(sim_pars, paste0(pars$varying, " = 0"))

  return(sim_pars)
}


#' Vectorized R-side run of the full model.
#'
#' Used internally in mcp.
#'
#' @aliases simulate_vectorized
#' @keywords internal
#' @inheritParams mcp
#' @inheritParams simulate.mcpfit
#' @param ... Parameter names (e.g., `cp_1 = c(4.3, 4.5, 6.2), Intercept_1 = c(11.2, 12.1, 10.9)`, etc.)
#'   and predictor columns from `rhs_tables$matrix_data` prefixed with ".pred_" (e.g., `.pred_Intercept_1 = c(1, 1, 1)`).
#' @return Vector with same length as inputs.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
simulate_vectorized = function(fit, ..., .type = "predict", .rate = FALSE, .which_y = "mu", .arma = TRUE, .scale = "response") {
  assert_types(fit, "mcpfit")
  rhs_table = fit$.internal$rhs_table  # Shorthand

  # Assert that the ellipsis contains the expected argument names
  param_pars = get_sim_pars(rhs_table, fit$pars)
  pred_pars = paste0(".pred_", rhs_table$code_name)
  data_pars = c(fit$pars$x, fit$pars$trials)
  expected_arg_names = c(param_pars, pred_pars, data_pars)

  args = list(...)
  missing_args = dplyr::setdiff(expected_arg_names, names(args))
  if (length(missing_args) > 0)
    stop_github("Missing the following arguments: ", and_collapse(missing_args))

  # Other args
  assert_value(.type, allowed = c('predict', 'fitted'), len = 1)
  assert_logical(.rate, len = 1)
  assert_value(.which_y, allowed = paste0(rhs_table$dpar, tidyr::replace_na(rhs_table$order, "")))  # "mu" "ar1" "ar2" "sigma", etc.
  assert_logical(.arma, len = 1)
  assert_value(.scale, allowed = c('response', 'linear'), len = 1)
  if (.scale == 'linear' && .type == 'predict')
    stop("Only `type = 'fitted'` is meaningful when `scale = 'linear'`")


  #############################
  # BUILD AND EXECUTE FORMULA #
  #############################
  # Get model parameters from ellipsis (args)
  #param_args = args[names(args) %in% param_pars]
  #data_args = args[names(args) %in% data_pars]

  # rhs_matrix_ from ellipsis.
  pred_args = args[names(args) %in% pred_pars]
  rhs_matrix_ = do.call(cbind, pred_args)
  rhs_matrix_ = rhs_matrix_[, match(pred_pars, colnames(rhs_matrix_))]  # Same order as rhs_table$code_name no matter order of args

  # Run it!
  cp_0 = -Inf
  assign(paste0("cp_", length(fit$model)), Inf)  # e.g., cp_3 = Inf
  eval(parse(text = fit$.internal$formula_r))


  ################
  # RETURN STUFF #
  ################
  .which_y = paste0(.which_y, '_')
  if (.scale == 'response') {
    assign(.which_y, fit$family$linkinv(get(.which_y)))
  }

  # Simply return for fitted non-mu params
  if (.which_y != 'mu_' & .type == 'fitted') {
    return(get(.which_y))
  } else if (.which_y != 'mu_' & .type != 'fitted') {
    stop("`type = 'fitted'` is the only option when `which_y != 'mu'`")
  }

  # Shortcuts for readability
  is_arma = any(rhs_table$dpar == "ar")

  # Return functions here
  if (fit$family$family == "gaussian") {
    # If ARMA, build resid_ and return with that
    if (is_arma && .arma == TRUE) {
      ar_code = get_ar_code(
        ar_order = get_arma_order(fit$pars$arma),
        family = fit$family,
        is_R = TRUE,
        xvar = fit$pars$x,
        yvar = fit$pars$y
      )

      eval(parse(text = ar_code))
    }

    # If fitted or no data
    if (.type == 'fitted') {
      return(mu_)
    } else if (.type == 'predict') {
      if (any(fit$family$linkinv(sigma_) < 0))
        stop("Modelled negative sigma. First detected at ", fit$pars$x, " = ", min(get(fit$pars$x)[fit$family$linkinv(sigma_) < 0]))

      # Complex code if ARMA. Simple if not. resid_sigma_ was generated from sigma_ so no need to do an extra rnorm().
      if (is_arma == TRUE) {
        if (arma == TRUE) {
          return(fit$family$linkinv(fit$family$linkfun(mu_ + resid_sigma_)))
        } else {
          return(family$linkinv(fit$family$linkfun(mu_) + rnorm(length(mu_), 0, sigma_)))
        }
      } else {
        return(rnorm(length(mu_), mu_, sigma_))
      }
    }

    # OTHER FAMILIES ---------------------
  } else if (fit$family$family == "binomial") {
    if (.type == 'predict') {
      if (.rate == FALSE) return(rbinom(length(mu_), get(pars$trials), mu_))
      if (.rate == TRUE)  return(rbinom(length(mu_), get(pars$trials), mu_) / get(pars$trials))
    } else if (.type == 'fitted') {
      if (.rate == FALSE) return(", pars$trials, " * mu_)
      if (.rate == TRUE)  return(mu_)
    }
  } else if (fit$family$family == "bernoulli") {
    if (.type == 'predict') return(rbinom(length(mu_), 1, mu_))
    if (.type == 'fitted') return(mu_)
  } else if (fit$family$family == "poisson") {
    if (.type == 'predict') {
      if ((.scale == 'response' && any(mu_ > 2146275819)) || (.scale == 'linear' && any(fit$family$linkinv(mu_) > 2146275819)))
        stop("Modelled extremely large value: ", fit$family$linkinv_str, "(", fit$pars$y, ") > 2146275819.")
      return(rpois(length(mu_), mu_))
    } else if (.type == 'fitted') {
      return(mu_)
    }
  } else if (fit$family$family == "exponential") {
    if (.type == 'predict') return(rexp(length(mu_), mu_))
    if (.type == 'fitted') return(1 / mu_)
  }
}


#' User-friendy interface to simulate_vectorized
#'
#' Does the following:
#'  * Converts factors in `newdata` to dummies with correct levels.
#'  * Checks and binds parameter values (`...`)
#'  * Call `simulate_vectorized`
#'  * Add "simulated" attribute and returns
#'
#' @aliases simulate_atomic
#' @keywords internal
#' @inheritParams mcp
#' @inheritParams simulate.mcpfit
#' @param newdata A data.frame or tibble that should contain the same columns as `fit$data`
#' @param ... Parameter values of length 1, e.g., `cp_1 = 80, Intercept_1 = -22.5` etc.
#' @return Numeric vector with attribute "simulated".
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
simulate_atomic = function(fit,
                        newdata,
                        ...,
                        .type = 'predict',
                        .rate = FALSE,
                        .which_y = 'mu',
                        .arma = TRUE,
                        .scale = 'response') {

  # Check inputs.
  mcp:::assert_types(fit, "mcpfit")
  mcp:::assert_types(newdata, "data.frame", "tibble")
  args = list(...)
  expected_args = c(
    fit$pars$reg[stringr::str_starts(fit$pars$reg, "cp_")],  # cp_*
    fit$.internal$rhs_table$code_name  # mu, sigma, ar, etc.
  )
  if (is.null(names(args)) | any(names(args) == ""))  # A more helpful error than assert_ellipsis()
    stop("All arguments must be named.")
  mcp:::assert_ellipsis(..., allowed = expected_args)
  lapply(args, mcp:::assert_numeric, len = 1)
  # Other values are asserted in func_fast

  # Get permutations
  predictors = add_rhs_predictors(newdata, fit)
  pred_param_grid = cbind(predictors, args)  # Use tidyr::expand_grid() if any args have length > 1.

  # Get y
  return_me = pred_param_grid %>%
    dplyr::mutate(
      .return_me = rlang::exec(simulate_vectorized,
                               fit,
                               !!!.,
                               .type = .type,
                               .rate = .rate,
                               .which_y = .which_y,
                               .arma = .arma,
                               .scale = .scale)
    ) %>%
    dplyr::pull(.return_me)

  # add_simulated etc. and return
  attr(return_me, 'simulated') = args  # Set as attribute
  class(attr(return_me, 'simulated')) = c("mcplist", "list")  # for nicer printing
  return(return_me)
}


#' Wrapper for `simulate_atomic()` with named arguments
#'
#' Auto-completion of named arguments  makes it easier to call from, e.g., RStudio.
#' This is typically stored as `fit$simulate()`.
#'
#' @aliases get_fitsimulate
#' @keywords internal
#' @param pars A list of model parameters, typically from `fit$pars`
#' @return An R function.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_fitsimulate = function(pars) {
  # List of argument names
  pars$reg = pars$reg[!stringr::str_ends(pars$reg, "_sd")]  # Remove hyperparameter on varying effects from pars$reg since it is not used for simulation
  y_args = ifelse(length(pars$arma > 0), paste0(pars$y, " = NULL, .ydata = NULL"), "")
  varying_args = ifelse(length(pars$varying) > 0, paste0(pars$varying, " = 0"), "")
  required_args = c(pars$reg, pars$sigma, pars$arma, y_args, varying_args)
  required_args = required_args[required_args != ""]  # remove empty strings

  # Build function
  fitsimulate_code = paste0("function(fit, newdata, ", paste0(required_args, collapse = ", "), ",
  .type = 'predict',
  .rate = FALSE,
  .which_y = 'mu',
  .arma = TRUE,
  .scale = 'response') {

  result = simulate_atomic(fit, newdata, ", paste0(required_args, " = ", required_args, collapse = ", "), ", .type = .type, .rate = .rate, .which_y = .which_y, .arma = .arma, .scale = .scale)
  return(result)
}")

  fitsimulate_func = eval(parse(text = fitsimulate_code))
  return(fitsimulate_func)
}
