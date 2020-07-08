#' Call `get_formula_str` for multiple ytypes and paste strings
#'
#' Currently used to differentiate between the JAGS model (use all) and the
#' fit$simulate model (do not include arma).
#'
#' @aliases get_all_formulas
#' @keywords internal
#' @inheritParams get_formula_str
#' @inheritParams mcp
#' @param ytypes A character vector of ytypes to including in model building
#' @return A string with JAGS code.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_all_formulas = function(ST, prior, par_x, ytypes = c("ct", "sigma", "arma")) {
  # Initiate with central tendency
  if ("ct" %in% ytypes) {
    formula_str = get_formula_str(ST, par_x, ytype = "ct", init = TRUE)
  }

  # Add sigma
  if ("sigma" %in% ytypes) {
    if (any(stringr::str_starts(names(prior), "sigma_"))) {
      formula_str = paste0(formula_str, "\n\n",
                           get_formula_str(ST, par_x, ytype = "sigma"))
    }
  }

  # Add arma
  if ("arma" %in% ytypes) {
    all_arma = names(prior)[stringr::str_starts(names(prior), "(^ar|^ma)[0-9]")]  # All priors starting with ar[number] or ma[number]
    arma_bases = unique(sub("*([0-9]+)_.*$", "\\1", all_arma))  # extract these starts
    for (arma_base in arma_bases) {
      formula_str = paste0(formula_str, "\n\n",
                           get_formula_str(ST, par_x, ytype = arma_base))
    }
  }

  # Return
  return(formula_str)
}


#' Build an R formula (as string) given a segment table (ST)
#'
#' You will need to replace PAR_X for whatever your x-axis observation column
#' is called. In JAGS typically `x[i_]`. In R just `x`.
#'
#' @aliases get_formula_str
#' @keywords internal
#' @inheritParams mcp
#' @param ST Tibble. Returned by `get_segment_table`.
#' @param ytype One of "ct" (central tendency), "sigma", "ar1" (or another order), or "ma1" (or another order)
#' @param init TRUE/FALSE. Set to TRUE for the first call. Adds segment-relative
#'   X-codings and verbose commenting of one formula
#' @return A string with JAGS code.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_formula_str = function(ST, par_x, ytype = "ct", init = FALSE) {
  # Build this! Start empty...
  formula_str = ""

  ############
  # Initiate #
  ############
  # Optionally add X-helpers which code the X relative to the start of each segment.
  # Used to compute slopes so that they start in the beginning of each segment.
  if (init == TRUE) {
    for (i in seq_len(nrow(ST))) {
      S = ST[i, ]
      segment_start = ifelse(i > 1, yes = paste0(" - ", S$cp_code_form, ""), no = "")  #
      segment_end = ifelse(i < nrow(ST), yes = paste0(ST$cp_code_form[i + 1], ""), no = paste0("cp_", i))  # infinite if last segment.

      formula_str = paste0(formula_str, "\nX_", i, "_[i_] = min(", par_x, "[i_], ", segment_end, ")", segment_start)
    }
  }

  for (i in seq_len(nrow(ST))) {
    #################################
    # GET ST COLUMNS FOR THIS YTYPE #
    #################################
    # Extract int, slope_table, and slope_code from ST given ytype.
    # This is possible since the generated formula works the same for all of them.
    # Programmer's note: The code below is quite repetitive, but I found no simpler solutions.
    S = ST[i, ]
    if (ytype == "ct") {
      int = S$ct_int[[1]]
      slope_table = S$ct_slope[[1]]
      slope_code = S$ct_code[[1]]
    } else if (ytype == "sigma") {
      int = S$sigma_int[[1]]
      slope_table = S$sigma_slope[[1]]
      slope_code = S$sigma_code[[1]]
    } else if (stringr::str_starts(ytype, "ma")) {
      # Moving average: more involved
      ma_order = get_arma_order(ytype)  # e.g., "ar12" --> 12

      # Get int (intercept)
      S_int_order = length(S$ma_int[[1]])  # How many orders are recorded in the segment table
      if (ma_order <= S_int_order) {
        int = S$ma_int[[1]][[ma_order]]
      } else {
        int = NA
      }

      # Get slope
      S_slope_order = length(S$ma_slope[[1]])  # How many orders are recorded in the segment table
      if (ma_order <= S_slope_order) {
        slope_table = S$ma_slope[[1]][[ma_order]]
        slope_code = S$ma_code[[1]][[ma_order]]
      } else {
        slope_table = NA
        slope_code = NA
      }
    } else if (stringr::str_starts(ytype, "ar")) {
      # Autoregressive: more involved
      ar_order = get_arma_order(ytype)  # e.g., "ar12" --> 12

      # Get int (intercept)
      S_int_order = length(S$ar_int[[1]])  # How many orders are recorded in the segment table
      if (ar_order <= S_int_order) {
        int = S$ar_int[[1]][[ar_order]]
      } else {
        int = NA
      }

      # Get slope
      S_slope_order = length(S$ar_slope[[1]])  # How many orders are recorded in the segment table
      if (ar_order <= S_slope_order) {
        slope_table = S$ar_slope[[1]][[ar_order]]
        slope_code = S$ar_code[[1]][[ar_order]]
      } else {
        slope_table = NA
        slope_code = NA
      }
    } else {
      stop("Got wrong argument for ytype: ", ytype, ". Report an issue on GitHub if you see this.")
    }


    ########################
    # BUILD FORMULA STRING #
    ########################

    # Start formula string with an informative comment
    if (i == 1) {
      if (ytype == "ct") {
        formula_str = paste0(formula_str, "\n\n# Fitted value\ny_[i_] = \n")
      } else if (ytype == "sigma") {
        formula_str = paste0(formula_str, "# Fitted standard deviation\nsigma_[i_] = max(0, \n")  # Add max(0, [formula]) to prevent modeling negative sigmas. JAGS uses precision = 1 / sigma^2 which yields positive precisions for negative sigmas.
      } else if (stringr::str_detect(ytype, "ar[0-9]+")) {
        formula_str = paste0(formula_str, "# Autoregressive coefficient for all AR(", ar_order,")\n", ytype, "_[i_] = \n")
      } else if (stringr::str_detect(ytype, "ma[0-9]+")) {
        formula_str = paste0(formula_str, "# Moving Average: MA(", ma_order, ")\n", ytype, "_[i_] = \n")
      }
    }

    # FUTURE_REL will sometimes be replaced by a less-than indicator (ind_past).
    # cp_code_form includes varying effects
    ind_this = paste0("  (", par_x, "[i_] >= ", S$cp_code_form, ") * FUTURE_REL")
    ind_past = paste0("(", par_x, "[i_] < ", S$cp_code_form, ") * ")

    # Begin building formula
    if (init == TRUE) {  # Verbose for mean formula
      formula_str = paste0(formula_str, "\n  # Segment ", i, ": ", S$form, "\n")
    }


    # Add intercept
    if (!all(is.na(int))) {
      # For absolute intercepts, "remove" earlier stuff affecting the intercept
      # Multiply it with zero from this change point and on
      if (int$rel == FALSE) {
        formula_str = gsub("FUTURE_REL", ind_past, formula_str)
      }

      # Add intercept with indicator
      formula_str = paste0(formula_str, ind_this, int$name, " + \n")
    }

    # Add slope
    if (!all(is.na(slope_table))) {
      formula_str = paste0(formula_str, ind_this, slope_code, " + \n")
    }

    # If this is just a plateau (~ 0), i.e., the absence of intercept, slope, and varying effects
    if (init == TRUE) {  # Verbose for mean formula. Not for sigma.
      if (all(is.na(int)) & all(is.na(slope_table))) {
        formula_str = paste0(formula_str, "  0 + \n")
      }
    }

    # Finish up formula_str for the last segment
    if (i == nrow(ST)) {
      formula_str = substr(formula_str, 1, nchar(formula_str) - 3)
    }
  }

  # Remove all "unrealized" FUTURE_REL
  formula_str = gsub("FUTURE_REL", "", formula_str)

  if (ytype == "sigma") {
    formula_str = paste0(formula_str, ")")  # End max() statement
  }

  # Return
  formula_str
}


#' Turn formula_str into a proper R function
#'
#' @aliases get_simulate
#' @keywords internal
#' @inheritParams mcp
#' @param formula_str string. Returned by `get_formula`.
#' @param par_trials String. For binomial models: name of trials column.
#' @param pars_pop List of population parameters which the user must provide.
#' @param pars_varying List varying parameters. They will default to zero
#'   (optional for the user).
#' @param nsegments Positive integer. Number of segments, typically `nrow(ST)`.
#' @return A string with R code for the fit$simulate() function corresponding to the model.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_simulate = function(formula_str, pars, nsegments, family) {
  # First some substitutions
  formula_func = gsub("\\[i_\\]", "", formula_str)  # No explicit indexing needed for R function
  formula_func = gsub("min\\(", "pmin\\(", formula_func)  # vectorized min
  formula_func = gsub("max\\(", "pmax\\(", formula_func)  # vectorized max
  for (i in seq_len(nsegments)) {
    formula_func = gsub(paste0("CP_", i, "_INDEX"), "", formula_func)  # Use vector of data
  }

  # Remove hyperparameter on varying effects from pars$reg since it is not used for simulation
  pars$reg = pars$reg[!stringr::str_ends(pars$reg, "_sd")]

  # Helper
  is_arma = length(pars$arma) > 0

  # Helper to build the list of simulate() arguments. Simply comma separates correctly
  args_if_exists = function(args, postfix = "") {
    if(length(args) > 0) {
      args_str = paste0(args, postfix, collapse = ", ")
      return(paste0(args_str, ", "))
    }
    else {
      return("")
    }
  }

  # Now build the function R code
  # x_and_trials handles the special case that binomial "trials" is also used as a predictor.
  x_and_trials = ifelse(is.null(pars$trials),
                        yes = pars$x,  # just the x
                        no = paste0(unique(c(pars$x, pars$trials)), collapse = ", "))  # x and N, if they differ.
  out = paste0("
  function(",
    x_and_trials, ", ",
    args_if_exists(pars$reg),
    args_if_exists(pars$sigma),
    args_if_exists(pars$arma),
    args_if_exists(pars$varying, " = 0"),
    ifelse(is_arma, paste0("\n    ", pars$y, " = NULL, \n    .ydata = NULL, "), ""), "
    type = 'predict',
    quantile = FALSE,
    rate = FALSE,
    which_y = 'ct',
    add_attr = TRUE,
    arma = TRUE,
    ...) {

    # Return predictions or fitted values?
    if (!type %in% c('predict', 'fitted'))
      stop(\"'`type` must be one of 'predict' or 'fitted'\")

    # Use this for return to add simulation parameters to the output
    args_names = as.list(match.call())  # Which arguments this function was called with
    all_values = c(as.list(environment()), list(...))  # all vars and values in env
    add_simulated = function(x) {
      # Do not add simulated attribute
      if (add_attr == FALSE) {
        return(x)
      } else {
        # Add it!
        args_values = all_values[names(all_values) %in% names(args_names)]  # Only those coming from call
        args_values[['", pars$x ,"']] = NULL  # Remove x
        ", ifelse(length(pars$trials) > 0, yes = paste0("args_values[['", pars$trials, "']] = NULL  # Remove trials"), no = ""), "
        attr(x, 'simulated') = args_values  # Set as attribute
        return(x)
      }
    }

    # Helpers to simplify making the code for this function
    cp_0 = -Inf
    cp_", nsegments, " = Inf

    ", formula_func)

  # Return depends on family
  # GAUSSIAN ------------------------------
  if (family$family == "gaussian") {
    # If ARMA, build resid_ and return with that
    if (is_arma) {
      out = paste0(out, "

      if (arma == TRUE) {
        # Simulate AR residuals",
          get_ar_code(
            ar_order = get_arma_order(pars$arma),
            family = family,
            is_R = TRUE,
            xvar = pars$x,
            yvar = pars$y
          ), "
        # Finally, add these autocorrelated predictions to the fitted y_:
        y_ = y_ + resid_
      }
    ")
    }

    out = paste0(out, "
    # Use which_y to return something else than ct.
    # First, rename to internal parameter name
    if (which_y == 'sigma') {
      which_y = 'sigma_'
    } else if (stringr::str_detect(which_y, '^ar([0-9]+)$')) {
      which_y = paste0(which_y, '_')
    }
    if (which_y != 'ct') {
      if (!exists(which_y))
        stop(which_y, ' was not found. Try one of \\'sigma\\', \\'ar1\\', etc.')
      if (type != 'fitted')
        stop('`type = \\'fitted\\` is the only option when `which_y != \\'ct\\'`')

      return(add_simulated(get(which_y)))
    }")

    if (is_arma) {
      out = paste0(out, "
   if (is.null(", pars$y, ") & is.null(.ydata) & arma == TRUE)
     message('Returning without AR(N) since `", pars$y, "` was not in the data.')")
    }

    # If fitted or no data
    out = paste0(out, "
    if (type == 'fitted') {
      return(add_simulated(", family$linkinv_r, "(y_)))
    }")
    out = paste0(out, " else if (type == 'predict') {
      if (any(", family$linkinv_r, "(sigma_) < 0))
        stop('Modelled negative sigma. First detected at ", pars$x, " = ', min(", pars$x, "[", family$linkinv_r, "(sigma_) < 0]))
      return(add_simulated(rnorm(length(", pars$x, "), ", family$linkinv_r, "(y_), sigma_)))
    }")

    # return(add_simulated(qnorm(length(", pars$x, "), ", family$linkinv_r, "(y_), sigma_)))
    # if (is.numeric(quantile)) {
    #   add_simulated(qnorm(quantile, ", family$linkinv_r, "(y_), sigma_))
    # } else if (quantile == FALSE) {
    #   add_simulated(rnorm(length(", pars$x, "), ", family$linkinv_r, "(y_), sigma_))
    # } else {
    #   stop('Invalid `quantile` argument to simulate()')
    # }

  # OTHER FAMILIES ---------------------
  } else if (family$family == "binomial") {
    out = paste0(out, "
    if (type == 'predict') {
      if (rate == FALSE) return(add_simulated(rbinom(length(", pars$x, "), ", pars$trials, ", ", family$linkinv_r, "(y_))))
      if (rate == TRUE)  return(add_simulated(rbinom(length(", pars$x, "), ", pars$trials, ", ", family$linkinv_r, "(y_)) / ", pars$trials, "))
    } else if (type == 'fitted') {
      if (rate == FALSE) return(add_simulated(", pars$trials, " * ", family$linkinv_r, "(y_)))
      if (rate == TRUE)  return(add_simulated(", family$linkinv_r, "(y_)))
    }")
  } else if (family$family == "bernoulli") {
    out = paste0(out, "
    if (type == 'predict') return(add_simulated(rbinom(length(", pars$x, "), 1, ", family$linkinv_r, "(y_))))
    if (type == 'fitted') return(add_simulated(", family$linkinv_r, "(y_)))")
  } else if (family$family == "poisson") {
    out = paste0(out, "
    if (type == 'predict') {
      if (any(", family$linkinv_r, "(y_) > 2146275819))
        stop('Modelled extremely large value: ", family$linkinv_r, "(", pars$y, ") > 2146275819. First detected at ", pars$x, " = ', min(", pars$x, "[", family$linkinv_r, "(y_) > 2146275819]))
      return(add_simulated(rpois(length(", pars$x, "), ", family$linkinv_r, "(y_))))
    } else if (type == 'fitted') {
        return(add_simulated(", family$linkinv_r, "(y_)))
      }")
  } else if (family$family == "exponential") {
    out = paste0(out, "
    if (type == 'predict') return(add_simulated(rexp(length(", pars$x, "),", family$linkinv_r, "(y_))))
    if (type == 'fitted') return(add_simulated(1 / ", family$linkinv_r, "(y_)))
      ")
  }

  out = paste0(out, "
  }")

  # Return function
  eval(parse(text = out))
}



#' Gets code for ARMA terms, resulting in a "resid_"
#'
#' Developer note: Ensuring that this can be used in both simulate() and JAGS
#' got quite messy with a lot of if-statements. It works but some refactoring
#' may be good in the future.
#'
#' @aliases get_ar_code
#' @keywords internal
#' @param ar_order Positive integer. The order of ARMA
#' @param family An mcpfamily object
#' @param is_R Bool. Is this R code (TRUE) or JAGS code (FALSE)?
#' @return String with JAGS code for AR.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_ar_code = function(ar_order, family, is_R, xvar, yvar = NA) {
  mm = "\n"

  ##################
  # FOR SIMULATE() #
  ##################
  if (is_R) {
    # mm is the code string to be populated below
    mm = paste0("
      if(!is.null(.ydata)) ", yvar, " = .ydata  # Hack to get a consistent argument name
      ar0_ = sigma_[1:", ar_order, "] * 0 + 1

      # If got y. Use it to compute residuals
      if (!all(is.na(", yvar, "))) {
        if (!is.numeric(", yvar, "))
          stop('Wrong format of ", yvar, ". Should be numeric.')
        resid_sigma_ = ", yvar, " - y_
      } else {
        # No y, simulate residuals
        resid_sigma_ = rnorm(length(", xvar, "), 0, sigma_)
      }

      resid_ = numeric(length(", xvar, "))")

    # For data points lower than the full order
    for (i in seq_len(ar_order)) {
      mm = paste0(mm, "\n      resid_[", i, "] = ",
                  paste0("ar", 0:(i-1), "_[", i, " - ", 0:(i-1), "] * resid_sigma_[", i, " - ", 0:(i-1), "]", collapse = " + "))
    }

    # For full order
    mm = paste0(mm, "
      for (i_ in ", ar_order + 1, ":length(", xvar, ")) {
        if (all(is.na(", yvar, "))) {

          resid_[i_] = resid_sigma_[i_]")
    for (i in seq_len(ar_order)) {
      mm = paste0(mm, " + ar", i, "_[i_] * resid_[i_ - ", i, "]")
    }
    mm = paste0(mm, "
        } else {
          resid_[i_] = ")
    for (i in seq_len(ar_order)) {
      mm = paste0(mm, " + ar", i, "_[i_] * resid_sigma_[i_ - ", i, "]")
    }

    mm = paste0(mm, "
        }")

    # Finish up
    mm = paste0(mm, "\n
      }")
  } else {


    #################
    # FOR JAGS CODE #
    #################
    # For data points lower than the full order
    mm = paste0(mm, "
  # Apply autoregression to the residuals
  resid_[1] = 0")

    # For data points lower than the full order
    if (ar_order >= 2) {
      for (i in 2:ar_order) {
        mm = paste0(mm, "
  resid_[", i, "] = ", paste0("ar", 1:(i-1), "_[", i, " - ", 1:(i-1), "] * resid_sigma_[", i, " - ", 1:(i-1), "]", collapse = " +\n              "))
      }
    }

    # For full order
    mm = paste0(mm, "
  for (i_ in ", ar_order + 1, ":length(", xvar, ")) {
    resid_[i_] = 0")
    for (i in seq_len(ar_order)) {
      mm = paste0(mm, " + \n      ar", i, "_[i_] * resid_sigma_[i_ - ", i, "]")
    }

    # Finish up
    mm = paste0(mm, "
  }")
  }

  return(mm)
}
