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
        formula_str = paste0(formula_str, "# Fitted standard deviation\nsigma_[i_] = \n")
      } else if (stringr::str_detect(ytype, "ar[0-9]+")) {
        formula_str = paste0(formula_str, "# Autoregressive: AR(", ar_order,")\n", ytype, "_[i_] = \n")
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
#'
get_simulate = function(formula_str, par_x, par_trials = NA, pars_pop, pars_varying, nsegments, family) {
  # First some substitutions
  formula_func = gsub("\\[i_\\]", "", formula_str)  # No explicit indexing needed for R function
  #formula_func = gsub("PAR_X", par_x, formula_str)  # Proper par_x name
  formula_func = gsub("min\\(", "pmin\\(", formula_func)  # vectorized mean for function
  for (i in seq_len(nsegments)) {
    formula_func = gsub(paste0("CP_", i, "_INDEX"), "", formula_func)  # Use vector of data
  }

  # Now build the function R code
  #x_and_trials handles the special case that binomial "trials" is also used as a predictor.
  x_and_trials = ifelse(is.na(par_trials), par_x, no = paste0(unique(c(par_x, par_trials)), collapse = ", "))
  func_str = paste0("
  function(",
    x_and_trials, ", ",
    paste0(pars_pop, collapse = ", "), ", ",
    paste0(pars_varying, collapse = " = 0, "), ifelse(length(pars_varying) > 0, " = 0, ", ""),
    "type = 'predict', quantiles = FALSE, rate = FALSE, ...) {
    # Helpers to simplify making the code for this function
    cp_0 = -Inf
    cp_", nsegments, " = Inf

    ", formula_func, "


    # Return predictions or fitted values?
    if (!type %in% c('predict', 'fitted'))
      stop(\"'`type` must be one of 'predict' or 'fitted'\")
    ")

  # Return depends on family
  if (family$family == "gaussian") {
    func_str = paste0(func_str, "
    if (type == 'fitted') return(y_)
    if (type == 'predict') {
      if (quantiles == TRUE)
        quantiles = c(0.025, 0.975)
      if (is.numeric(quantiles)) {
        return(qnorm(quantiles, y_, sigma_))
      } else if (quantiles == FALSE) {
        return(rnorm(length(", par_x, "), y_, sigma_))
      } else {
        stop('Invalid quantiles argument to simulate()')
      }
    }")
  } else if (family$family == "binomial") {
    func_str = paste0(func_str, "
    if (type == 'predict') {
      if (rate == FALSE) return(rbinom(length(", par_x, "), ", par_trials, ", ilogit(y_)))
      if (rate == TRUE)  return(rbinom(length(", par_x, "), ", par_trials, ", ilogit(y_)) / ", par_trials, ")
    }
    if (type == 'fitted')
      if (rate == FALSE) return(", par_trials, " * ilogit(y_))
      if (rate == TRUE)  return(ilogit(y_))")
  } else if (family$family == "bernoulli") {
    func_str = paste0(func_str, "
    if (type == 'predict') return(rbinom(length(", par_x, "), 1, ilogit(y_)))
    if (type == 'fitted') return(ilogit(y_))")
  } else if (family$family == "poisson") {
    func_str = paste0(func_str, "
    if (type == 'predict') return(rpois(length(", par_x, "), exp(y_)))
    if (type == 'fitted') return(exp(y_))")
  }

  func_str = paste0(func_str, "
  }")

  # Return function
  eval(parse(text = func_str))
}
