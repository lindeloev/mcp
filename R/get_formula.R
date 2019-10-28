#' Build an R formula (as string) given a segment table (ST)
#'
#' You will need to replace PAR_X for whatever your x-axis observation column
#' is called. In JAGS typicaly `x[i_]`. In R just `x`.
#'
#' @aliases get_formula_str
#' @param ST Tibble. Returned by \code{get_segment_table}.
get_formula_str = function(ST) {

  formula_str = ""
  for (i in seq_len(nrow(ST))) {
    # Helper: Current segment.
    S = ST[i, ]

    # PAR_X will be replaced by the real name. FUTURE_REL will sometimes be
    # replaced by a less-than indicator (ind_past).
    # cp_code_form includes varying effects
    ind_this = paste0("(PAR_X >= ", S$cp_code_form, ") * FUTURE_REL")
    ind_past = paste0("(PAR_X < ", S$cp_code_form, ") * ")

    # Begin building formula
    formula_str = paste0(formula_str, "\n# Segment ", i, ": ", S$form, "\n")


    # Add intercept
    if (S$int == TRUE) {
      # For absolute intercepts, "remove" earlier stuff affecting the intercept
      # Multiply it with zero from this change point and on
      if (S$cp_int_rel == FALSE) {
        formula_str = gsub("FUTURE_REL", ind_past, formula_str)
      }

      # Add intercept with indicator
      formula_str = paste0(formula_str, ind_this, S$int_name, " + \n")
    }

    # Add slope
    if (!is.na(S$slope)) {
      # What is the start of this segment? Used to subtract so that slope is computed form segment start.
      subract_cp = ifelse(i > 1, yes = paste0(" - (", S$cp_code_form, ")"), no = "")
      next_cp = ifelse(i < nrow(ST), yes = paste0(ST$cp_code_form[i + 1], ""), no = paste0("cp_", i))  # infinite if last segment.

      formula_str = paste0(formula_str, ind_this, "(", S$slope_code, " * (min(PAR_X, ", next_cp, ")", subract_cp, ")) + \n")
    }

    # If this is just a plateau (~ 0), i.e., the absence of intercept, slope, and varying effects
    if (!S$int & is.na(S$slope) & !tibble::is_tibble(S$varying[[1]])) {
      formula_str = paste0(formula_str, "0 + \n")
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
#' @aliases get_func_y
#' @inheritParams mcp
#' @param formula_str string. Returned by \code{get_formula}.
#' @param par_x string. Same as for \link{mcp}.
#' @param par_trials String. For binomial models: name of trials column.
#' @param pars_pop List of population parameters which the user must provide.
#' @param pars_varying List varying parameters. They will default to zero
#'   (optional for the user).
#' @param nsegments Positive integer. Number of segments, typically \code{nrow(ST)}.
#'
get_func_y = function(formula_str, par_x, par_trials = NA, pars_pop, pars_varying, nsegments, family) {
  # First some substitutions
  formula_func = gsub("PAR_X", par_x, formula_str)  # Proper par_x name
  formula_func = gsub("min\\(", "pmin\\(", formula_func)  # vectorized mean for function
  for (i in seq_len(nsegments)) {
    formula_func = gsub(paste0("CP_", i, "_INDEX"), "", formula_func)  # Use vector of data
  }

  # Now build the function R code
  func_str = paste0("
  function(",
    par_x, ", ",
    ifelse(!is.na(par_trials), yes = paste0(par_trials, ", "), no = ""),
    paste0(pars_pop, collapse = ", "), ", ",
    paste0(pars_varying, collapse = " = 0, "), ifelse(length(pars_varying) > 0, " = 0, ", ""),
    "type = 'predict', rate = FALSE, ...) {
    # Helpers to simplify making the code for this function
    cp_0 = -Inf
    cp_", nsegments, " = Inf

    # Fitted value
    y = ", formula_func, "


    # Return predictions or fitted values?
    if (!type %in% c('predict', 'fitted'))
      stop(\"'`type` must be one of 'predict' or 'fitted'\")
    ")

  # Return depends on family
  if (family == "gaussian") {
    func_str = paste0(func_str, "
    if (type == 'predict') return(rnorm(length(", par_x, "), y, sigma))
    if (type == 'fitted') return(y)")
  } else if (family == "binomial") {
    func_str = paste0(func_str, "
    inverse_logit = function(x) exp(x) / (1 + exp(x))
    if (type == 'predict') {
      if (rate == FALSE) return(rbinom(length(", par_x, "), ", par_trials, ", inverse_logit(y)))
      if (rate == TRUE)  return(rbinom(length(", par_x, "), ", par_trials, ", inverse_logit(y)) / ", par_trials, ")
    }
    if (type == 'fitted')
      if (rate == FALSE) return(", par_trials, " * inverse_logit(y))
      if (rate == TRUE)  return(inverse_logit(y))")
  }

  func_str = paste0(func_str, "
  }")

  # Return function
  eval(parse(text = func_str))
}
