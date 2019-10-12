#' Understand the list of segments
#'
#' @aliases unpack_segments
#' @inheritParams mcp
#' @importFrom stats terms
#' @importFrom utils modifyList
#' @import dplyr
#' @return A list with prior (list), formula_jags (string), func_y (function),
#'   param_x (string), param_y(string).
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}


unpack_segments = function(segments, prior = list(), param_x = NULL) {
  ##############
  # PREPROCESS #
  ##############
  all_terms = unlist(lapply(segments, function(x) attr(terms(x), "term.labels")))

  # Try extracting param_x from segments.
  if(is.null(param_x)) {
    param_x = all_terms[!all_terms %in% c("rel(0)", "rel(1)")]
    param_x = gsub("rel\\(|\\)", "", param_x)
    param_x = unique(param_x)

    # Convert from character(0) to NULL if param_x was not found in formulas...
    if(identical(param_x, character(0))) {
      param_x = NULL
    }
  }

  # Detect param_y
  param_y = as.character(segments[[1]][[2]])



  ##########
  # CHECKS #
  ##########
  # Segments
  if(length(segments) == 0) {
    stop("At least one segment is needed")
  }
  if(any(stringr::str_detect(as.character(segments[[1]]), "rel\\("))) {
    stop("`rel` cannot be used in first segment. Nothing to be relative to.")
  }

  # Predictor(s)
  if(is.null(param_x)) {
    stop("This is a plateau-only model, but `param_x` is missing in `mcp`")
  }
  if("rel(0)" %in% all_terms) {
    stop("rel(0) is not supported in segment formulas.")
  }
  if(length(param_x) > 1) {
    stop(paste0("More than one slope variable in segment formulas: ", paste0(param_x, collapse=",")))
  }

  # Response
  test_form = function(x) length(as.character(x)) == 3  # TRUE: 1 ~ 1, FALSE: ~ 1
  if(!all(unlist(lapply(segments, test_form)))) {
    stop("One or more segments miss a response variable.")
  }


  ###############
  # BUILD MODEL #
  ###############
  default_prior_slope = "dt(0, SDY / (MAXX - MINX) * 3, 1)"
  default_prior_int = "dt(0, SDY * 3, 1)"


  # Build list of priors and a formula
  default_prior = list(sigma = "dnorm(0, SDY) T(0, )")
  formula_str = ""

  for(i in 1:length(segments)) {

    # RHS as list of chars. Used to detect intercepts and slopes
    RHS = stringr::str_trim(unlist(strsplit(as.character(segments[[i]])[[3]], "\\+")))
    LHS = stringr::str_trim(unlist(strsplit(as.character(segments[[i]])[[2]], "\\+")))


    # Settings
    has_intercept = any(c("1", "rel(1)") %in% RHS)
    rel_param_x = paste0("rel(", param_x, ")")
    has_slope = any(c(param_x, rel_param_x) %in% RHS)
    formula_str = paste0(formula_str, "\n# Segment ", i, ": ~ ", paste0(RHS, collapse = " + "), "\n")

    # Prepare change point
    if(i > 1) {
      if(as.character(segments[[i]])[2] != "1") {
        stop("Change point can currently only be 1.")
      }

      cp_name = paste0("cp_", i-1)
      min_val = ifelse(i > 2, paste0("cp_", i-2), "MINX")
      default_prior[[cp_name]] = paste0("dunif(", min_val, ", MAXX)")

      # If there is a user-supplied prior that is not dunif nor truncated, add truncation
      if(cp_name %in% names(prior)) {
        if(!stringr::str_detect(prior[[cp_name]], "dunif|T\\(")) {
          prior[[cp_name]] = paste0(prior[[cp_name]], " T(cp_", i-2, ", )")
        }
      }
    }

    # Prepare indicators. ind_abs = absolute, ind_rel = relative
    # PARAM_X will be replaced by the real name. FUTURE_REL will sometimes be
    # replaced by a less-than indicator (ind_past).

    ind_this = paste0("(PARAM_X >= cp_", i - 1 , ") * FUTURE_REL")
    ind_past = paste0("(PARAM_X < cp_", i - 1, ") * ")

    # Add intercept
    if(has_intercept) {
      int_name = paste0("int_", i)
      default_prior[[int_name]] = default_prior_int

      # "remove" earlier terms for absolute intemin(rcepts
      if("1" %in% RHS) {
        formula_str = gsub("FUTURE_REL", ind_past, formula_str)
      }

      formula_str = paste0(formula_str, ind_this, int_name, " + \n")
    }

    # Add slope
    if(has_slope) {
      slope_name = paste0(param_x, "_", i)
      default_prior[[slope_name]] = default_prior_slope

      # If "slope", slope_str = this_slope. If "rel(slope)", slope_str = (this_slope - last_slope)
      slope_str = slope_name  # candidate until proven wrong
      if(rel_param_x %in% RHS) {  # If rel(slope) in formula
        past_slope_name = paste0(param_x, "_", i-1)
        if(past_slope_name %in% names(default_prior)) {  # if the last segment had a slope
          slope_str = paste0("(", slope_name, " + ", param_x, "_", i-1, ")")
        }
      }

      #formula_str = paste0(formula_str, indicator, "((min(", slope_name, ", cp_", i, ") - cp_", i-1, ") + \n")
      subract_cp = ifelse(i == 1, "",  paste0("- cp_", i-1))
      formula_str = paste0(formula_str, ind_this, "(", slope_str, " * (min(PARAM_X, cp_", i, ")", subract_cp, ")) + \n")
    }

    if(!has_intercept & !has_slope) {
      formula_str = paste0(formula_str, "0 + \n")
    }

    # Finish up formula_str
    if(i == length(segments)) {
      formula_str = substr(formula_str, 1, nchar(formula_str) - 3)
    }
  }

  # Overwrite default with user's prior. Add truncation to cp_ if not already
  prior = modifyList(default_prior, prior)

  # Remove all "unrealized" FUTURE_REL
  formula_str = gsub("FUTURE_REL", "", formula_str)

  # Substitute to get a formula to insert in JAGS
  formula_jags = gsub("PARAM_X", paste0(param_x, "[i_]"), formula_str)

  # Turn formula_str into a proper R function
  formula_func = gsub("PARAM_X", param_x, formula_str)  # Proper param_x name
  formula_func = gsub("min\\(", "pmin\\(", formula_func)  # vectorized mean for function

  func_str = paste0("
  function(", param_x, ", ", paste0(names(default_prior), collapse=", "), ", type=\"predict\") {
    # Helpers to simplify making the code for this function
    cp_0 = -Inf
    cp_", length(segments), " = Inf

    # Fitted value
    y = ", formula_func, "

    # Return predictions or fitted values?
    if(type == 'predict') return(rnorm(length(", param_x, "), y, sigma))
    if(type == 'fitted') return(y)
    if(!type %in% c('predict', 'fitted')) stop(\"`type` must be one of 'predict' or 'fitted'\")
  }")

  func_y = eval(parse(text = func_str))


  return(list(
    prior = modifyList(default_prior, prior),
    formula_jags = formula_jags,
    func_y = func_y,
    param_x = param_x,
    param_y = param_y
  ))
}
