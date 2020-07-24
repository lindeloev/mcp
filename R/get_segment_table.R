#' Takes a formula and returns a string representation of y, cp, and rhs
#' @aliases unpack_tildes
#' @keywords internal
#' @inheritParams unpack_rhs
#' @param segment A formula
#' @return A one-row tibble with columns:
#'   * `form`: String. The full formula for this segment.
#'   * `form_y`: String. The expression for y (without tilde)
#'   * `form_cp`: String. The formula for the change point
#'   * `form_rhs`: String. The formula for the right-hand side
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
unpack_tildes = function(segment, i) {
  has_LHS = attributes(stats::terms(segment))$response == 1
  if (!has_LHS & i == 1) {
    stop("No response variable in segment 1.")
  } else if (!has_LHS & i > 1) {
    # If no LHS, add a change point "intercept"
    form_str = paste("1 ", format(segment))
  } else if (has_LHS) {
    # Make regular formula into string
    form_str = format(segment)
  }

  # Check for rel(0)
  if (stringr::str_detect(form_str, "rel\\(0\\)"))
    stop("rel(0) in segment ", i, " is not meaningful. Just write 0 if you want a plateau starting where the earlier segment ended.")

  # List of strings for each section
  chunks = stringr::str_trim(strsplit(form_str, "~")[[1]])


  if (length(chunks) == 2) {
    # Only one tilde. This is the first segment or y is implicit from earlier segment(s)
    return(tibble::tibble(
      form = form_str,
      form_y = ifelse(i == 1, chunks[1], NA),
      form_cp = ifelse(i == 1, NA, paste0(" ~ ", chunks[1])),
      form_rhs = paste0(" ~ ", chunks[2])
    ))
  } else if (length(chunks) == 3) {
    if (i == 1)
      stop("The first segment must have exactly one tilde. Got two.")

    return(tibble::tibble(
      form = form_str,
      form_y = chunks[1],
      form_cp = paste0(" ~ ", chunks[2]),
      form_rhs = paste0(" ~ ", chunks[3])
    ))
  } else {
    stop("Error in segment ", i, ": Got none or more than two tildes in a segment formula.")
  }
}

#' Checks if all terms are in the data
#'
#' @aliases check_terms_in_data
#' @keywords internal
#' @inheritParams unpack_rhs
#' @param form Formula or character (tilde will be prefixed if it isn't already)
#' @param n_terms Int >= 1. Number of expected terms. Will raise error if it doesn't match.
#' @return NULL
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
check_terms_in_data = function(form, data, i, n_terms = NULL) {
  # To formula if it isn't.
  form = to_formula(form)

  # Match varnames to data
  found = all.vars(form) %in% colnames(data)
  if (!all(found))
    stop("Error in segment ", i, ": Term '", paste0(all.vars(form)[!found], collapse="' and '"), "' found in formula but not in data.")

  # Check n_terms
  if (!is.null(n_terms)) {
    check_integer(n_terms, "n_terms", lower = 1)
    if (n_terms != length(found))
      stop("Expected ", n_terms, " terms but got ", length(found), ". Specifically, got: ", paste0(all.vars(form), collapse = ", "))
  }
}


#' Unpacks y variable name
#'
#' @aliases unpack_y
#' @keywords internal
#' @inheritParams unpack_rhs
#' @param form_y Character representation of formula
#' @return A one-row tibble with the columns
#'   * `y`: string. The y variable name.
#'   * `trials`: string. The trials variable name.
#'   * `weights`: string. The weights variable name.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
unpack_y = function(form_y, i, family) {
  # If NA and not segment 1, just return empty
  if (is.na(form_y)) {
    if (i == 1)
      stop("A response must be defined in segment 1, e.g., 'y ~ 1'")

    return(tibble::tibble(
      y = NA,
      trials = NA,
      weights = NA
    ))
  }


  # Split by |
  y_split = strsplit(form_y, "\\|")[[1]]
  if (length(y_split) > 2)
    stop("There can only be zero or one pipe in response. Got '", form_y, "' in segment ", i)

  # RESPONSE
  lhs = y_split[1]
  y_col = attr(stats::terms(to_formula(lhs)), "term.labels")
  if (length(y_col) != 1)
    stop("There should be exactly one response variable. Got ", length(y_col), " in segment ", i)

  # Unpack the stuff on the RHS of the pipe, if it was detected:
  if (length(y_split) == 2) {
    rhs = y_split[2]
    term_labels = attr(stats::terms(to_formula(rhs)), "term.labels")
    ok_terms = stringr::str_detect(term_labels, "trials\\(|weights\\(")
    if (all(ok_terms) == FALSE)
      stop("Only terms `trials()` and `weights()` allowed after the pipe in the response variable. Got '", rhs, "'.")

    # BINOMIAL TRIALS
    trials_term_index = stringr::str_detect(term_labels, "trials\\(")  # Which terms?
    got_trials = sum(trials_term_index) > 0

    # trials(N) is reserved and required for binomial models
    if (family$family == "binomial" & !got_trials)
      stop("Error in response of segment ", i, ": need a valid trials() specification, e.g., 'y | trials(N) ~ 1 + x'")

    if (family$family != "binomial" & got_trials)
      stop("Response format `y | trials(N)` only meaningful for family = binomial(); not for ", family$family, "()")

    # Unpack trials_col
    if (got_trials == TRUE) {
      trials_term = term_labels[trials_term_index]  # Extract terms
      trials_content = get_term_content(trials_term)
      trials_col = attr(stats::terms(trials_content), "term.labels")
      if (length(trials_col) != 1)
        stop("There must be exactly one term inside trials(). Got ", trials_term, " in segment ", i)
    } else {
      trials_col = NA
    }

    # WEIGHTS (same procedure as for trials(N))
    weights_term_index = stringr::str_detect(term_labels, "weights\\(")
    got_weights = sum(weights_term_index) > 0
    if (got_weights == TRUE) {
      if (family$family != "gaussian")
        stop("Weights are currently only implemented for `family = gaussian()`. Raise an issue on GitHub if you need it for other families.")
      weights_term = term_labels[weights_term_index]
      weights_content = get_term_content(weights_term)
      weights_col = attr(stats::terms(weights_content), "term.labels")
      if (length(weights_col) != 1)
        stop("There must be exactly one term inside weights(). Got ", weights_term, " in segment ", i)
    } else {
      weights_col = NA
    }
  } else {
    # No pipe in response formula:
    trials_col = NA
    weights_col = NA
  }

  # Finally:
  return(
    tibble::tibble(
      y = y_col,  # Char
      trials = trials_col,  # Char or NA
      weights = weights_col
    )
  )
}


#' Takes a cp formula (as a string) and returns its properties
#'
#' @aliases unpack_cp
#' @keywords internal
#' @inheritParams unpack_rhs
#' @param form_cp Segment formula as string.
#' @return A one-row tibble with columns:
#'   * `cp_int`: bool. Whether there is an intercept change in the change point.
#'   * `cp_in_rel`: bool. Is this intercept change relative?
#'   * `cp_ran_int`: bool or NA. Is there a random intercept on the change point?
#'   * `cp_group_col`: char or NA. Which column in data define the random intercept?
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_cp = function(form_cp, i) {
  if (is.na(form_cp)) {
    return(tibble::tibble(
      cp_int = FALSE,
      cp_int_rel = FALSE,
      cp_ran_int = FALSE,
      cp_group_col = NA
    ))
  }
  form_cp = stats::as.formula(form_cp)

  # Varying effects
  form_varying = remove_terms(form_cp, "population")

  if (!is.null(form_varying)) {
    varying_terms = attr(stats::terms(form_varying), "term.labels")
    if (length(varying_terms) > 1)
      stop("Error in segment", i, " (change point): only one varying effect allowed. Found ", form_cp)

    varying_parts = strsplit(gsub(" ", "", varying_terms), "\\|")[[1]]
    if (!varying_parts[1] == "1")
      stop("Error in segment ", i, " (change point): Only plain intercepts are allowed in varying effects, e.g., (1|id).", i)

    if (!grepl("^[A-Za-z._0-9]+$", varying_parts[2]))
      stop("Error in segment ", i, " (change point): invalid format of grouping variable in varying effects. Got: ", varying_parts[2])
  }

  # Fixed effects
  attrs = attributes(stats::terms(remove_terms(form_cp, "varying")))
  is_int_rel = attrs$term.labels == "rel(1)"
  if (any(is_int_rel))
    attrs$term.labels = attrs$term.labels[-is_int_rel]  # code as no term

  if (length(attrs$term.labels) > 0)
    stop("Error in segment ", i, " (change point): Only intercepts (1 or rel(1)) are allowed in population-level effects.")

  if (is.null(form_varying) & attrs$intercept == 0)
    stop("Error in segment ", i, " (change point): no intercept or varying effect. You can do e.g., ~ 1 or ~ (1 |id).")

  # Return as list.
  if (!is.null(form_varying)) {
    # If there is a varying effect
    return(tibble::tibble(
      cp_int = attrs$intercept == 1,
      cp_int_rel = any(is_int_rel),  # the intercept is relative
      cp_ran_int = ifelse(varying_parts[1] == "1", TRUE, NA),  # placeholder for later
      cp_group_col = varying_parts[2]
    ))
  } else {
    # If there is no varying effect
    return(tibble::tibble(
      cp_int = attrs$intercept == 1,
      cp_int_rel = any(is_int_rel),
      cp_ran_int = FALSE,
      cp_group_col = NA
    ))
  }
}




#' Unpack right-hand side
#'
#' This is a pretty big function. It includes unpacking sigma(), ar(), etc.
#'
#' @aliases unpack_rhs
#' @keywords internal
#' @param form_rhs A character representation of a formula
#' @param i The segment number (integer)
#' @param family An mcpfamily object returned by `mcp_family()`.
#' @param data A data.frame or tibble
#' @param last_segment The last row in the segment table, made in `get_segment_table()`
#' @return A one-row tibble with three columns for each of `ct`. `sigma`, `ar`, and `ma`:
#'   * `_int`: NA or a one-row tibble describing the intercept.
#'   * `_slope`: NA or a tibble with a row for each slope term.
#'   * `_code`: NA or a char with the JAGS/R code to implement the slope.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
unpack_rhs = function(form_rhs, i, family, data, last_segment) {
  # Get general format
  form_rhs = stats::as.formula(form_rhs)
  attrs = attributes(stats::terms(remove_terms(form_rhs, "varying")))
  term_labels = attrs$term.labels


  #########
  # SIGMA #
  #########
  # Extract sigma terms
  sigma_term_index = stringr::str_detect(term_labels, "sigma\\(")  # Which terms?
  sigma_term = term_labels[sigma_term_index]  # Extract terms
  sigma_form = get_term_content(sigma_term)
  sigma_int = unpack_int(sigma_form, i, "sigma")
  sigma_slope = unpack_slope(sigma_form, i, "sigma", last_segment$sigma_slope[[1]])
  term_labels = term_labels[!sigma_term_index]  # Remove from list of all terms

  # If not specified, sigma_1 is implicit in segment 1.
  if (all(is.na(sigma_int)) == TRUE & all(is.na(sigma_slope)) == TRUE & i == 1 & family$family == "gaussian") {
    sigma_int = tibble::tibble(
      name = "sigma_1",
      rel = FALSE
    )
  }

  #############################
  # CONVERT ARMA TO AR AND MA #
  #############################

  # Split arma into ma and ar for later handling
  arma_term_index = stringr::str_detect(term_labels, "arma\\(")  # Which terms?
  arma_term = term_labels[arma_term_index]  # Extract terms
  term_labels = term_labels[!arma_term_index]  # Remove from list of all terms
  term_labels = c(term_labels,
                  gsub("arma\\(", "ar\\(", arma_term),
                  gsub("arma\\(", "ma\\(", arma_term))


  ######
  # AR #
  ######
  ar_term_index = stringr::str_detect(term_labels, "ar\\(")  # Which terms?
  ar_term = term_labels[ar_term_index]  # Extract terms
  ar_stuff = unpack_arma(ar_term)  # $order and $form_str
  ar_form = get_term_content(ar_stuff$form_str)

  # Populate each of these for each order of AR
  ar_int = list()  # Populate this with intercepts for each order
  ar_slope = list()  # Populate this with slopes for each order
  ar_code = list()
  if (!is.na(ar_stuff$order)) {
    for (order in seq_len(ar_stuff$order)) {
      # Set last_segment = NA if the AR order was lower order in the last segment, thus not defined...
      if (i == 1) {
        this_last = NA
      } else if (order > length(last_segment$ar_slope[[1]])) {
        this_last = NA
      } else {
        this_last = last_segment$ar_slope[[1]][[order]]
      }

      # Get intercept and slope like we're used to
      ar_par_name = paste0("ar", order)  # ar1, ar2, ...
      ar_int[[order]] = unpack_int(ar_form, i, ar_par_name)
      tmp = unpack_slope(ar_form, i, ar_par_name, this_last)
      ar_slope[[order]] = tmp$slope
      ar_code[[order]] = tmp$code
    }
  } else {
    ar_int[[1]] = NA
    ar_slope[[1]] = NA
    ar_code[[1]] = NA
  }
  term_labels = term_labels[!ar_term_index]  # Remove from list of all terms


  ######
  # MA #
  ######
  # Same strategy as for AR
  ma_term_index = stringr::str_detect(term_labels, "ma\\(")  # Which terms?
  ma_term = term_labels[ma_term_index]  # Extract terms
  ma_stuff = unpack_arma(ma_term)  # $order and $form_str
  ma_form = get_term_content(ma_stuff$form_str)

  # Populate each of these for each order of MA
  ma_int = list()  # Populate this with intercepts for each order
  ma_slope = list()  # Populate this with slopes for each order
  ma_code = list()

  if (!is.na(ma_stuff$order)) {
    for (order in seq_len(ma_stuff$order)) {
      # Set last_segment = NA if the MA order was lower order in the last segment, thus not defined...
      if (i == 1) {
        this_last = NA
      } else if (order > length(last_segment$ma_slope[[1]])) {
        this_last = NA
      } else {
        this_last = last_segment$ma_slope[[1]][[order]]
      }

      # Get intercept and slope like we're used to
      ma_par_name = paste0("ma", order)  # ma1, ma2, ...
      ma_int[[order]] = unpack_int(ma_form, i, ma_par_name)
      tmp = unpack_slope(ma_form, i, ma_par_name, this_last)
      ma_slope[[order]] = tmp$slope
      ma_code[[order]] = tmp$code
    }
  } else {
    ma_int[[1]] = NA
    ma_slope[[1]] = NA
    ma_code[[1]] = NA
  }
  term_labels = term_labels[!ma_term_index]  # Remove from list of all terms


  #############################
  # CENTRAL TENDENCIES (MEAN) #
  #############################
  # Start by building it as a string: "ce(1 + x + ...)" to bring it into a compatible format
  if (length(term_labels > 0)) {
    term_labels[1] = paste0(attrs$intercept, " + ", term_labels[1])
    ct_terms = paste0(term_labels, collapse = " + ")  # for use in fit$model and in summary()
    ct_terms = paste0("ct(", ct_terms, ")")  # Get it in "standard" format
  } else {
    ct_terms = paste0("ct(", attrs$intercept, ")")  # Plateau model: "ct(0)" or "ct(1)"
  }
  ct_form = get_term_content(ct_terms)
  ct_int = unpack_int(ct_form, i, "ct")
  ct_slope = unpack_slope(ct_form, i, "ct", last_segment$ct_slope[[1]])


  #################
  # EXTRACT PAR_X #
  #################
  par_x = stats::na.omit(unique(c(
    ifelse(!is.na(sigma_slope$slope), sigma_slope$slope$par_x, NA),
    ifelse(!is.na(ar_slope[[1]]), ar_slope[[1]]$par_x, NA),  # if it's in order 2+, it's also in order1
    ifelse(!is.na(ma_slope[[1]]), ma_slope[[1]]$par_x, NA),  # if it's in order 2+, it's also in order1
    ifelse(!is.na(ct_slope$slope), ct_slope$slope$par_x, NA)
  )))
  if (length(par_x) > 1) {
    stop("Found more than one x-variable in segment ", i, ": '", paste(par_x, collapse ="', '"), "'")
  } else if (length(par_x) == 0) {
    par_x = NA
  }


  ##########
  # RETURN #
  ##########
  return(tibble::tibble(
    x = par_x,

    # Central tendency stuff
    ct_int = list(ct_int),
    ct_code = ct_slope$code,  # must be before ct_slope is defined locally below
    ct_slope = list(ct_slope$slope),
    #ct_varying = list(ct_varying),

    # Sigma stuff
    sigma_int = list(sigma_int),
    sigma_code = sigma_slope$code,  # must be before sigma_slope is defined locally below
    sigma_slope = list(sigma_slope$slope),

    # AR stuff
    ar_int = list(ar_int),
    ar_slope = list(ar_slope),
    ar_code = list(ar_code),

    # MA stuff
    ma_int = list(ma_int),
    ma_slope = list(ma_slope),
    ma_code = list(ma_code)
  ))
}



#' Get formula inside a wrapper
#'
#' @aliases get_term_content
#' @keywords internal
#' @param term E.g., "ct(1 + x)", "sigma(0 + rel(x) + I(x^2))", etc.
#' @return char formula with the content inside the brackets.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_term_content = function(term) {
  # Handle cases of no input or several inputs
  if (length(term) == 0) {
    return(NA)
  } else if (length(term) > 1) {
    #term_type = paste0(substr(term, 0, content_start), ")")
    stop("Only one ", term, " allowed in each formula.")
  } else if (is.na(term)) {
    return(NA)
  } else if (length(term) == 1) {
    # Get formula inside wrapper
    content_start = stringr::str_locate(term, "\\(") + 1  # Location of first character in contents
    content_end = stringr::str_length(term) - 1  # Location of last character in contents
    content = substr(term, content_start, content_end)

    # To formula
    if (content == "")
      stop("Empty terms not allowed in the formulas. Found '", term, "'.")
    form = stats::as.formula(paste0("~", content), env = globalenv())
    return(form)
  }
}


#' Get the intercept of a formula
#'
#' @aliases unpack_int
#' @keywords internal
#' @inheritParams unpack_slope
#' @return A one-row tibble describing the intercept.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_int = function(form, i, ttype) {
  # It there is no formula for this, return NA and don't go further
  if (!rlang::is_formula(form)) {
    return(NA)
  }

  # Got a formula... continue
  attrs = attributes(stats::terms(remove_terms(form, "varying")))
  int_rel = attrs$term.labels == "rel(1)"

  if (any(int_rel) == TRUE & i == 1)
    stop("rel() cannot be used in segment 1. There is nothing to be relative to.")

  if (any(int_rel) | attrs$intercept == TRUE) {
    # Different naming schemes for central tendency (int_i) and others (e.g., sigma_1; ma_1)
    name = ifelse(ttype == "ct", yes = paste0("int_", i), no = paste0(ttype, "_", i))
    int = tibble::tibble(
      name = paste0(name),
      rel = any(int_rel)
    )
  } else {
    int = NA
  }

  return(int)
}


#' Unpack the slope of a formula
#'
#' Makes A list of terms and applies unpack_slope_term() to each of them. Then builds the code for this segment's slope
#' @aliases unpack_slope
#' @keywords internal
#' @param form A formula
#' @param i Segment number (integer)
#' @param ttype The term type. One of "ct" (central tendency), "sigma" (variance),
#'   or "ar" (autoregressive)
#' @param last_slope The element in the slope column for this ttype in the previous
#'   segment. I.e., probably what this function returned last time it was called!
#' @return A "one-row" list with code (char) and a tibble of slopes.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
unpack_slope = function(form, i, ttype, last_slope) {
  # If there is no formula information, return NA
  if (!rlang::is_formula(form)) {
    return(list(
      slope = NA,
      code = NA
    ))
  }

  # Get a list of terms as strings
  attrs = attributes(stats::terms(remove_terms(form, "varying")))
  term_labels = attrs$term.labels[attrs$term.labels != "rel(1)"]

  # Population-level slopes
  if (length(term_labels) > 0) {
    slope = lapply(
      term_labels,
      FUN = unpack_slope_term,
      i = i,
      last_slope = last_slope,
      ttype = ifelse(ttype == "ct", "", paste0(ttype, "_"))  # No prefix for ce ("x_1"). Others are "sigma_x_1"
    ) %>%
      dplyr::bind_rows()  # Make it a proper table-like tibble

    # Build code. Add parentheses if there are more slopes
    code = paste(slope$code, collapse = " + ")
    if (stringr::str_detect(code, "\\+"))
      code = paste0("(", code, ")")

  } else {
    slope = NA
    code = NA
  }

  return(list(
    slope = slope,
    code = code
  ))
}


#' Unpacks a single term
#'
#' Returns a row for `unpack_slope()`.
#'
#' @aliases unpack_slope_term
#' @keywords internal
#' @inheritParams unpack_slope
#' @param term A character, e.g., "x", "I(x^2)", or "log(x)".
#' @return A "one-row" list describing a slope term.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
unpack_slope_term = function(term, i, last_slope, ttype = "") {
  # Remove the "rel()" and "I()" in that order. Must be done before checking starts.
  if (stringr::str_starts(term, "rel\\(")) {
    term = gsub("^rel\\(|\\)$", "", term)
    rel = TRUE
  } else {
    rel = FALSE
  }
  if (stringr::str_starts(term, "I\\("))
    term = gsub("^I\\(|\\)$", "", term)

  # Do some checks on rel()
  if (rel == TRUE) {
    # Check if the relative slope is legit
    if (i == 1)
      stop("rel() cannot be used in segment 1. There is nothing to be relative to.")
    if (all(is.na(last_slope)) == TRUE)
      stop("Found rel(", term, ") in segment ", i, " but there are no corresponding slopes in segment ", i-1)
    if (!term %in% last_slope$term) {
      stop("Found rel(", term, ") in segment ", i, " does not have a corresponding term to be relative to in segment ", i-1)
    }
  }
  if (stringr::str_detect(term, "^log\\(|^sqrt\\(") & i > 1)
    stop("log() or sqrt() detected in segment 2+. This would fail because mcp models earlier segments as negative x values, and sqrt()/log() cannot take negative values.")

  # Regular expressions. Only recognize stuff that is identical between JAGS and base R
  func_list = c("abs", "cos", "exp", "log", "sin", "sqrt", "tan")
  funcs_regex = paste0("^", func_list, "\\(", collapse = "|")  # ^abs(|^cos(|...
  exponent_regex = "\\^\\(?([0-9.-]+)\\)?$"  # something^[number] where [number] can be e.g., "(-1.5)
  end_regex = "\\)$"

  # Find par_x by removing everything associated with accepted functions
  par_x = gsub(paste0(c(funcs_regex, exponent_regex, end_regex), collapse = "|"), "", term)

  # Give the term a valid variable name
  rel_x_code = paste0("X_", i, "_[i_]")  # x relative to segment start. Must match that inserted in get_formula()
  if (stringr::str_detect(term, exponent_regex)) {
    # Exponential
    exponent = sub(".*\\^\\(?([0-9.-]+)\\)?$", "\\1", term)
    check_integer(as.numeric(exponent), term, lower = 0)  # JAGS does not allow for negative or decimal exponents. (TO DO: check if implementing stan version)
    name = paste0(par_x, "_", i, "_E", sub("-", "m", exponent))  # e.g., x_i_E2
    term_recode = gsub(paste0("^", par_x, "\\^"), paste0(rel_x_code, "\\^"), term)
  } else if (stringr::str_detect(term, funcs_regex)) {
    # A simple function
    name = paste0(par_x, "_", i, "_", sub("\\(.*", "\\2", term))  # e.g., x_i_log
    term_recode = gsub(paste0("\\(", par_x), paste0("(", rel_x_code), term)
  } else if (term == par_x) {
    # No function; just vanilla :-)
    name = paste0(par_x, "_", i)
    term_recode = rel_x_code
  } else {
    stop("mcp failed for term ", term, ". If there is no obvious reason why, please raise an issue at GitHub.")
  }

  # Add sigma_name, ma_name, or ar_name
  name = paste0(ttype, name)

  # Add the corresponding term from the last segment if this slope is relative
  if (rel == FALSE) {
    name_cumul = name
    code = paste0(name, " * ", term_recode)
  } else if (rel == TRUE) {
    last_name_cumul = last_slope$name_cumul[which(last_slope$term == term)]
    name_cumul = paste0(last_name_cumul, " + ", name)
    code = paste0("(", name_cumul, ") * ", term_recode)
  }

  return(list(term = term,
              par_x = par_x,
              name = name,
              name_cumul = name_cumul,
              rel = rel,
              code = code))
}



#' Unpack arma order and formula
#'
#' @aliases unpack_arma
#' @keywords internal
#' @param form_str_in A character. These are allowed: "ar(number)" or "ar(number, formula)"
#' @return A list with $order and $form_str (e.g., "ar(formula)"). The formula is ar(1) or ma(1) if no formula is given
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
unpack_arma = function(form_str_in) {
  if (length(form_str_in) == 0) {
    return(list(
      order = NA,
      form_str = NA
    ))
  } else if (length(form_str_in) > 1) {
    stop("Only one of these allowed per segment: ", form_str_in)
  }

  # GET ORDER
  order_start = stringr::str_locate(form_str_in, "\\(") + 1  # Location of first character in contents
  order_end = stringr::str_locate(form_str_in, ",") - 1  # Where is comma? If no comma, this returns NA, NA
  has_formula = !all(is.na(order_end))  # Is there a formula (a comma?)
  if (!has_formula)
    order_end = stringr::str_length(form_str_in) - 1  # No formula; just remove the end parenthesis
  order = suppressWarnings(as.numeric(substr(form_str_in, order_start, order_end)))

  # Check the order
  if (is.na(order))
    stop("Wrong specification of order in '", form_str_in, "'. Must be ar(order) or ar(order, formula) where order is a positive integer.")
  check_integer(order, form_str_in, lower = 1)

  # GET FORMULA
  if (has_formula) {
    # If there is a formula, remove the order
    form_str = gsub(paste0(order, ", "), "", form_str_in)
  } else {
    # If there is no formula, return "ar(1)" or "ma(1)"
    term_type = substr(form_str_in, 1, order_start - 2)  # "ar" or "ma"
    form_str = paste0(term_type, "(1)")
  }

  return(list(
    order = order,
    form_str = form_str
  ))
}



#' Unpack varying effects
#'
#' @aliases unpack_varying_term
#' @keywords internal
#' @inheritParams unpack_slope_term
#' @return A "one-row" list describing a varying intercept.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
unpack_varying_term = function(term, i) {
  parts = stringr::str_trim(strsplit(term, "\\|")[[1]])  # as c("lhs", "group")

  # Check that there is just one grouping term
  if (!grepl("^[A-Za-z._0-9]+$", parts[2]))
    stop("Error in segment ", i, " (linear): Grouping variable in varying effects for change points.")

  # Check that nothing is relative
  if (any(stringr::str_detect(parts, "rel\\(")))
    stop("Error in segment ", i, " (linear): rel() not supported in varying effects.")

  # LHS: Split intercepts and variable
  preds = strsplit(gsub(" ", "", parts[1]), "\\+")[[1]]
  slope = preds[!preds %in% c("0", "1")]
  if (length(slope) > 1)
    stop("Error in segment ", i, " (linear): More than one variable in LHS of varying effect.")
  else if (length(slope) == 0)
    # If not slope
    slope = NA

  # Return
  return(list(
    int = !"0" %in% preds,  # bool. Is intercept present?
    slope = slope,
    group = parts[2]))
}



#' Build a table describing a list of segments
#'
#' Used internally for most mcp functions.
#'
#' @aliases get_segment_table
#' @keywords internal
#' @inheritParams mcp
#' @return A tibble with one row describing each segment and the corresponding code.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @importFrom magrittr %>%
#' @importFrom stats gaussian binomial
#' @export
#' @examples
#' model = list(
#'   y ~ 1 + x,
#'   1 + (1|id) ~ 1
#' )
#' get_segment_table(model)

get_segment_table = function(model, data = NULL, family = gaussian(), par_x = NULL) {
  #####################################################
  # BUILD "SEGMENT TABLE (ST)" FROM ISOLATED SEGMENTS #
  #####################################################
  ST = tibble::tibble()
  for (i in seq_along(model)) {
    # Get ready...
    segment = model[[i]]
    if (!is.null(data))
      check_terms_in_data(segment, data, i)

    # Go! Unpack this segment
    row = tibble::tibble(segment = i)
    row = dplyr::bind_cols(row, unpack_tildes(segment, i))
    row = dplyr::bind_cols(row, unpack_y(row$form_y, i, family))
    row = dplyr::bind_cols(row, unpack_cp(row$form_cp, i))
    row = dplyr::bind_cols(row, unpack_rhs(row$form_rhs, i, family, data, ST[i-1,]))

    ST = dplyr::bind_rows(ST, row)
  }

  # Fill y and trials, where not explicit.
  # Build "full" formula (with explicit intercepts) and insert instead of the old
  ST = ST %>%
    tidyr::fill(.data$y, .data$form_y, .data$trials, .data$weights, .direction="downup") %>%  # Usually only provided in segment 1
    dplyr::mutate(form = paste0(.data$form_y, ifelse(segment == 1, "", .data$form_cp), .data$form_rhs)) %>%  # build full formula
    dplyr::select(-.data$form_y, -.data$form_cp, -.data$form_rhs)  # Not needed anymore




  ###########################
  # CHECK SEGMENTS AND DATA #
  ###########################

  # Check segment 1: rel() not possible here.
  if (any(ST[1, c("cp_int", "cp_int_rel", "cp_ran_int", "cp_group_col")] != FALSE, na.rm = T))
    stop("Change point defined in first segment. This should not be possible. Submit bug report in the GitHub repo.")

  if (nrow(ST) > 1) {
    if (ST$cp_int_rel[2] == TRUE)
      stop("rel() cannot be used for change points in segment 2. There are no earlier change points to be relative to. Relative changepoints work from segment 3 and on.")
  }

  # Set ST$x (what is the x-axis dimension?)
  derived_x = unique(stats::na.omit(ST$x))
  if (length(derived_x) == 1) {
    # One x derived from segments
    if (is.null(par_x))  # par_x not provided. Rely on derived
      ST$x = derived_x
    else if (par_x == derived_x)  # par_x provided and matches devided.
      ST$x = derived_x
    else  # contradicting derived and provided x
      stop("par_x provided but it does not match the predictor found in segment right-hand side")
  } else if (length(derived_x) == 0) {
    # Zero x derived from segments. Rely on par_x?
    if (all(is.na(ST$x) & is.character(par_x)))
      ST$x = par_x
    else
      stop("This is a plateau-only model so no x-axis variable could be derived from the segment formulas. Use argument 'par_x' to set it explicitly")
  } else if (length(derived_x) > 1)
    # More than one...
    stop("More than one predictor found: '", paste0(unique(stats::na.omit(ST$x)), collapse = "' and '"), "'")

  # Response variables
  derived_y = unique(stats::na.omit(ST$y))
  if (length(derived_y) != 1)
    stop("There should be exactly one response variable. Found '", paste0(derived_y, collapse="' and '"), "' across segments.")

  # Weights
  derived_weights = unique(stats::na.omit(ST$weights))
  if (length(derived_weights) > 1)
    stop("There should be exactly zero or one column used for weights(). Found '", paste0(derived_weights, collapse = "' and '"), "' across segments.")

  # Varying effects
  derived_varying = unique(stats::na.omit(ST$cp_group_col))

  # Sigma
  if (any(c(!is.na(ST$sigma_int)), !is.na(ST$sigma_slope)) & family$family != "gaussian")
    stop("sigma() is only meaningful for family = gaussian()")

  # Check data types
  if (!is.null(data)) {
    # Check x and y
    if (!is.numeric(data[, ST$x[1]]))
      stop("Data column '", ST$x[1], "' has to be numeric.")
    if (!is.numeric(data[, ST$y[1]]))
      stop("Data column '", ST$y[1], "' has to be numeric.")
    if (any(is.na(data[, ST$x[1]])))
      stop("NA not allowed in predictor: '", ST$x[1], "'")
    if (any(is.na(data[, ST$y[1]])))
      message("NA values detected in '", ST$y[1], "'. They will be imputed using the posterior predictive.")

    # Check varying
    if (length(derived_varying) > 0) {
      for (varying_col in derived_varying) {
        data_varying = data[, varying_col]
        if (!is.character(data_varying) & !is.factor(data_varying))
          if (!all(data_varying == floor(data_varying)))
            stop("Varying group '", varying_col, "' has to be integer, character, or factor.")
      }
    }

    # Check y and trials if binomial
    if (family$family == "binomial") {
      check_integer(data[, ST$y[1]], ST$y[1], lower = 0)
      check_integer(data[, ST$trials[1]], ST$trials[1], lower = 1)
    } else if (family$family == "bernoulli") {
      if (any(!data[, ST$y[1]] %in% c(0, 1)))
        stop("Only responses 0 and 1 are allowed for family = bernoulli() in column '", ST$y[1], "'")
    } else if (family$family == "poisson") {
      check_integer(data[, ST$y[1]], ST$y[1], lower = 0)
    }

    # Check weights
    if (length(derived_weights) == 1) {
      if (!is.numeric(data[, derived_weights]))
        stop("Data column '", derived_weights, "' has to be numeric.")
      if (any(data[, derived_weights] <= 0))
        stop("All weights must be greater than zero.")
    }
  }


  ###################################
  # MUTATE ST WITH NEW HANDY VALUES #
  ###################################

  # Recode relative columns so 0 = not relative. N > 0 is the number of consecutive "relatives"
  ST = ST %>%

    # Add variable names
    dplyr::mutate(
      cp_name = paste0("cp_", .data$segment - 1),
      cp_sd = ifelse(.data$cp_ran_int == TRUE, paste0(.data$cp_name, "_sd"), NA),
      cp_group = ifelse(.data$cp_ran_int == TRUE, paste0(.data$cp_name, "_", .data$cp_group_col), NA)
    ) %>%

    # Add "cumulative" cp_code_form for concecutive relative intercepts
    dplyr::group_by(cumsum(!.data$cp_int_rel)) %>%
    dplyr::mutate(
      cp_code_prior = cumpaste(.data$cp_name, " + "),
      cp_code_form = ifelse(!is.na(.data$cp_group), yes = paste0(.data$cp_code_prior, " + ", .data$cp_group, "CP_", .data$segment, "_INDEX"), no = .data$cp_code_prior),
      cp_code_form = format_code(.data$cp_code_form, na_col = .data$cp_name),
      cp_code_prior = format_code(.data$cp_code_prior, na_col = .data$cp_name)
    ) %>%
    dplyr::ungroup() %>%

    # Finish up
    dplyr::select(-dplyr::starts_with("cumsum"))

  # Return
  ST
}


#' Cumulative pasting of character columns
#'
#' @aliases cumpaste
#' @keywords internal
#' @param x A column
#' @param .sep A character to append between pastes
#' @return string.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk} but Inspired by
#'   https://stackoverflow.com/questions/24862046/cumulatively-paste-concatenate-values-grouped-by-another-variable
#'
cumpaste = function(x, .sep = " ")
  Reduce(function(x1, x2) paste(x1, x2, sep = .sep), x, accumulate = TRUE)


#' Format code with one or multiple terms
#'
#' Take a value like "a + b" and
#' (1) replace it with NA if na_col == NA.
#' (2) Change to "(a + b)" if there is a "+"
#' (3) Return itself otherwise, e.g., "a" --> "a".
#'
#' @aliases format_code
#' @keywords internal
#' @param col A column
#' @param na_col If this column is NA, return NA
#' @return string
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
format_code = function(col, na_col) {
  # Replace with NA
  col = ifelse(is.na(na_col), NA, col)

  # Add parenthesis if this is an equation
  col = ifelse(stringr::str_detect(col, "\\+"), paste0("(", col, ")"), col)
  col
}
