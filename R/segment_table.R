#' Takes a formula and returns a string representation of y, cp, and rhs
#' @aliases unpack_tildes
#' @keywords internal
#' @param form A formula
#' @param i The segment number
#' @return A one-row tibble with columns:
#'   * `form`: String. The full formula for this segment.
#'   * `form_y`: String. The expression for y (without tilde)
#'   * `form_cp`: String. The formula for the change point.
#'   * `form_rhs`: String. The formula for the RHS. Only used to build the formula representation in `summary.mcpfit()`.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_tildes = function(form, i) {
  has_LHS = attributes(stats::terms(form))$response == 1
  form_str = formula_to_char(form)
  if (has_LHS == FALSE && i == 1) {
    stop("No response variable in segment 1.")
  } else if (has_LHS == FALSE && i > 1) {
    # If no LHS, add a change point "intercept"
    form_str = paste("1", form_str)
  }

  # List of strings for each section
  chunks = stringr::str_trim(strsplit(form_str, "~")[[1]])

  if (length(chunks) == 2) {
    # Only one tilde. This is the first segment or y is implicit from earlier segment(s)
    return(tibble::tibble(
      form = form_str,
      form_y = ifelse(i == 1, chunks[1], NA),
      form_cp = ifelse(i == 1, NA, paste0(" ~ ", chunks[1])),
      form_rhs = chunks[2]
    ))
  } else if (length(chunks) == 3) {
    if (i == 1)
      stop("The first segment must have exactly one tilde. Got two.")

    return(tibble::tibble(
      form = form_str,
      form_y = chunks[1],
      form_cp = paste0(" ~ ", chunks[2]),
      form_rhs = chunks[3]
    ))
  } else {
    stop("Error in segment ", i, ": Got none or more than two ~ in a segment formula.")
  }
}


#' Unpacks y variable name
#'
#' @aliases unpack_y
#' @keywords internal
#' @inheritParams mcp
#' @param form_y Character representation of formula
#' @param i Segment number
#' @return A one-row tibble with the columns
#'   * `y`: string. The y variable name.
#'   * `trials`: string. The trials variable name.
#'   * `weights`: string. The weights variable name.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
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
    if (family$family == "binomial" && !got_trials)
      stop("Error in response of segment ", i, ": need a valid trials() specification, e.g., 'y | trials(N) ~ 1 + x'")

    if (family$family != "binomial" && got_trials)
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
  tibble::tibble(
    y = y_col,  # Char
    trials = trials_col,  # Char or NA
    weights = weights_col
  )
}


#' Takes a cp formula (as a string) and returns its properties
#'
#' @aliases unpack_cp
#' @keywords internal
#' @param form_cp Segment formula as string.
#' @param i segment number
#' @return A one-row tibble with columns:
#'   * `cp_int`: bool. Whether there is an intercept change in the change point.
#'   * `cp_ran_int`: bool or NA. Is there a random intercept on the change point?
#'   * `cp_group_col`: char or NA. Which column in data define the random intercept?
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_cp = function(form_cp, i) {
  if (is.na(form_cp)) {
    return(tibble::tibble(
      cp_int = FALSE,
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
  if (length(attrs$term.labels) > 0)
    stop("Error in segment ", i, " (change point): Only intercepts (1) are allowed in population-level effects.")

  if (is.null(form_varying) && attrs$intercept == 0)
    stop("Error in segment ", i, " (change point): no intercept or varying effect. You can do e.g., ~ 1 or ~ (1 |id).")

  # Return as list.
  if (!is.null(form_varying)) {
    # If there is a varying effect
    return(tibble::tibble(
      cp_int = attrs$intercept == 1,
      cp_ran_int = ifelse(varying_parts[1] == "1", TRUE, NA),  # placeholder for later
      cp_group_col = varying_parts[2]
    ))
  } else {
    # If there is no varying effect
    return(tibble::tibble(
      cp_int = attrs$intercept == 1,
      cp_ran_int = FALSE,
      cp_group_col = NA
    ))
  }
}



# #' Unpack varying effects
# #'
# #' @aliases unpack_varying_term
# #' @keywords internal
# #' @param term
# #' @return A "one-row" list describing a varying intercept.
# #' @encoding UTF-8
# #' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
# unpack_varying_term = function(term, i) {
#   parts = stringr::str_trim(strsplit(term, "\\|")[[1]])  # as c("lhs", "group")
#
#   # Check that there is just one grouping term
#   if (!grepl("^[A-Za-z._0-9]+$", parts[2]))
#     stop("Error in segment ", i, " (right-hand side): Grouping variable in varying effects for change points.")
#
#   # LHS: Split intercepts and variable
#   preds = strsplit(gsub(" ", "", parts[1]), "\\+")[[1]]
#   slope = preds[(preds %in% c("0", "1")) == FALSE]
#   if (length(slope) > 1)
#     stop("Error in segment ", i, " (right-hand side): More than one variable in LHS of varying effect.")
#   else if (length(slope) == 0)
#     # If not slope
#     slope = NA
#
#   # Return
#   list(
#     int = "0" %notin% preds,  # bool. Is intercept present?
#     slope = slope,
#     group = parts[2]
#   )
# }



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
#' @export
#' @examples
#' model = list(
#'   y ~ 1 + x,
#'   1 + (1|id) ~ 1
#' )
#' get_segment_table(model, par_x = "x")

get_segment_table = function(model, data = NULL, family = gaussian(), par_x) {
  assert_types(par_x, "character", len = 1)

  #####################################################
  # BUILD "SEGMENT TABLE (ST)" FROM ISOLATED SEGMENTS #
  #####################################################
  ST = tibble::tibble()
  for (i in seq_along(model)) {
    row = tibble::tibble(segment = i)
    row = dplyr::bind_cols(row, unpack_tildes(model[[i]], i))
    row = dplyr::bind_cols(row, unpack_y(row$form_y, i, family))
    row = dplyr::bind_cols(row, unpack_cp(row$form_cp, i))

    ST = dplyr::bind_rows(ST, row)
  }

  # Fill y and trials, where not explicit.
  # Build "full" formula (with explicit intercepts) and insert instead of the old
  ST = ST %>%
    tidyr::fill(.data$y, .data$form_y, .data$trials, .data$weights, .direction = "downup") %>%  # Usually only provided in segment 1
    dplyr::mutate(form = ifelse(.data$segment == 1, .data$form, paste0(.data$form_y, .data$form_cp, " ~ ", .data$form_rhs))) %>%  # build full formula
    dplyr::select(-.data$form_y, -.data$form_cp, -.data$form_rhs)  # Not needed anymore

  ST$x = par_x  # TO DO: TEMPORARY UNTIL MULTIPLE REGRESSION IS IMPLEMENTED



  ###########################
  # CHECK SEGMENTS AND DATA #
  ###########################

  # Check segment 1: change point not possible here
  if (any(ST[1, c("cp_int", "cp_ran_int", "cp_group_col")] != FALSE, na.rm = T))
    stop("Change point defined in first segment. This should not be possible. Submit bug report in the GitHub repo.")

  # Response variables
  derived_y = unique(stats::na.omit(ST$y))
  if (length(derived_y) != 1)
    stop("There should be exactly one response variable. Found ", and_collapse(derived_y), " across segments.")

  # Weights
  derived_weights = unique(stats::na.omit(ST$weights))
  if (length(derived_weights) > 1)
    stop("There should be exactly zero or one column used for weights(). Found ", and_collapse(derived_weights), " across segments.")

  # Varying effects
  derived_varying = unique(stats::na.omit(ST$cp_group_col))

  # Check data types
  if (!is.null(data)) {
    # Check y
    if (!is.numeric(data[, ST$y[1]]))
      stop("Data column '", ST$y[1], "' has to be numeric.")
    if (any(is.na(data[, ST$y[1]])))
      message("NA values detected in '", ST$y[1], "'. They will be imputed using the posterior predictive.")

    # Check varying
    if (length(derived_varying) > 0) {
      for (varying_col in derived_varying) {
        data_varying = data[, varying_col]
        if (!is.character(data_varying) && !is.factor(data_varying))
          if (!all(data_varying == floor(data_varying)))
            stop("Varying group '", varying_col, "' has to be integer, character, or factor.")
      }
    }

    # Check y and trials if binomial
    if (family$family == "binomial") {
      assert_integer(data[, ST$y[1]], ST$y[1], lower = 0)
      assert_integer(data[, ST$trials[1]], ST$trials[1], lower = 1)
    } else if (family$family == "bernoulli") {
      if (any(data[, ST$y[1]] %notin% c(0, 1)))
        stop("Only responses 0 and 1 are allowed for family = bernoulli() in column '", ST$y[1], "'")
    } else if (family$family == "poisson") {
      assert_integer(data[, ST$y[1]], ST$y[1], lower = 0)
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
  ST %>%

    # Add variable names
    dplyr::mutate(
      cp_name = paste0("cp_", .data$segment - 1),
      cp_sd = ifelse(.data$cp_ran_int == TRUE, paste0(.data$cp_name, "_sd"), NA),
      cp_group = ifelse(.data$cp_ran_int == TRUE, paste0(.data$cp_name, "_", .data$cp_group_col), NA)
    ) %>%

    # Code varying change points
    dplyr::mutate(
      cp_code_form = ifelse(!is.na(.data$cp_group), yes = paste0(.data$cp_name, " + ", .data$cp_group, "CP_", .data$segment, "_INDEX"), no = .data$cp_name),
      cp_code_form = format_code(.data$cp_code_form, na_col = .data$cp_name)
    )
}


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
format_code = function(col, na_col) {
  # Replace with NA
  col = ifelse(is.na(na_col), NA, col)

  # Add parenthesis if this is an equation
  col = ifelse(stringr::str_detect(col, "\\+"), paste0("(", col, ")"), col)
  col
}
