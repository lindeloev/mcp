# TO DO:
# * Nest cp varying effects too?

source("R/lme4_utils.R")


# Takes a formula and returns a string representation of y, cp, and rhs
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
  if ("rel(0)" %in% attributes(stats::terms(segment))$term.labels)
    stop("rel(0) is not (currently) supported in segment formulas.")

  # List of strings for each section
  chunks = stringr::str_trim(strsplit(form_str, "~")[[1]])

  if (length(chunks) == 2) {
    # Only one tilde. This is the first segment or y is implicit from earlier segment(s)
    if (i == 1 & stringr::str_detect(chunks[1], "^[-0-9.]+"))  # check illigal terms
      stop("y must be a variable.")
    return(tibble::tibble(
      form = form_str,
      form_y = ifelse(i == 1, chunks[1], NA),
      form_cp = ifelse(i == 1, NA, paste0(" ~ ", chunks[1])),
      form_rhs = paste0(" ~ ", chunks[2])
    ))
  } else if (length(chunks) == 3) {
    return(tibble::tibble(
      form = form_str,
      form_y = chunks[1],
      form_cp = paste0(" ~ ", chunks[2]),
      form_rhs = paste0(" ~ ", chunks[3])
    ))
  } else {
    stop("Error in segment ", i, ": None or more than two tildes in a segment formula.")
  }
}

# Checks if all terms are in the data
check_terms_in_data = function(form, data, i) {
  form_terms = all.vars(form)
  found = form_terms %in% colnames(data)
  if (!all(found))
    stop("Error in segment ", i, ": Term '", paste0(form_terms[!found], collapse="' and '"), "' found in formula but not in data.")
}


# Unpacks y variable name
unpack_y = function(form_y, i, family) {
  # Segment 1 must have y specified
  if (i == 1) {
    if (is.na(form_y))
      stop("No response defined in segment 1.")

    if (family == "binomial" & !stringr::str_detect(gsub(" ", "", form_y), "\\|trials\\("))
      stop("Error in response of segment 1: need a valid trials() specification, e.g., 'y | trials(N) ~ 1 + x'")
  }

  if (!is.na(form_y)) {
    form_y = gsub(" ", "", form_y)  # remove whitespace

    # If binomial, get trials
    chunks = gsub(")", "", strsplit(gsub(" ", "", form_y), "\\|trials\\(")[[1]])

    if (!all(grepl("^[A-Za-z._0-9]+$", chunks)))
      stop("Error in segment ", i, ": Invalid format for response variable. Only single column names are allowed.")

    if (length(chunks) == 1) {
      return(tibble::tibble(
        y = chunks[1],
        trials = NA))
    } else if (length(chunks) == 2) {
      return(tibble::tibble(
        y = chunks[1],
        trials = chunks[2]))
    } else {
      stop("Invalid value '", form_y, "' for response. More than one pipe?")
    }
  } else {
    return(tibble::tibble(
      y = NA,
      trials = NA
    ))
  }
}


#' Takes a cp formula (as astring) and returns its properties
#' @param form_cp Segment formula as string.
#' @param i Positive integer. Segment number.
unpack_cp = function(form_cp, i) {
  if (is.na(form_cp)) {
    return(tibble::tibble(
      cp_int = FALSE,
      cp_int_rel = FALSE,
      cp_ran_int = FALSE,
      cp_group_col = NA
    ))
  }
  form_cp = as.formula(form_cp)

  # Varying effects
  varying = findbars(form_cp)
  if (!is.null(varying)) {
    if (length(varying) > 1)
      stop("Error in segment", i, " (change point): only one varying effect allowed. Found ", form_cp)

    varying_parts = as.character(varying[[1]])
    if (!varying_parts[2] %in% "1")
      stop("Error in segment ", i, " (change point): Only plain intercepts are allowed in varying effects, e.g., (1|id).", i)

    if (!grepl("^[A-Za-z._0-9]+$", varying_parts[2]))
      stop("Error in segment ", i, " (change point): Grouping variable in varying effects.")
  }

  # Fixed effects
  population = attributes(stats::terms(nobars(form_cp)))
  is_int_rel = population$term.labels == "rel(1)"
  if (any(is_int_rel))
    population$term.labels = population$term.labels[-is_int_rel]  # code as no term

  if (length(population$term.labels) > 0)
    stop("Error in segment ", i, " (change point): Only intercepts (1 or rel(1)) are allowed in population-level effects.")

  if (is.null(varying) & population$intercept == 0)
    stop("Error in segment ", i, " (change point): no intercept or varying effect. You can do e.g., ~ 1 or ~ (1 |id).")

  # Return as list.
  if (!is.null(varying)) {
    # If there is a varying effect
    return(tibble::tibble(
      cp_int = population$intercept == 1,
      cp_int_rel = any(is_int_rel),  # the intercept is relative
      cp_ran_int = ifelse(varying_parts[2] == "1", TRUE, NA),  # placeholder for later
      cp_group_col = varying_parts[3]
    ))
  } else {
    # If there is no varying effect
    return(tibble::tibble(
      cp_int = population$intercept == 1,
      cp_int_rel = any(is_int_rel),
      cp_ran_int = FALSE,
      cp_group_col = NA
    ))
  }
}

unpack_rhs = function(form_rhs, i) {
  form_rhs = formula(form_rhs)

  # Varying effects
  varying = findbars(form_rhs)
  if (!is.null(varying)) {
    # For each varying effect...
    V = tibble::tibble()
    for (term in varying) {
      parts = as.character(term)

      # Check that there is just one grouping term
      if (!grepl("^[A-Za-z._0-9]+$", parts[3]))
        stop("Error in segment ", i, " (linear): Grouping variable in varying effects for change points.")

      # Check that nothing is relative
      if (any(stringr::str_detect(parts, "rel\\(")))
        stop("Error in segment ", i, " (linear): rel() not supported in varying effects.")

      # LHS: Split intercepts and variable
      preds = strsplit(gsub(" ", "", parts[2]), "\\+")[[1]]
      slope = preds[!preds %in% c("0", "1")]
      if (length(slope) > 1)
        stop("Error in segment ", i, " (linear): More than one variable in LHS of varying effect.")
      else if (length(slope) == 0)
        # If not slope
        slope = NA

      # Return
      new_row = tibble::tibble(
        int = !"0" %in% preds | parts[2] == "",  # bool. Is intercept present?
        slope = slope,
        group = parts[3])
      V = dplyr::bind_rows(V, new_row)
    }
  } else V = NA

  # Population-level intercepts
  population = attributes(stats::terms(nobars(form_rhs)))
  int_rel = population$term.labels == "rel(1)"
  population$term.labels = population$term.labels[!population$term.labels %in% "rel(1)"]  # code as no term
  if (any(int_rel))
    population$intercept = 1


  # Population-level slopes
  n_slopes = length(population$term.labels)
  if (n_slopes > 1) {
    stop("Error in segment ", i, " (linear): Only one slope allowed in population-level effects but found: '", paste0(population$term.labels, collapse="' and '"), "'")
  } else if (n_slopes == 1) {
    # One slope. Code as slope_rel = TRUE/FALSE and x = str
    slope_rel = stringr::str_detect(population$term.labels, "rel\\(")
    population$term.labels = gsub("rel\\(|\\)", "", population$term.labels)
    slope = population$term.labels
  } else if (n_slopes == 0) {
    slope = NA
    slope_rel = FALSE
  }

  # Return as list.
  return(tibble::tibble(
    int = population$intercept == 1,
    int_rel = any(int_rel),
    slope = slope,
    slope_rel = slope_rel,
    varying = list(V)
  ))
}


#' Build a table describing a list of segments
#'
#' Used internally for most mcp functions.
#'
#' @aliases get_segment_table
#' @inheritParams mcp
#' @return A tibble.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @importFrom magrittr %>%
#' @importFrom stats gaussian binomial
#' @export
#' @examples
#' \dontrun{
#' segments = list(
#'   y ~ 1 + x,
#'   1 ~ 1
#' )
#' get_segment_table(segments)
#' }

get_segment_table = function(segments, data = NULL, family = gaussian(), par_x = NULL) {
  #####################################################
  # BUILD "SEGMENT TABLE (ST)" FROM ISOLATED SEGMENTS #
  #####################################################
  ST = tibble::tibble()
  for (i in seq_len(length(segments))) {
    # Get ready...
    segment = segments[[i]]
    if (!is.null(data))
      check_terms_in_data(segment, data, i)

    # Go! Unpack this segment
    row = tibble::tibble(segment = i)
    row = dplyr::bind_cols(row, unpack_tildes(segment, i))
    row = dplyr::bind_cols(row, unpack_y(row$form_y, i, family = family))
    row = dplyr::bind_cols(row, unpack_cp(row$form_cp, i))
    row = dplyr::bind_cols(row, unpack_rhs(row$form_rhs, i))

    ST = dplyr::bind_rows(ST, row)
  }

  # Return the cols we need
  ST = dplyr::select(ST, -.data$form_y, -.data$form_cp, -.data$form_rhs) %>%
    tidyr::fill(.data$y, .data$trials, .direction="downup")  # Usually only provided in segment 1


  ###########################
  # CHECK SEGMENTS AND DATA #
  ###########################

  # Check segment 1: rel() not possible here.
  if (any(ST[1, c("cp_int", "cp_int_rel", "cp_ran_int", "cp_group_col")] != FALSE, na.rm = T))
    stop("Change point defined in first segment. This should not be possible. Submit bug report in the GitHub repo.")
  if (any(ST[1, c("int_rel", "slope_rel")] != FALSE, na.rm = T))
    stop("rel() cannot be used in segment 1. There is nothing to be relative to.")

  # Check rel() in segment 2+
  rel_slope_after_plateau = dplyr::lag(is.na(ST$slope), 1) & ST$slope_rel != 0
  if (any(rel_slope_after_plateau))
    stop("rel(slope) is not meaningful after a plateau segment (without a slope). Use absolute slope to get the same behavior. Found in segment ", which(rel_slope_after_plateau))
  if (nrow(ST) > 1) {
    if (ST$cp_int_rel[2] == TRUE)
      stop("rel() cannot be used for change points in segment 2. There are no earlier change points to be relative to. Relative changepoints work from segment 3 and on.")
  }

  # Set ST$x (what is the x-axis dimension?)
  derived_x = unique(stats::na.omit(ST$slope))
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
    if (all(is.na(ST$slope) & is.character(par_x)))
      ST$x = par_x
    else
      stop("This is a plateau-only model so no x-axis variable could be derived from the segment formulas. Use argument 'par_x' to set it explicitly")
  } else if (length(derived_x) > 1)
    # More than one...
    stop("More than one predictor found: '", paste0(unique(stats::na.omit(ST$slope)), collapse = "' and '"), "'")

  # Response variables
  derived_y = unique(stats::na.omit(ST$y))
  if (length(derived_y) != 1)
    stop("There should be exactly one response variable. Found '", paste0(derived_y, collapse="' and '", "'."))

  if (!is.na(ST$trials[1]) & family != "binomial")
    stop("response format `y | trials(N)` only meaningful for family = binomial(); not for ", family, "()")

  # Varying effects
  derived_varying = unique(stats::na.omit(ST$cp_group_col))

  # Check data types
  if (!is.null(data)) {
    # Convert to data.frame. Makes it easier to test column types.
    # Tibble will still be used in the rest of mcp
    if (tibble::is_tibble(data))
      data = data.frame(data)

    # Check x and y
    if (!is.numeric(data[, ST$x[1]]))
      stop("Data column '", ST$x[1], "' has to be numeric.")
    if (!is.numeric(data[, ST$y[1]]))
      stop("Data column '", ST$y[1], "' has to be numeric.")

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
    if (family == "binomial") {
      check_integer(data[, ST$y[1]], ST$y[1], lower = 0)
      check_integer(data[, ST$trials[1]], ST$trials[1], lower = 1)
    } else if (family == "bernoulli") {
      if (any(!data[, ST$y[1]] %in% c(0, 1)))
        stop("Only responses 0 and 1 are allowed for family = bernoulli() in column '", ST$y[1], "'")
    } else if (family == "poisson") {
      check_integer(data[, ST$y[1]], ST$y[1], lower = 0)
    }
  }


  ###################################
  # MUTATE ST WITH NEW HANDY VALUES #
  ###################################

  # Recode relative columns so 0 = not relative. N > 0 is the number of consecutive "relatives"
  ST = ST %>%

    # Add variable names
    dplyr::mutate(
      int_name = ifelse(.data$int, yes = paste0("int_", .data$segment), no = NA),
      slope_name = ifelse(!is.na(.data$slope), yes = paste0(.data$slope, "_", .data$segment), no = NA),
      slope_code = .data$slope_name,  # Will be modified in next step
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

    # Same for slope_code
    dplyr::group_by(cumsum(!.data$slope_rel)) %>%
    dplyr::mutate(
      slope_code = cumpaste(.data$slope_name, " + "),
      slope_code = format_code(.data$slope_code, na_col = .data$slope_name)
    ) %>%
    dplyr::ungroup()

  # Return
  ST
}


# Inspired by https://stackoverflow.com/questions/24862046/cumulatively-paste-concatenate-values-grouped-by-another-variable
cumpaste = function(x, .sep = " ")
  Reduce(function(x1, x2) paste(x1, x2, sep = .sep), x, accumulate = TRUE)


# Take a value like "a + b" and
# (1) replace it with NA if na_col == NA.
# (2) Change to "(a + b)" if there is a "+"
# (3) Return itself otherwise, e.g., "a" --> "a".
format_code = function(col, na_col) {
  # Replace with NA
  col = ifelse(is.na(na_col), NA, col)

  # Add parenthesis if this is an equation
  col = ifelse(stringr::str_detect(col, "\\+"), paste0("(", col, ")"), col)
  col
}


#' Throws an error if a number/vector contains non-numeric, decimal, or less-than-lower
#'
#' The expected behavior of is.integer, with informative error messages.
#'
#' @aliases check_integer
#' @param x Numeric value or vector
#' @param name Name to show in error message.
#' @param lower the smallest allowed value. lower = 1 checks for positive integers.
#'
check_integer = function(x, name, lower = -Inf) {
  if (!is.numeric(x))
    stop("Only integers >= ", lower, " allowed for '", name, "'")
  if (!all(x == floor(x)) | !all(x >= lower))  # any decimals or negative
    stop("Only integers >= ", lower, " allowed for '", name, "'")

  TRUE
}
