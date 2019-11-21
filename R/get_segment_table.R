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

# Checks if all terms are in the data
check_terms_in_data = function(form, data, i) {
  found = all.vars(form) %in% colnames(data)
  found = found[stringr::str_starts(found, "sigma")]  # Sigma need not be in data
  if (!all(found))
    stop("Error in segment ", i, ": Term '", paste0(all.vars(form)[!found], collapse="' and '"), "' found in formula but not in data.")
}


# Unpacks y variable name
unpack_y = function(form_y, i, family) {
  # If NA and not segment 1, just return empty
  if (is.na(form_y)) {
    if ( i == 1)
      stop("A response must be defined in segment 1, e.g., 'y ~ 1'")

    return(tibble::tibble(
      y = NA,
      trials = NA
    ))
  }

  # Codings for binomial and variance-change
  form_y = gsub(" ", "", form_y)  # remove whitespace
  got_trials = stringr::str_detect(form_y, "\\|trials\\(")

  # Segment 1 must have y specified (and correctly so)
  if (i == 1) {
    if (family == "binomial" & !got_trials)
      stop("Error in response of segment 1: need a valid trials() specification, e.g., 'y | trials(N) ~ 1 + x'")

    if (family != "binomial" & got_trials)
      stop("y | trials(N) not meaningful for non-binomial models.")
  }


  # Split the response into its constituents, if there are any
  if (got_trials) {
    # If binomial, split into y (chunks[1]) and trials (chunks[2])
    chunks = gsub(")", "", strsplit(form_y, "\\|trials\\(")[[1]])
  } else {
    chunks = form_y
  }

  # Check that the constituents have valid names
  if (!all(grepl("^[A-Za-z._0-9]+$|\\+|\\-|\\/|\\*", chunks)))
    stop("Error in segment ", i, ": Invalid format for response variable.")

  # Return!
  return(
    tibble::tibble(
      y = chunks[1],  # Char
      trials = ifelse(got_trials == TRUE, chunks[2], NA)  # Char or NA
    )
  )
}


#' Takes a cp formula (as a string) and returns its properties
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



unpack_rhs = function(form_rhs, i, family, data, last_segment) {
  form_rhs = formula(form_rhs)
  population = attributes(stats::terms(nobars(form_rhs)))

  #############
  # INTERCEPT #
  #############
  int_rel = population$term.labels == "rel(1)"
  if (any(int_rel) == TRUE) {
    # Not a term. Remove
    population$term.labels = population$term.labels[!int_rel]

    if (i == 1)
      stop("rel() cannot be used in segment 1. There is nothing to be relative to.")
  }

  if (population$intercept == TRUE | any(int_rel)) {
    int = tibble::tibble(
      name = paste0("int_", i),
      rel = any(int_rel)
    )
  } else {
    int = NA
  }


  #########
  # SIGMA #
  #########
  # The sigma code is before the "slope" code, because sigma terms must be
  # removed so as to not be coded as slopes on the mean.
  sigma_terms_index = stringr::str_detect(population$term.labels, "sigma\\(")  # Which terms are sigma?
  sigma_terms_str = population$term.labels[sigma_terms_index]  # Extract sigma terms
  population$term.labels = population$term.labels[!sigma_terms_index]  # Remove sigma from terms

  if (length(sigma_terms_str) > 0) {
    # Extract the contents of (several) sigma(x) into one formula and get it's attributes
    sigma_formula_str = paste0(gsub("^sigma\\(|\\)$", "", sigma_terms_str), collapse = "+")
    sigma_formula = as.formula(paste0("~", sigma_formula_str))
    sigma_attributes = attributes(terms(sigma_formula))
    sigma_terms = sigma_attributes$term.labels

    # INTERCEPT
    # Sigma intercept (this is basically copy-paste of intercept code above)
    sigma_int_rel = sigma_terms == "rel(1)"
    if (any(sigma_int_rel) == TRUE) {
      # Not a term. Remove
      sigma_terms = sigma_terms[!sigma_int_rel]
      if (i == 1)
        stop("rel() cannot be used in segment 1. There is nothing to be relative to.")
    }

    # Got an intercept!
    if (sigma_attributes$intercept == TRUE | any(sigma_int_rel)) {
      sigma_int = tibble::tibble(
        name = paste0("sigma_", i),
        rel = any(sigma_int_rel)
      )
    }

    # SLOPE
    # Sigma slope (this is basically copy-paste of slope code below)
    if (length(sigma_terms) > 0) {
      sigma_slope = lapply(
        sigma_terms,
        FUN = unpack_slope_term,
        i = i,
        last_slope = last_segment$sigma_slope,
        prefix = "sigma_") %>%
        dplyr::bind_rows()  # Make it a proper table-like tibble

      # Build code. Add parentheses if there are more slopes
      sigma_code = paste(sigma_slope$code, collapse = " + ")
      if (stringr::str_detect(sigma_code, "\\+"))
        sigma_code = paste0("(", sigma_code, ")")
    }
  }

  # Sigma was not detected in form_rhs. But it is implicitly started in segment 1.
  if(i == 1 & family == "gaussian") {
    sigma_int = tibble::tibble(
      name = "sigma_1",
      code = "sigma_1",
      rel = FALSE
    )
  }

  # The known unknown :-)
  if(!exists("sigma_int"))
    sigma_int = NA
  if(!exists("sigma_slope"))
    sigma_slope = NA
  if(!exists("sigma_code"))
    sigma_code = NA



  ############
  # SLOPE(S) #
  ############
  # Population-level slopes
  if (length(population$term.labels) > 0) {
    slope = lapply(
      population$term.labels,
      FUN = unpack_slope_term,
      i = i,
      last_slope = last_segment$slope) %>%
      dplyr::bind_rows()  # Make it a proper table-like tibble

    # # Checks
    # par_x = unique(slope$par_x)  # What is the data predictor?
    # if (length(par_x) > 1)
    #   stop("Found more than one x-variable in segment ", i, ": '", paste(slope$par_x, collapse ="', '"), "'")

    # Build code. Add parentheses if there are more slopes
    slope_code = paste(slope$code, collapse = " + ")
    if (stringr::str_detect(slope_code, "\\+"))
      slope_code = paste0("(", slope_code, ")")

  } else {
    slope = NA
    slope_code = NA
  }


  #################
  # EXTRACT PAR_X #
  #################
  par_x = stats::na.omit(unique(c(
    ifelse(!is.na(slope), slope$par_x, NA),
    ifelse(!is.na(sigma_slope), sigma_slope$par_x, NA)
  )))
  if (length(par_x) > 1) {
    stop("Found more than one x-variable in segment ", i, ": '", paste(par_x, collapse ="', '"), "'")
  } else if (length(par_x) == 0) {
    par_x = NA
  }


  ###########
  # VARYING #
  ###########
  varying_terms = as.character(findbars(form_rhs))
  if (length(varying_terms) > 0) {
    V = lapply(varying_terms, unpack_varying_term, i = i) %>%
      dplyr::bind_rows()
  } else {
    V = NA
  }


  # Return as list.
  return(tibble::tibble(
    int = list(int),

    slope = list(slope),
    slope_code = slope_code,

    sigma_int = list(sigma_int),
    sigma_slope = list(sigma_slope),
    sigma_code = sigma_code,

    x = par_x,
    varying = list(V)
  ))
}


unpack_slope_term = function(term, i, last_slope, prefix = "") {
  # Remove the "rel()" and "I()" in that order.
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
    if(is.null(last_slope))
      stop("rel(", term, ") in segment ", i, " but no slopes in the previous segment.")
    if(!term %in% last_slope[[1]]$term) {
      stop("rel(", term, ") in segment ", i, " does not have a corresponding term to be relative to in the previous segment.")
    }
  }
  if(stringr::str_detect(term, "^log\\(|^sqrt\\(") & i > 1)
    stop("log() or sqrt() detected in segment 2+. This would fail because mcp models earlier segments as negative x values, and sqrt()/log() cannot take negative values.")

  # Regular expressions. Only recognize stuff that is identical between JAGS and base R
  func_list = c("abs", "cos", "exp", "log", "sin", "sqrt", "tan")
  funcs_regex = paste0("^", func_list, "\\(", collapse = "|")  # ^func1(|^func2(|...
  exponent_regex = "\\^[0-9.]+$"  # something^[number]
  end_regex = "\\)$"

  # Find par_x by removing everything associated with accepted functions
  par_x = gsub(paste0(c(funcs_regex, exponent_regex, end_regex), collapse = "|"), "", term)

  # Give exponents a valid variable name
  rel_x_code = paste0("X_", i, "_[i_]")  # x relative to segment start. Must match that inserted in get_formula()
  if (stringr::str_detect(term, exponent_regex)) {
    # Exponential
    slope_base = gsub("\\^", "E", term)
    term_recode = gsub(paste0("^", par_x, "\\^"), paste0(rel_x_code, "\\^"), term)
  } else if (stringr::str_detect(term, funcs_regex)) {
    # A simple function
    slope_base = gsub("\\(", "_", term)  # Replace first parenthesis with underscore
    slope_base = gsub("\\)$", "", slope_base)  # Remove second (last) parenthesis
    term_recode = gsub(paste0("\\(", par_x), paste0("(", rel_x_code), term)
  } else if (term == par_x) {
    # No function; just vanilla :-)
    slope_base = par_x
    term_recode = rel_x_code
  } else {
    stop("mcp failed for term ", term, ". If there is no obvious reason why, please raise an issue at GitHub.")
  }

  # Add last segment if this slope is relative
  name = paste0(prefix, slope_base, "_", i)
  if (rel == FALSE) {
    name_cumul = name
    code = paste0(name, " * ", term_recode)
  } else if (rel == TRUE) {
    last_name_cumul =last_slope[[1]]$name_cumul[which(last_slope[[1]]$term == term)]
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
#'   ~ 1
#' )
#' get_segment_table(segments)
#' }

get_segment_table = function(segments, data = NULL, family = gaussian()$family, par_x = NULL) {
  #####################################################
  # BUILD "SEGMENT TABLE (ST)" FROM ISOLATED SEGMENTS #
  #####################################################
  ST = tibble::tibble()
  for (i in seq_along(segments)) {
    # Get ready...
    segment = segments[[i]]
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
    tidyr::fill(.data$y, .data$form_y, .data$trials, .direction="downup") %>%  # Usually only provided in segment 1
    dplyr::mutate(form = paste0(.data$form_y, ifelse(segment == 1, "", .data$form_cp), .data$form_rhs)) %>%  # build full formula
    dplyr::select(-.data$form_y, -.data$form_cp, -.data$form_rhs)  # Not needed anymore




  ###########################
  # CHECK SEGMENTS AND DATA #
  ###########################

  # Check segment 1: rel() not possible here.
  if (any(ST[1, c("cp_int", "cp_int_rel", "cp_ran_int", "cp_group_col")] != FALSE, na.rm = T))
    stop("Change point defined in first segment. This should not be possible. Submit bug report in the GitHub repo.")
  #if (any(ST[1, c("int_rel", "slope_rel")] != FALSE, na.rm = T))
  #  stop("rel() cannot be used in segment 1. There is nothing to be relative to.")

  # Check rel() in segment 2+
  #rel_slope_after_plateau = dplyr::lag(is.na(ST$slope), 1) & ST$slope_rel != 0
  #if (any(rel_slope_after_plateau))
  #  stop("rel(slope) is not meaningful after a plateau segment (without a slope). Use absolute slope to get the same behavior. Found in segment ", which(rel_slope_after_plateau))
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
    stop("There should be exactly one response variable. Found '", paste0(derived_y, collapse="' and '", "'."))

  if (!is.na(ST$trials[1]) & family != "binomial")
    stop("Response format `y | trials(N)` only meaningful for family = binomial(); not for ", family, "()")

  # Varying effects
  derived_varying = unique(stats::na.omit(ST$cp_group_col))

  # Sigma
  if (any(c(!is.na(ST$sigma_int)), !is.na(ST$sigma_slope)) & family != "gaussian")
    stop("sigma() is only meaningful for family = gaussian()")

  # Check data types
  if (!is.null(data)) {
    # Convert to data.frame. Makes it easier to test column types.
    # Tibble will still be used in the rest of mcp
    if (tibble::is_tibble(data))
      data = data.frame(data)

    # Check x and y
    if (!is.numeric(dplyr::pull(data, ST$x[1])))
      stop("Data column '", ST$x[1], "' has to be numeric.")
    if (!is.numeric(dplyr::pull(data, ST$y[1])))
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
      #int_name = ifelse(.data$int, yes = paste0("int_", .data$segment), no = NA),
      #slope_name = ifelse(!is.na(.data$slope), yes = paste0(.data$slope, "_", .data$segment), no = NA),
      #slope_code = .data$slope_name,  # Will be modified later
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
    # dplyr::group_by(cumsum(!.data$slope_rel)) %>%
    # dplyr::mutate(
    #   slope_code = cumpaste(.data$slope_name, " + "),
    #   slope_code = format_code(.data$slope_code, na_col = .data$slope_name)
    # ) %>%
    # dplyr::ungroup() %>%

    # Finish up
    dplyr::select(-dplyr::starts_with("cumsum"))

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
  greater_than = ifelse(lower == -Inf, " ", paste0(" >= ", lower, " "))
  if (!is.numeric(x))
    stop("Only integers", greater_than, "allowed for '", name, "'")
  if (!all(x == floor(x)) | !all(x >= lower))
    stop("Only integers", greater_than, "allowed for '", name, "'")

  TRUE
}
