#' Detects par_x and verifies model-data fit
#'
#' @aliases get_par_x
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return The column name of par_x.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_par_x = function(model, data, par_x = NULL) {
  assert_types(model, "mcpmodel")
  assert_types(data, "data.frame")
  assert_types(par_x, "null", "character", len = c(0, 1))

  # Just check par_x
  if (is.character(par_x)) {
    if ((par_x %in% colnames(data)) == FALSE)
      stop("par_x = '", par_x, "' not found in data.")
    if (is_continuous(data[, par_x]) == FALSE)
      stop("par_x = '", par_x, "' has to be continuous. Is it binary or categorical?")
  }

  # Check for exactly one continuous
  rhs_vars = get_rhs_vars(model)
  data_in_rhs = data %>% dplyr::select(dplyr::all_of(rhs_vars), par_x)
  continuous_cols = lapply(data_in_rhs, is_continuous) %>% unlist()
  par_x_candidates = names(continuous_cols)[continuous_cols]
  if (is.character(par_x)) {
    if (length(par_x_candidates) == 0) {
      return(par_x)
    } else if (par_x %in% par_x_candidates) {
      return(par_x)
    } else {
      stop("Got par_x = '", par_x, "' but it does not seem to be both continuous, present in the data, and in the model. mcp identified '", par_x_candidates, "' as the only viable change point dimension(s) as the data and model is set up now.")
    }
  } else if (is.null(par_x)) {
    if (length(par_x_candidates) == 0) {
      stop("No continuous column for change points found in the formulas. Either provide mcp(..., par_x = 'my_col') or update the model.")
    } else if (length(par_x_candidates) > 1) {
      stop("Could not automatically determine the change point dimension (multiple candidates: ", and_collapse(par_x_candidates), "). Set it explictly using mcp(..., par_x = 'my_col').")
    } else if (length(par_x_candidates) == 1) {
      return(par_x_candidates)
    }
  } else {
    stop_github("Reached the end of get_par_x() without returning par_x")
  }
}


#' Get parameter table for a particular RHS dpar
#'
#' This function extracts an `par_x`-less design matrix.
#' `par_x` will be relative to the segment onset, so it will be multiplied in in the formula
#' (`jags_code` and `fit$simulate()`).
#'
#' @aliases get_rhs_table_dpar
#' @keywords internal
#' @inheritParams mcp
#' @param form The RHS formula for a particular dpar of a segment.
#' @param form_rhs The full RHS formula of a segment, including one or several `form`s.
#' @param segment Integer. The segment number
#' @param dpar One of `c("mu", "sigma", "ar")`.
#' @param order Currently only applies to `dpar == "ar"`.
#' @param check_rank Boolean. Whether to stop on rank deficiency.
#' @return A tibble with one row per model parameter and the columns
#'   - `dpar`: character.
#'   - `segment`: the segment number (positive integer).
#'   - `display_name`: user-facing parameter name used in summary functions.
#'   - `code_name`: parameter name used in JAGS and internally in mcp.
#'   - `par_type`: One of "Intercept", "dummy", or "slope". Used for setting priors and for change point indicator func.
#'   - `order`: positive integer or NA. Currently only relevant for `dpar == "ar"`.
#'   - `matrix_data`: column of the design matrix less the `par_x` term.
#'
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs_table_dpar = function(data, form_rhs, segment, dpar, par_x, order = NULL, check_rank = TRUE) {
  # EMpty segments return no rows
  if (all(as.character(form_rhs) == c("~", "0")))
    return(data.frame(dpar = character(0), segment = numeric(0)))

  assert_types(data, "data.frame", "tibble")
  assert_types(form_rhs, "formula", len = 2)
  assert_integer(segment, lower = 1, len = 1)
  assert_types(dpar, "character", len = 1)
  assert_types(par_x, "character", len = 1)
  assert_types(order, "null", "integer", len = c(0, 1))
  if (is.null(order) == FALSE)
    assert_integer(order, lower = 1)

  # Variable names for non-mu terms are prefixed with the term type.
  if (dpar == "mu") {
    dpar_prefix = ""
  } else {
    dpar_prefix = paste0(dpar, order, "_")
  }

  mat = stats::model.matrix(form_rhs, data)
  if (check_rank == TRUE)
    assert_rank(mat, segment, dpar)


  #######################
  # GET PARAMATER NAMES #
  #######################
  pars = colnames(mat)

  # Check that all contents with par_x within parantheses contain only a single term
  split_terms = lapply(pars, function(x) stringr::str_split(x, ":")) %>% unlist()
  split_contains_multiple_terms = split_terms %>%
    stringr::str_extract("(?<=\\().*(?=\\))") %>%
    stringr::str_detect("[+:*]")  # TO DO: will fail to detect I(x:b) because the string is already split by :
  split_contains_x = term_contains(par_x, split_terms)
  is_bad = split_contains_x & split_contains_multiple_terms
  if (any(stats::na.omit(is_bad) == TRUE))
    stop("mcp does not currently support 2+ terms within a formula function when one of them is par_x = '", par_x, "'. Found: ", and_collapse(split_terms[which(is_bad)]))

  # Replace I(...) with ...
  I_contents = stringr::str_extract(pars, "(?<=I\\().*(?=\\))")
  pars = stringr::str_replace(pars, "I\\(.*\\)", I_contents)

  # Replace (Intercept) with Intercept
  is_intercept = pars == "(Intercept)"
  intercept_name = ifelse(dpar == "mu", "Intercept", "")
  pars[is_intercept] = intercept_name

  # display_name
  display_name = gsub("\\(|\\)", "", pars)
  display_name = gsub("^", "E", display_name, fixed = TRUE)
  display_name = gsub("-", "M", display_name, fixed = TRUE)
  display_name = paste0(dpar_prefix, display_name, "_", segment)
  display_name = gsub("__", "_", display_name, fixed = TRUE)

  # code_name
  code_name = gsub("[: +]", "", display_name)

  # is_dummy
  is_dummy = apply(mat, 2, function(x) all(x %in% c(0, 1)))


  ################
  # GET X_FACTOR #
  ################

  # Detect terms with par_x and extract this multiplicative part
  all_patterns = c("x", "x\\^[\\+\\-0-9]+", "abs\\(x\\)", "sin\\(x\\)", "cos\\(x\\)", "tan\\(x\\)", "exp\\(x\\)", "log\\(x\\)", "sqrt\\(x\\)")
  pars_terms = lapply(all_patterns, extract_expr, pars, par_x)

  # Now convert par_x to "x"
  pattern_convert_to_x = gsub("x", par_x, "^x$|^x(?=\\^)|(?<=\\()x", fixed = TRUE)
  pars_terms = lapply(pars_terms, function(x) stringr::str_replace(x, pattern_convert_to_x, "x"))

  # Finally, multiply-merge to vector
  data_subterms = tibble::tibble(as.data.frame(pars_terms))
  colnames(data_subterms) = gsub("\\", "", all_patterns, fixed = TRUE)
  x_factor = tidyr::unite(data_subterms, "x_factor", sep = "*", na.rm = TRUE) %>% dplyr::pull(x_factor)  # same length as pars

  # Check exponent
  exponent = stringr::str_extract(dplyr::pull(data_subterms, 2), "(?<=\\^)[-0-9]+$")
  sapply(as.numeric(exponent), assert_integer, lower = 0, name = ". Got exponents in formula.")

  # Independent of x?
  is_independent_of_x = term_contains(par_x, pars) == FALSE
  if (any(x_factor[is_independent_of_x] != ""))
    stop_github("Internal mcp error: coded term as both dependent and independent of x")
  x_factor[is_independent_of_x] = "1"
  if (any(x_factor == ""))
    stop_github("Internal mcp error: did not code x_factor for term ", pars[x_factor == ""])



  ################################
  # COMPUTE X-LESS DESIGN MATRIX #
  ################################

  # Divide design matrix cols with x_factor.
  # Evaluate x_factor funcs on par_x and divide it out of the design matrix.
  x = data[, par_x]
  x_factor_local = gsub("1", "rep(1, length(x))", x_factor)  # Make intercepts ("1") have the correct dimension.
  mat_factor_x = eval(str2lang(paste0("as.matrix(data.frame(", paste0(x_factor_local, collapse = ", "), "))")))
  mat_without_x = mat / mat_factor_x
  mat_without_x[mat == 0 & mat_factor_x == 0] = 1  # 0 / 0 means "identityt", i.e., = 1.

  rhs_table = data.frame(
    dpar = dpar,
    segment = segment,
    display_name,
    code_name = code_name,
    par_type = dplyr::case_when(
      is_intercept == TRUE ~ "Intercept",
      is_dummy == TRUE ~ "dummy",
      TRUE ~ "slope"
    ),
    order = ifelse(is.null(order), NA, order),
    x_factor = x_factor,
    matrix_col = seq_len(ncol(mat_without_x)),
    stringsAsFactors = FALSE
  ) %>%
    # Add data
    dplyr::rowwise() %>%
    dplyr::mutate(
      matrix_data = list(mat_without_x[, .data$matrix_col])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$matrix_col)

  # Return
  rhs_table
}


#' Get expressions as they appear with/without interactions
#'
#' @aliases extract_expr
#' @keywords internal
#' @noRd
#' @param expr The expression to search for, e.g., "x" or "sin(x)". This is the needle.
#' @param pars This is the haystack.
#' @param par_x The parameter to substitute for x, e.g., "myvar".
#' @return A character vector of length `pars` with matches to `expr`.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
extract_expr = function(expr, pars, par_x) {
  regex_exact = gsub("x", par_x, gsub("expr", expr, "^expr$|^expr(?=:)|(?<=:)expr(?=:)|(?<=:)expr$", fixed = TRUE), fixed = TRUE)  # Alone or as part of an interaction. Prevents something like "this_is_not_x^2" being detected as "x^2"
  stringr::str_extract(pars, regex_exact)
}


#' Detect terms that contain a particular variable
#'
#' Finds it whether it's in an interaction, in an expression, etc. without false positives.
#'
#' @aliases term_contains
#' @keywords internal
#' @noRd
#' @param par_x The parameter to search for (character).
#' @param terms A character vector of terms.
#' @return A logical vector of length `terms`
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
term_contains = function(par_x, terms) {
  regex_contains_par_x = gsub("x", par_x, "^x$|^x[: +*^]|[: +*^\\(]x[: +*^\\)]|[: +*^]x$")
  stringr::str_detect(terms, regex_contains_par_x)
}


#' @aliases get_rhs_table_segment
#' @keywords internal
#' @noRd
#' @describeIn get_rhs_table_dpar Apply `get_rhs_table_dpar` to each formula in a segment
get_rhs_table_segment = function(form_rhs, segment, family, data, par_x, check_rank = TRUE) {
  assert_types(form_rhs, "formula", len = c(1, 3))
  assert_integer(segment, lower = 1, len = 1)
  assert_types(family, "mcpfamily")
  assert_types(data, "data.frame", "tibble")
  assert_types(par_x, "character", len = 1)

  # Get general format
  form_rhs = stats::as.formula(form_rhs)
  attrs = attributes(stats::terms(remove_terms(form_rhs, "varying")))
  term_labels = attrs$term.labels



  ######
  # MU #
  ######
  # Start by building it as a string: "mu(1 + x + ...)" to bring it into a compatible format
  mu_terms = term_labels[stringr::str_detect(attrs$term.labels, "sigma\\(|ar\\(") == FALSE]

  if (length(mu_terms > 0)) {
    mu_terms[1] = paste0(attrs$intercept, " + ", mu_terms[1])
    mu_term = paste0(mu_terms, collapse = " + ")  # for use in fit$model and in summary()
    mu_term = paste0("mu(", mu_term, ")")  # Get it in "standard" format
  } else {
    mu_term = paste0("mu(", attrs$intercept, ")")  # Plateau model: "mu(0)" or "mu(1)"
  }
  mu_form = get_term_content(mu_term)
  mu_pars = get_rhs_table_dpar(data, mu_form, segment, "mu", par_x, NULL, check_rank)



  #########
  # SIGMA #
  #########
  # Extract sigma terms
  sigma_term = term_labels[stringr::str_detect(term_labels, "sigma\\(")]  # Which terms?

  # If not specified, sigma_1 is implicit in segment 1.
  if (length(sigma_term) == 0 && family$family == "gaussian" && segment == 1) {
    sigma_form = ~1
    sigma_pars = get_rhs_table_dpar(data, sigma_form, segment, dpar = "sigma", par_x)
  } else if (length(sigma_term) > 0) {
    if (family$family != "gaussian")
      stop("sigma() is only meaningful for family = gaussian()")

    sigma_form = get_term_content(sigma_term)
    sigma_pars = get_rhs_table_dpar(data, sigma_form, segment, dpar = "sigma", par_x, NULL, check_rank)
  } else {
    sigma_pars = NULL
  }


  ######
  # AR #
  ######
  ar_term = term_labels[stringr::str_detect(term_labels, "ar\\(")]  # Extract terms
  ar_stuff = unpack_arma(ar_term)  # $order and $form_str
  ar_form = get_term_content(ar_stuff$form_str)

  # Populate each of these for each order of AR
  ar_pars = list()
  if (!is.na(ar_stuff$order)) {
    for (order in seq_len(ar_stuff$order)) {
      ar_pars = rbind(ar_pars, get_rhs_table_dpar(data, ar_form, segment, "ar", par_x, order, check_rank))
    }
  }



  ##########
  # RETURN #
  ##########
  rbind(
    mu_pars,
    sigma_pars,
    ar_pars
  )
}




#' Get formula inside a wrapper
#'
#' @aliases get_term_content
#' @keywords internal
#' @noRd
#' @param term E.g., "mu(1 + x)", "sigma(0 + I(x^2))", etc.
#' @return char formula with the content inside the brackets.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_term_content = function(term) {
  # Handle cases of no input or several inputs
  if (length(term) == 0) {
    return(NA)
  } else if (length(term) > 1) {
    #dpar = paste0(substr(term, 0, content_start), ")")
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


#' Unpack arma order and formula
#'
#' @aliases unpack_arma
#' @keywords internal
#' @noRd
#' @param form_str_in A character. These are allowed: "ar(number)" or "ar(number, formula)"
#' @return A list with `$order` and `$form_str` (e.g., "ar(formula)"). The formula is ar(1) if no formula is given
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
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
  assert_integer(order, form_str_in, lower = 1, len = 1)

  # GET FORMULA
  if (has_formula) {
    # If there is a formula, remove the order
    form_str = gsub(paste0(order, ", "), "", form_str_in)
  } else {
    # If there is no formula, return "ar(1)"
    dpar = substr(form_str_in, 1, order_start - 2)  # "ar"
    form_str = paste0(dpar, "(1)")
  }

  # Return
  list(
    order = order,
    form_str = form_str
  )
}


#' @aliases get_rhs_table
#' @keywords internal
#' @describeIn get_rhs_table_dpar Apply `get_rhs_table_segment` to all segments of a model.
get_rhs_table = function(model, data, family, par_x, check_rank = TRUE) {
  rhs = lapply(model, get_rhs)

  rhs_table = lapply(seq_along(rhs), function(segment) get_rhs_table_segment(rhs[[segment]], segment, family, data, par_x, check_rank)) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(.data$dpar, .data$segment) %>%
    dplyr::mutate(matrix_col = dplyr::row_number())

  # Code next_intercept: Which segment has the next intercept?
  # Strategy: (1) select one row for segments with intercepts for each dpar (filter)
  #           (2) save this segment number in the last segment that had an intercept (lag)
  #           (3) left-join this into rhs_table
  #           (4) fill downwards into intermittent segments without intercepts
  df_next_intercept = rhs_table %>%
    dplyr::arrange(.data$dpar, .data$order, .data$segment) %>%
    dplyr::group_by(.data$dpar, .data$order) %>%
    dplyr::filter(.data$par_type == "Intercept") %>%
    dplyr::mutate(next_intercept = as.integer(dplyr::lead(.data$segment))) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$dpar, .data$segment, .data$order, .data$next_intercept)

  # Return: left-join and fill-down. NA means "there is no next intercept-segment"
  rhs_table %>%
    dplyr::left_join(df_next_intercept, by = c("dpar", "segment", "order")) %>%
    dplyr::group_by(.data$dpar, .data$order) %>%
    tidyr::fill(.data$next_intercept, .direction = "down") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(next_intercept = dplyr::if_else(.data$segment >= .data$next_intercept, NA_integer_, .data$next_intercept))
}
