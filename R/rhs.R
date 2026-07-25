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
  data_in_rhs = data %>% dplyr::select(dplyr::all_of(rhs_vars), dplyr::all_of(par_x))
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
#' @param form_rhs The full RHS formula of a segment, including one or several `form`s.
#' @param segment Integer. The segment number
#' @param dpar A distributional parameter or an `ar`/`ma` component.
#' @param order Applies to `dpar %in% c("ar", "ma")`.
#' @param check_rank Boolean. Whether to stop on rank deficiency.
#' @return A tibble with one row per model parameter and the columns
#'   - `dpar`: character.
#'   - `segment`: the segment number (positive integer).
#'   - `matrix_name`: original column name from the model matrix. Used to
#'     diagnose collisions after parameter names are converted for JAGS.
#'   - `display_name`: user-facing parameter name used in summary functions.
#'   - `code_name`: parameter name used in JAGS and internally in mcp.
#'   - `par_type`: One of "Intercept", "dummy", or "slope". Used for setting priors and for change point indicator func.
#'   - `order`: positive integer or NA. Only relevant for `ar` and `ma`.
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

  # Check raw terms before model.matrix() evaluates expressions such as I(x:b).
  formula_terms = attr(stats::terms(form_rhs), "term.labels")
  contains_multiple_terms = formula_terms %>%
    stringr::str_extract("(?<=\\().*(?=\\))") %>%
    stringr::str_detect("[+:*]")
  contains_x = term_contains(par_x, formula_terms)
  is_bad = contains_x & contains_multiple_terms
  if (any(stats::na.omit(is_bad) == TRUE))
    stop("mcp does not currently support 2+ terms within a formula function when one of them is par_x = '", par_x, "'. Found: ", and_collapse(formula_terms[which(is_bad)]))

  mat = stats::model.matrix(form_rhs, data)
  if (check_rank == TRUE)
    assert_rank(mat, segment, dpar)


  #######################
  # GET PARAMATER NAMES #
  #######################
  pars = colnames(mat)
  matrix_name = pars

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
    matrix_name = matrix_name,
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
    dplyr::select(-"matrix_col")

  # Return
  rhs_table
}


#' Check that model terms have unique internal parameter names
#'
#' Formula punctuation is removed from model-matrix column names because the
#' resulting parameter names must be valid in JAGS. Distinct terms can
#' therefore occasionally produce the same name, for example `a:b` and `ab`.
#'
#' @keywords internal
#' @noRd
#' @param rhs_table A data frame returned by `get_rhs_table_dpar()`.
#' @return `rhs_table`, invisibly. Stops with an informative error on collision.
assert_unique_rhs_names = function(rhs_table) {
  duplicated_name = duplicated(rhs_table$code_name) |
    duplicated(rhs_table$code_name, fromLast = TRUE)

  if (!any(duplicated_name))
    return(invisible(rhs_table))

  collision_names = unique(rhs_table$code_name[duplicated_name])
  collision_lines = vapply(collision_names, function(code_name) {
    rows = rhs_table[rhs_table$code_name == code_name, , drop = FALSE]
    sources = paste0(
      "`", rows$matrix_name, "` (", rows$dpar,
      ", segment ", rows$segment, ")"
    )
    paste0("  `", code_name, "`: ", and_collapse(sources))
  }, character(1))

  stop(
    "Model terms produce the same parameter name:\n",
    paste0(collision_lines, collapse = "\n"),
    "\nRename one predictor column and refit."
  )
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

  # Formula wrappers belonging to distributional parameters are declared by
  # the family. AR and MA are separate model components because they carry an
  # order as well as a formula.
  model_dpars = family$dpar_specs$dpar[family$dpar_specs$dpar != "mu"]
  dpar_patterns = paste0("^", model_dpars, "\\(")
  arma_components = c("ar", "ma")
  arma_pattern = paste0("^(", paste0(arma_components, collapse = "|"), ")\\(")

  # Give a family-specific error when a recognized dpar wrapper is unavailable.
  used_dpar_wrappers = known_dpar_wrappers()[vapply(
    known_dpar_wrappers(),
    function(dpar) any(stringr::str_detect(term_labels, paste0("^", dpar, "\\("))),
    logical(1)
  )]
  unsupported_dpars = setdiff(used_dpar_wrappers, model_dpars)
  if (length(unsupported_dpars) > 0) {
    dpar_calls = paste0("`", unsupported_dpars, "()`", collapse = " and ")
    family_call = paste0(family$family, "()")
    stop(
      dpar_calls, " is not a distributional parameter for family = ", family_call, ". ",
      "See available parameters with `mcpfamily(", family_call, ")$dpars`."
    )
  }



  ######
  # MU #
  ######
  # Start by building it as a string: "mu(1 + x + ...)" to bring it into a compatible format
  is_dpar_term = rep(FALSE, length(term_labels))
  for (pattern in dpar_patterns)
    is_dpar_term = is_dpar_term | stringr::str_detect(term_labels, pattern)
  is_arma_term = stringr::str_detect(term_labels, arma_pattern)
  mu_terms = term_labels[!is_dpar_term & !is_arma_term]

  if (length(mu_terms > 0)) {
    mu_terms[1] = paste0(attrs$intercept, " + ", mu_terms[1])
    mu_term = paste0(mu_terms, collapse = " + ")  # for use in fit$model and in summary()
    mu_term = paste0("mu(", mu_term, ")")  # Get it in "standard" format
  } else {
    mu_term = paste0("mu(", attrs$intercept, ")")  # Plateau model: "mu(0)" or "mu(1)"
  }
  mu_form = get_term_content(mu_term)
  mu_pars = get_rhs_table_dpar(data, mu_form, segment, "mu", par_x, NULL, check_rank)



  #############################
  # DISTRIBUTIONAL PARAMETERS #
  #############################
  dpar_pars = list()
  for (dpar in model_dpars) {
    spec = get_dpar_spec(family, dpar)
    dpar_term = term_labels[stringr::str_detect(term_labels, paste0("^", dpar, "\\("))]

    # An implicit dpar receives an intercept in segment 1 and then continues
    # across later segments until the user supplies another dpar intercept.
    if (length(dpar_term) == 0 && spec$implicit && segment == 1) {
      dpar_form = ~1
      dpar_pars[[dpar]] = get_rhs_table_dpar(
        data, dpar_form, segment, dpar = dpar, par_x, NULL, check_rank
      )
    } else if (length(dpar_term) > 0) {
      dpar_form = get_term_content(dpar_term)
      dpar_pars[[dpar]] = get_rhs_table_dpar(
        data, dpar_form, segment, dpar = dpar, par_x, NULL, check_rank
      )
    }
  }


  #########
  # AR/MA #
  #########
  arma_pars = list()
  for (component in arma_components) {
    component_term = term_labels[stringr::str_detect(term_labels, paste0("^", component, "\\("))]
    component_stuff = unpack_arma(component_term)

    if (!is.na(component_stuff$order)) {
      component_form = get_term_content(component_stuff$form_str)
      # Expand one formula into a separate regression parameter for each lag.
      arma_pars[[component]] = lapply(
        seq_len(component_stuff$order),
        function(order) get_rhs_table_dpar(
          data, component_form, segment, component, par_x, order, check_rank
        )
      ) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(boundary = component_stuff$boundary)
    }
  }

  # AR and MA use the same transformed observation, so their boundary must be
  # shared within a segment. One explicitly supplied value applies to both.
  supplied_boundaries = unique(stats::na.omit(unlist(lapply(arma_pars, function(x) x$boundary))))
  if (length(supplied_boundaries) > 1)
    stop("ar() and ma() must use the same `boundary` within a segment.")
  
  # Most users need no boundary argument; resolve the common default here.
  segment_boundary = if (length(supplied_boundaries) == 1) supplied_boundaries else 0.1
  arma_pars = lapply(arma_pars, dplyr::mutate, boundary = segment_boundary)



  ##########
  # RETURN #
  ##########
  dplyr::bind_rows(
    mu_pars,
    dplyr::bind_rows(dpar_pars),
    dplyr::bind_rows(arma_pars)
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
#' @param form_str_in A character such as `"ar(number)"`, `"ma(number)"`, or
#'   either component with a second formula argument.
#' @return A list with `$order`, `$form_str` (e.g., `"ar(formula)"`), and
#'   `$boundary`. The component formula is 1 if no formula is given.
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

  component_call = str2lang(form_str_in)
  component = as.character(component_call[[1]])
  assert_value(component, allowed = c("ar", "ma"))

  component_args = as.list(component_call)[-1]
  component_arg_names = names(component_args)
  if (is.null(component_arg_names))
    component_arg_names = rep("", length(component_args))

  if (length(component_args) == 0 || component_arg_names[1] %notin% c("", "order"))
    stop("The first argument to ", component, "() must be its order.")

  boundary_index = which(component_arg_names == "boundary")
  if (length(boundary_index) > 1)
    stop("Only one `boundary` value is allowed in ", component, "().")

  formula_index = setdiff(seq_along(component_args), c(1, boundary_index))
  if (length(formula_index) > 1 || any(component_arg_names[formula_index] %notin% c("", "formula")))
    stop(component, "() accepts only `order`, an optional formula, and `boundary`.")

  # GET ORDER
  order_str = paste(deparse(component_args[[1]], width.cutoff = 500), collapse = "")
  order = suppressWarnings(as.numeric(order_str))

  # Check the order
  if (is.na(order))
    stop("Wrong specification of order in '", form_str_in, "'. Must be ", component, "(order) or ", component, "(order, formula) where order is a positive integer.")
  assert_integer(order, form_str_in, lower = 1, len = 1)

  # GET FORMULA AND BOUNDARY
  if (length(formula_index) == 1) {
    formula_str = paste(deparse(component_args[[formula_index]], width.cutoff = 500), collapse = "")
    form_str = paste0(component, "(", formula_str, ")")
  } else {
    # If there is no formula, use an intercept-only component formula.
    form_str = paste0(component, "(1)")
  }

  if (length(boundary_index) == 1) {
    boundary_str = paste(deparse(component_args[[boundary_index]], width.cutoff = 500), collapse = "")
    boundary = suppressWarnings(as.numeric(boundary_str))
    if (length(boundary) != 1 || is.na(boundary) || !is.finite(boundary) || boundary <= 0 || boundary >= 1)
      stop("`boundary` in ", component, "() must be one number between 0 and 1.")
  } else {
    # Defer the default until AR and MA have been parsed together for a segment.
    boundary = NA_real_
  }

  # Return
  list(
    order = order,
    form_str = form_str,
    boundary = boundary
  )
}


#' @aliases get_rhs_table
#' @keywords internal
#' @describeIn get_rhs_table_dpar Apply `get_rhs_table_segment` to all segments of a model.
get_rhs_table = function(model, data, family, par_x, check_rank = TRUE) {
  family = add_dpar_specs(family)
  rhs = lapply(model, get_rhs)

  rhs_table = lapply(seq_along(rhs), function(segment) get_rhs_table_segment(rhs[[segment]], segment, family, data, par_x, check_rank)) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(.data$dpar, .data$segment) %>%
    dplyr::mutate(matrix_col = dplyr::row_number())
  if ("boundary" %notin% names(rhs_table))
    rhs_table$boundary = NA_real_

  assert_unique_rhs_names(rhs_table)

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
    dplyr::select("dpar", "segment", "order", "next_intercept")

  # Return: left-join and fill-down. NA means "there is no next intercept-segment"
  rhs_table %>%
    dplyr::left_join(df_next_intercept, by = c("dpar", "segment", "order")) %>%
    dplyr::group_by(.data$dpar, .data$order) %>%
    tidyr::fill("next_intercept", .direction = "down") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(next_intercept = dplyr::if_else(.data$segment >= .data$next_intercept, NA_integer_, .data$next_intercept))
}
