# ABOUT: Construction of the population predictor design matrix: parsing a
# segment's predictor formula, evaluating `model.matrix()` per distributional
# parameter, and assembling the resulting per-segment and per-model tables.
# ------------------------------------------------------------

#' Build and retain an R model-matrix specification
#'
#' The terms stored on a model frame contain fitted calls produced by
#' `makepredictcall()`, including centers, scales, polynomial coefficients,
#' and spline knots. Reusing them makes design matrices independent of
#' `newdata` and later changes to contrast options.
#'
#' @keywords internal
#' @noRd
#' @param form A one-sided predictor formula.
#' @param data A data frame used to fit the design.
#' @param spec An optional fitted specification returned by this function.
#' @return A list with `matrix` and `spec`.
get_fitted_design = function(form = NULL, data, spec = NULL) {
  if (is.null(spec)) {
    frame = stats::model.frame(form, data)
    fitted_terms = attr(frame, "terms")
    matrix = stats::model.matrix(fitted_terms, frame)
    factor_cols = vapply(frame, is.factor, logical(1))
    spec = list(
      terms = fitted_terms,
      xlevels = lapply(frame[factor_cols], levels),
      contrasts = attr(matrix, "contrasts"),
      columns = colnames(matrix)
    )
  } else {
    frame = stats::model.frame(spec$terms, data, xlev = spec$xlevels)
    matrix = stats::model.matrix(
      spec$terms, frame, contrasts.arg = spec$contrasts
    )
    if (!identical(colnames(matrix), spec$columns))
      stop("The model matrix for `newdata` does not match the fitted model.")
  }

  list(matrix = matrix, spec = spec)
}


#' Collect fitted design specifications carried by predictor rows
#'
#' Each call to `get_predictors_dpar()` temporarily repeats its fitted design
#' specification on its coefficient rows so ordinary `bind_rows()` operations
#' can carry it through the parser. This helper deduplicates those temporary
#' columns into the named list stored once on the fitted model.
#'
#' @keywords internal
#' @noRd
#' @param ... Predictor tables that may contain `design_id` and `design_spec`.
#' @return A named list of fitted design specifications.
collect_design_specs = function(...) {
  tables = list(...)
  tables = lapply(tables, function(table) {
    if (!all(c("design_id", "design_spec") %in% names(table)))
      return(NULL)
    dplyr::select(table, "design_id", "design_spec")
  })
  rows = dplyr::bind_rows(tables) %>%
    dplyr::filter(!is.na(.data$design_id)) %>%
    dplyr::distinct(.data$design_id, .keep_all = TRUE)

  stats::setNames(rows$design_spec, rows$design_id)
}


#' Express a design matrix relative to the segment onset
#'
#' @keywords internal
#' @noRd
#' @param matrix A model matrix.
#' @param x_factor One factor for each matrix column.
#' @param x Values of the change-point predictor.
#' @return `matrix` with the local change-point factors divided out.
remove_x_factors = function(matrix, x_factor, x) {
  factor_code = x_factor
  factor_code[factor_code == "1"] = "rep(1, length(x))"
  factor_matrix = eval(str2lang(paste0(
    "as.matrix(data.frame(", paste0(factor_code, collapse = ", "), "))"
  )))
  out = matrix / factor_matrix
  out[matrix == 0 & factor_matrix == 0] = 1
  out
}

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
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")
  checkmate::assert_data_frame(data)
  checkmate::assert_string(par_x, null.ok = TRUE)

  # Just check par_x
  if (is.character(par_x)) {
    if ((par_x %in% colnames(data)) == FALSE)
      stop("par_x = '", par_x, "' not found in data.")
    if (is_continuous(data[, par_x]) == FALSE)
      stop("par_x = '", par_x, "' has to be continuous. Is it binary or categorical?")
  }

  # Check for exactly one continuous
  rhs_vars = setdiff(get_rhs_vars(model), get_rhs_group_vars(model))
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


#' Get predictors for one distributional parameter
#'
#' This function extracts an `par_x`-less design matrix.
#' `par_x` will be relative to the segment onset, so it will be multiplied in in the formula
#' (`jags_code` and `fit$simulate()`).
#'
#' @aliases get_predictors_dpar
#' @keywords internal
#' @inheritParams mcp
#' @param form_rhs The full predictor formula of a segment, including one or
#'   several distributional terms.
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
#'   - `explicit`: whether the distributional parameter was supplied in the formula.
#'   - `design_id`: key of the fitted component formula that produced the row.
#'   - `design_col`: column occupied by the row in that component's model matrix.
#'   - `matrix_data`: column of the design matrix less the `par_x` term.
#'
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_predictors_dpar = function(data, form_rhs, segment, dpar, par_x, order = NULL,
                               check_rank = TRUE, design_id = NULL) {
  # EMpty segments return no rows
  if (all(as.character(form_rhs) == c("~", "0")))
    return(tibble::tibble(
      dpar = character(),
      segment = integer(),
      matrix_name = character(),
      display_name = character(),
      code_name = character(),
      par_type = character(),
      order = integer(),
      x_factor = character(),
      design_id = character(),
      design_col = integer(),
      design_spec = list(),
      matrix_data = list()
    ))

  checkmate::assert_data_frame(data)
  checkmate::assert_formula(form_rhs)
  checkmate::assert_int(length(form_rhs), lower = 2, upper = 2, .var.name = "length(form_rhs)")
  checkmate::assert_int(segment, lower = 1)
  checkmate::assert_string(dpar)
  checkmate::assert_string(par_x)
  checkmate::assert_string(design_id)
  checkmate::assert_integer(order, max.len = 1, null.ok = TRUE)
  if (is.null(order) == FALSE)
    checkmate::assert_integerish(order, lower = 1)

  # Variable names for non-mu terms are prefixed with the term type.
  if (dpar == "mu") {
    dpar_prefix = ""
  } else {
    dpar_prefix = paste0(dpar, order, "_")
  }

  # Complex expressions involving par_x do not yet have stable parameter names.
  formula_terms = attr(stats::terms(form_rhs), "term.labels")
  contains_multiple_terms = formula_terms %>%
    stringr::str_extract("(?<=\\().*(?=\\))") %>%
    stringr::str_detect("[+:*]")
  contains_x = term_contains(par_x, formula_terms)
  is_bad = contains_x & contains_multiple_terms
  if (any(stats::na.omit(is_bad) == TRUE))
    stop("mcp does not currently support 2+ terms within a formula function when one of them is par_x = '", par_x, "'. Found: ", and_collapse(formula_terms[which(is_bad)]))

  design = get_fitted_design(form_rhs, data)
  mat = design$matrix
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

  # Multi-column bases often include arguments and namespaces in their matrix
  # names. Name both the position within the basis and the segment explicitly.
  is_basis = grepl("[,=]|::", matrix_name)
  term_index = attr(mat, "assign")[is_basis]
  basis_name = formula_terms[term_index]
  basis_name = sub("^.*::", "", basis_name)
  basis_name = gsub("[^A-Za-z0-9]+", "_", basis_name)
  basis_name = gsub("^_|_$", "", basis_name)
  basis_col = ave(term_index, term_index, FUN = seq_along)
  display_name[is_basis] = paste0(
    dpar_prefix, basis_name, "_basis", basis_col, "_", segment
  )

  # code_name
  code_name = gsub("[: +]", "", display_name)

  # is_dummy
  is_dummy = apply(mat, 2, function(x) all(x %in% c(0, 1)))


  ################
  # GET X_FACTOR #
  ################

  # Bare par_x and polynomial bases are relative to the segment onset. Other
  # transformations stay evaluated on the original input in the model matrix.
  local_patterns = c("x", "x\\^[\\+\\-0-9]+")
  pars_terms = lapply(local_patterns, extract_expr, pars, par_x)

  # Now convert par_x to "x"
  pattern_convert_to_x = gsub("x", par_x, "^x$|^x(?=\\^)|(?<=\\()x", fixed = TRUE)
  pars_terms = lapply(pars_terms, function(x) stringr::str_replace(x, pattern_convert_to_x, "x"))

  # Finally, multiply-merge to vector
  data_subterms = tibble::tibble(as.data.frame(pars_terms))
  colnames(data_subterms) = gsub("\\", "", local_patterns, fixed = TRUE)
  x_factor = tidyr::unite(data_subterms, "x_factor", sep = "*", na.rm = TRUE) %>% dplyr::pull(x_factor)  # same length as pars

  # Check exponent
  exponent = stringr::str_extract(dplyr::pull(data_subterms, 2), "(?<=\\^)[-0-9]+$")
  checkmate::assert_integerish(
    as.numeric(exponent),
    lower = 0,
    .var.name = "exponents in formula"
  )

  # Everything not recognized as a local basis remains in the model matrix.
  x_factor[x_factor == ""] = "1"



  ################################
  # COMPUTE X-LESS DESIGN MATRIX #
  ################################

  # Divide design matrix cols with x_factor.
  # Evaluate x_factor funcs on par_x and divide it out of the design matrix.
  mat_without_x = remove_x_factors(mat, x_factor, data[, par_x])

  predictors = data.frame(
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
    design_id = design_id,
    design_col = seq_len(ncol(mat_without_x)),
    matrix_col = seq_len(ncol(mat_without_x)),
    stringsAsFactors = FALSE
  ) %>%
    # Add data
    dplyr::rowwise() %>%
    dplyr::mutate(
      design_spec = list(design$spec),
      matrix_data = list(mat_without_x[, .data$matrix_col])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-"matrix_col")

  # Return
  predictors
}


#' Check that model terms have unique internal parameter names
#'
#' Formula punctuation is removed from model-matrix column names because the
#' resulting parameter names must be valid in JAGS. Distinct terms can
#' therefore occasionally produce the same name, for example `a:b` and `ab`.
#'
#' @keywords internal
#' @noRd
#' @param predictors A data frame returned by `get_predictors_dpar()`.
#' @return `predictors`, invisibly. Stops with an informative error on
#'   collision.
assert_unique_predictor_names = function(predictors) {
  duplicated_name = duplicated(predictors$code_name) |
    duplicated(predictors$code_name, fromLast = TRUE)

  if (!any(duplicated_name))
    return(invisible(predictors))

  collision_names = unique(predictors$code_name[duplicated_name])
  collision_lines = vapply(collision_names, function(code_name) {
    rows = predictors[predictors$code_name == code_name, , drop = FALSE]
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
#' @param expr The expression to search for, e.g., "x" or "I(x^2)". This is the needle.
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


#' @aliases get_predictors_segment
#' @keywords internal
#' @noRd
#' @describeIn get_predictors_dpar Apply `get_predictors_dpar` to
#'   each formula in a segment
get_predictors_segment = function(form_rhs, segment, family, data, par_x, check_rank = TRUE) {
  checkmate::assert_formula(form_rhs)
  checkmate::assert_int(segment, lower = 1)
  checkmate::assert_true(is.mcpfamily(family), .var.name = "family")
  checkmate::assert_data_frame(data)
  checkmate::assert_string(par_x)

  # Get general format. Top-level group terms belong to mu; group terms inside
  # distributional wrappers are removed from those formulas below.
  form_rhs = stats::as.formula(form_rhs)
  attrs = attributes(stats::terms(form_rhs))
  term_labels = attrs$term.labels
  top_level_group = vapply(term_labels, is_group_term, logical(1))
  term_labels = term_labels[!top_level_group]

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
  mu_pars = get_predictors_dpar(
    data, mu_form, segment, "mu", par_x, NULL, check_rank,
    design_id = paste("population", "mu", segment, sep = ":")
  ) %>%
    dplyr::mutate(explicit = TRUE)



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
      dpar_pars[[dpar]] = get_predictors_dpar(
        data, dpar_form, segment, dpar = dpar, par_x, NULL, check_rank,
        design_id = paste("population", dpar, segment, sep = ":")
      ) %>%
        dplyr::mutate(explicit = FALSE)
    } else if (length(dpar_term) > 0) {
      dpar_form = get_term_content(dpar_term)
      dpar_form = remove_terms(dpar_form, "varying")
      dpar_pars[[dpar]] = get_predictors_dpar(
        data, dpar_form, segment, dpar = dpar, par_x, NULL, check_rank,
        design_id = paste("population", dpar, segment, sep = ":")
      ) %>%
        dplyr::mutate(explicit = TRUE)
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
      if (length(get_group_terms(component_form)) > 0)
        stop(
          "Group-level effects inside ", component,
          "() are not currently supported. Found one in segment ", segment, "."
        )
      # Expand one formula into a separate regression parameter for each lag.
      arma_pars[[component]] = lapply(
        seq_len(component_stuff$order),
        function(order) get_predictors_dpar(
          data, component_form, segment, component, par_x, order, check_rank,
          design_id = paste("population", component, order, segment, sep = ":")
        )
      ) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(boundary = component_stuff$boundary, explicit = TRUE)
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


#' @aliases get_predictor_tables
#' @keywords internal
#' @describeIn get_predictors_dpar Apply `get_predictors_segment`
#'   to all segments of a model.
get_predictor_tables = function(model, data, family, par_x, check_rank = TRUE) {
  rhs = lapply(model, get_rhs)

  predictors = lapply(seq_along(rhs), function(segment) get_predictors_segment(rhs[[segment]], segment, family, data, par_x, check_rank)) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(.data$dpar, .data$segment) %>%
    dplyr::mutate(matrix_col = dplyr::row_number())
  if ("boundary" %notin% names(predictors))
    predictors$boundary = rep(NA_real_, nrow(predictors))

  assert_unique_predictor_names(predictors)

  # Code next_intercept: Which segment has the next intercept?
  # Strategy: (1) select one row for segments with intercepts for each dpar (filter)
  #           (2) save this segment number in the last segment that had an intercept (lag)
  #           (3) left-join this into predictors
  #           (4) fill downwards into intermittent segments without intercepts
  df_next_intercept = predictors %>%
    dplyr::arrange(.data$dpar, .data$order, .data$segment) %>%
    dplyr::group_by(.data$dpar, .data$order) %>%
    dplyr::filter(.data$par_type == "Intercept") %>%
    dplyr::mutate(next_intercept = as.integer(dplyr::lead(.data$segment))) %>%
    dplyr::ungroup() %>%
    dplyr::select("dpar", "segment", "order", "next_intercept")

  # Population predictors: left-join and fill-down. NA means "there is no next
  # intercept segment".
  predictors = predictors %>%
    dplyr::left_join(df_next_intercept, by = c("dpar", "segment", "order")) %>%
    dplyr::group_by(.data$dpar, .data$order) %>%
    tidyr::fill("next_intercept", .direction = "down") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(next_intercept = dplyr::if_else(.data$segment >= .data$next_intercept, NA_integer_, .data$next_intercept))

  # Predictor group-level effects have an independent segment lifetime. A
  # later definition for the same (dpar, grouping factor) replaces the whole
  # current coefficient block; `(0 | group)` is represented as an inactive
  # definition that ends it.
  definitions = lapply(
    seq_along(rhs),
    function(segment) get_predictor_group_definitions_segment(
      rhs[[segment]], segment, family, data, par_x, check_rank
    )
  ) %>%
    dplyr::bind_rows()

  predictor_group_effects = definitions
  if (nrow(definitions) > 0) {
    lifetimes = definitions %>%
      dplyr::distinct(.data$dpar, .data$group_col, .data$segment) %>%
      dplyr::arrange(.data$dpar, .data$group_col, .data$segment) %>%
      dplyr::group_by(.data$dpar, .data$group_col) %>%
      dplyr::mutate(next_segment = as.integer(dplyr::lead(.data$segment))) %>%
      dplyr::ungroup()

    predictor_group_effects = definitions %>%
      dplyr::left_join(
        lifetimes,
        by = c("dpar", "group_col", "segment")
      ) %>%
      dplyr::filter(.data$active) %>%
      dplyr::mutate(
        population_name = dplyr::if_else(
          .data$population_name %in% predictors$code_name,
          .data$population_name,
          NA_character_
        ),
        part = "predictor",
        matrix_col = nrow(predictors) + dplyr::row_number()
      ) %>%
      dplyr::select(
        "population_name", "name", "part", "group_col", "segment", "dpar",
        "sd_name", "par_type", "matrix_name", "display_name", "order",
        "x_factor", "design_id", "design_col", "matrix_col", "matrix_data",
        "next_segment", "correlated", "design_spec"
      )
  }

  # Store each fitted component specification once, not on every coefficient.
  design_specs = collect_design_specs(predictors, definitions)
  predictors = dplyr::select(predictors, -"design_spec")
  if ("design_spec" %in% names(predictor_group_effects))
    predictor_group_effects = dplyr::select(predictor_group_effects, -"design_spec")

  list(
    predictors = predictors,
    group_effects = predictor_group_effects,
    design_specs = design_specs
  )
}


#' @aliases get_predictors
#' @keywords internal
#' @describeIn get_predictors_dpar Return only the population predictor table.
get_predictors = function(model, data, family, par_x, check_rank = TRUE) {
  get_predictor_tables(model, data, family, par_x, check_rank)$predictors
}
