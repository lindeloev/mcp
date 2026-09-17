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
  # Construct model.frame. na.pass prevents silent row omission so we can throw an informative error below
  if (is.null(spec)) {
    frame = stats::model.frame(form, data, na.action = stats::na.pass)
  } else {
    frame = stats::model.frame(spec$terms, data, xlev = spec$xlevels, na.action = stats::na.pass)
  }

  # Transformations can introduce bad values on otherwise OK non-transformed data
  has_na_inf = vapply(frame, function(x) any(if (is.numeric(x)) !is.finite(x) else is.na(x)), logical(1))
  if (any(has_na_inf))
    stop("Predictor transformation resulted in NA or non-finite values: ", and_collapse(names(frame)[has_na_inf]), ".")

  if (is.null(spec)) {
    fitted_terms = attr(frame, "terms")
    matrix = stats::model.matrix(fitted_terms, frame)
    offset = stats::model.offset(frame)
    factor_cols = vapply(frame, function(x) is.factor(x) || is.character(x), logical(1))
    spec = list(
      terms = fitted_terms,
      xlevels = lapply(frame[factor_cols], function(x) if (is.factor(x)) levels(x) else levels(factor(x))),
      contrasts = attr(matrix, "contrasts"),
      columns = colnames(matrix),
      has_offset = !is.null(offset)
    )
  } else {
    matrix = stats::model.matrix(
      spec$terms, frame, contrasts.arg = spec$contrasts
    )
    offset = stats::model.offset(frame)
    if (!identical(colnames(matrix), spec$columns))
      stop("The model matrix for `newdata` does not match the fitted model.")
  }

  # Informative error now instead of downstream in JAGS or prediction
  if (any(!is.finite(matrix)))
    stop("Evaluated design matrix contains non-finite values.")
  if (!is.null(offset) && any(!is.finite(offset)))
    stop("Evaluated offset contains non-finite values.")

  list(matrix = matrix, offset = offset, spec = spec)
}


# Convert `offset(0)` to `offset(0 * par_x)` as a turn-off syntax.
# model.frame() rejects scalar offsets; this ensures newdata works cleanly.
rewrite_zero_offset = function(form, par_x) {
  form_str = deparse1(form)
  rewritten = gsub(
    "\\boffset\\(0\\)", paste0("offset(0 * ", par_x, ")"), form_str
  )
  if (identical(rewritten, form_str))
    return(form)
  stats::as.formula(rewritten, env = environment(form))
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


# Find component in next segment to determine lifetime (exclusive upper boundary)
# - definitions: A data frame with a `segment` column.
# - by: Columns identifying one replaceable component (e.g., c("dpar", "order")).
get_definition_lifetimes = function(definitions, by) {
  if (nrow(definitions) == 0) {
    definitions = definitions[, unique(c(by, "segment")), drop = FALSE]
    definitions$next_segment = integer()
    return(definitions)
  }

  definitions %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(by, "segment")))) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(c(by, "segment")))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::mutate(next_segment = as.integer(dplyr::lead(.data$segment))) %>%
    dplyr::ungroup()
}


#' Rewrite supported segment-local uses of the change-point axis
#'
#' Only bare `x` and polynomial powers expressed as `I(x^k)` are converted to
#' segment-local coordinates. Other functions, transformations, and basis
#' expansions (e.g., `poly(x, 2, raw = TRUE)`) remain on the global scale,
#' altering segment joining rather than simply changing coefficient parameterization.
#'
#' @keywords internal
#' @noRd
#' @param form A one-sided predictor formula.
#' @param par_x Name of the change-point axis.
#' @return The rewritten formula and the local degree of each formula term.
rewrite_local_x = function(form, par_x) {
  local_name = ".mcp_local_x"

  # Recognize only exact x and I(x^k) factors; transformations stay global.
  local_degree = function(x) {
    if (x %in% c(par_x, paste0("I(", par_x, ")")))
      return(1L)
    prefix = paste0("I(", par_x, "^")
    if (!startsWith(x, prefix) || !endsWith(x, ")"))
      return(0L)
    power = substr(x, nchar(prefix) + 1L, nchar(x) - 1L)
    if (grepl("^[+-]?[0-9]+$", power)) as.integer(power) else 0L
  }

  # terms() has expanded formula operators, leaving colon-separated factors.
  rewrite_term = function(x) {
    factors = strsplit(x, ":", fixed = TRUE)[[1]]
    degrees = vapply(factors, local_degree, integer(1))
    is_local = degrees != 0L
    factors[is_local] = sub(par_x, local_name, factors[is_local], fixed = TRUE)
    list(label = paste0(factors, collapse = ":"), degree = sum(degrees))
  }

  # Rewrite the expanded term labels, then restore the original intercept and offset settings.
  terms = stats::terms(form)
  rewritten = lapply(attr(terms, "term.labels"), rewrite_term)
  labels = vapply(rewritten, `[[`, character(1), "label")
  # terms() omits offset() from term.labels; preserve it when reformulating
  if (!is.null(attr(terms, "offset"))) {
    offset_terms = vapply(attr(terms, "offset"), function(i) deparse1(attr(terms, "variables")[[i + 1]]), character(1))
    labels = c(labels, offset_terms)
  }
  form = stats::reformulate(
    labels, intercept = attr(terms, "intercept"), env = environment(form)
  )
  list(
    form = form,
    degree = vapply(rewritten, `[[`, integer(1), "degree"),
    name = local_name
  )
}


# Make a shared R/JAGS-safe coefficient name
make_code_name = function(x) {
  out = gsub("[^A-Za-z0-9_.]", "_", x)
  out = gsub("_+", "_", out)
  needs_prefix = !grepl("^[A-Za-z]", out)
  out[needs_prefix] = paste0("b_", out[needs_prefix])
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

  # Check for exactly one continuous predictor; exclude grouping and offset variables
  rhs_vars = setdiff(get_rhs_non_offset_vars(model), get_rhs_group_vars(model))
  data_in_rhs = data %>% dplyr::select(dplyr::all_of(rhs_vars), dplyr::all_of(par_x))
  continuous_cols = lapply(data_in_rhs, is_continuous) %>% unlist()
  par_x_candidates = names(continuous_cols)[continuous_cols]
  resolved_par_x = if (is.character(par_x)) {
    if (length(par_x_candidates) == 0) {
      par_x
    } else if (par_x %in% par_x_candidates) {
      par_x
    } else {
      stop("Got par_x = '", par_x, "' but it does not seem to be both continuous, present in the data, and in the model. mcp identified '", par_x_candidates, "' as the only viable change point dimension(s) as the data and model is set up now.")
    }
  } else if (is.null(par_x)) {
    if (length(par_x_candidates) == 0) {
      stop("No continuous column for change points found in the formulas. Either provide mcp(..., par_x = 'my_col') or update the model.")
    } else if (length(par_x_candidates) > 1) {
      stop("Could not automatically determine the change point dimension (multiple candidates: ", and_collapse(par_x_candidates), "). Set it explicitly using mcp(..., par_x = 'my_col').")
    } else if (length(par_x_candidates) == 1) {
      par_x_candidates
    }
  } else {
    stop_github("Reached the end of get_par_x() without returning par_x")
  }

  # Warn if par_x is used in offset()
  offset_vars = get_rhs_offset_vars(model)
  if (!is.null(resolved_par_x) && resolved_par_x %in% offset_vars) {
    warning(
      "The change-point predictor '", resolved_par_x, "' is also used in offset().\n",
      "If '", resolved_par_x, "' represents cumulative exposure/duration from time 0, note that mcp applies ",
      "the segment's rate to the entire duration rather than integrating piecewise rates over time. ",
      "Offsets are valid when rows represent discrete intervals with local exposure.",
      call. = FALSE
    )
  }

  return(resolved_par_x)
}


#' Get predictors for one distributional parameter
#'
#' This function extracts a `par_x`-less design matrix.
#' `par_x` will be relative to the segment onset, so it will be multiplied in the formula
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
#' @param check_rank Logical scalar. Whether to stop on rank deficiency.
#' @return A tibble with one row per model parameter and the columns
#'   - `dpar`: character.
#'   - `segment`: the segment number (positive integer).
#'   - `matrix_name`: original column name from the model matrix. Used to
#'     diagnose collisions after parameter names are converted for JAGS.
#'   - `display_name`: user-facing parameter name used in summary functions.
#'   - `code_name`: parameter name used in JAGS and internally in mcp.
#'   - `term_key`: identifier for the formula term which generated the
#'     coefficient. Multi-column terms share one key.
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
      term_key = character(),
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

  form_rhs = rewrite_zero_offset(form_rhs, par_x)

  # Variable names for non-mu terms are prefixed with the term type.
  if (dpar == "mu") {
    dpar_prefix = ""
  } else {
    dpar_prefix = paste0(dpar, order, "_")
  }

  # Disallow multiple terms within functions involving par_x
  formula_terms = attr(stats::terms(form_rhs), "term.labels")
  contains_multiple_terms = formula_terms %>%
    stringr::str_extract("(?<=\\().*(?=\\))") %>%
    stringr::str_detect("[+:*]")
  contains_x = term_contains(par_x, formula_terms)
  is_bad = contains_x & contains_multiple_terms
  if (any(stats::na.omit(is_bad) == TRUE))
    stop("mcp does not currently support 2+ terms within a formula function when one of them is par_x = '", par_x, "'. Found: ", and_collapse(formula_terms[which(is_bad)]))

  # Build the source design on the original data
  source_design = get_fitted_design(form_rhs, data)
  matrix_name = colnames(source_design$matrix)

  # model.matrix() maps every coefficient column to its originating formula
  # term. Multi-column factors and bases therefore share a lifetime key.
  term_assign = attr(source_design$matrix, "assign")
  term_key = c("(Intercept)", formula_terms)[term_assign + 1L]

  # Rewrite formula with a placeholder for the segment-local change-point axis
  local = rewrite_local_x(form_rhs, par_x)
  if (local$name %in% names(data))
    stop("Data column '", local$name, "' is reserved for mcp's formula compiler.")

  # Compile the fitted design specification with segment and offset metadata
  local_data = data
  local_data[[local$name]] = data[[par_x]]
  design = get_fitted_design(local$form, local_data)
  design$spec$local_x_name = local$name
  design$spec$dpar = dpar
  design$spec$segment = segment
  design$spec$order = ifelse(is.null(order), NA_integer_, as.integer(order))
  if (!is.null(source_design$offset)) {
    design$spec$has_offset = TRUE
    design$spec$offset_name = paste0("offset_", dpar, ifelse(is.null(order) || is.na(order), "", order), "_", segment, "_")
    design$spec$offset_data = as.numeric(source_design$offset)
  }

  # Verify that rewriting local_x preserved matrix values and numerical rank
  mat = design$matrix
  if (!identical(dim(mat), dim(source_design$matrix)) ||
      !isTRUE(all.equal(unname(mat), unname(source_design$matrix))))
    stop_github("Rewriting the segment-local change-point axis changed the model matrix.")
  if (check_rank == TRUE)
    assert_rank(source_design$matrix, segment, dpar)


  #######################
  # GET PARAMATER NAMES #
  #######################
  pars = matrix_name

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
  code_name = make_code_name(gsub("[: +]", "", display_name))

  # is_dummy
  is_dummy = apply(mat, 2, function(x) all(x %in% c(0, 1)))


  ################
  # GET X_FACTOR #
  ################

  # Bare par_x and supported powers are relative to the segment onset. The
  # model-matrix assign vector maps expanded factor columns back to terms.
  local_degree = c(0L, local$degree)[attr(mat, "assign") + 1L]
  checkmate::assert_integerish(local_degree, lower = 0, .var.name = "exponents in formula")
  x_factor = ifelse(
    local_degree == 0L, "1",
    ifelse(local_degree == 1L, "x", paste0("x^", local_degree))
  )



  ################################
  # COMPUTE X-LESS DESIGN MATRIX #
  ################################

  # Evaluate the same fitted design with the local factor set to one. This is
  # exact at zeros and lets model.matrix handle interactions and contrasts.
  local_data[[local$name]] = 1
  mat_without_x = get_fitted_design(data = local_data, spec = design$spec)$matrix

  if (ncol(mat_without_x) == 0) {
    if (isTRUE(design$spec$has_offset)) {
      # Offset-only segment with no estimated coefficients (e.g., ~ 0 + offset(z));
      # return a placeholder row to carry design_spec into model_tables$design_specs
      return(tibble::tibble(
        dpar = dpar,
        segment = segment,
        matrix_name = NA_character_,
        display_name = NA_character_,
        code_name = NA_character_,
        term_key = "(offset)",
        par_type = "offset",
        order = ifelse(is.null(order), NA_integer_, as.integer(order)),
        x_factor = "1",
        design_id = design_id,
        design_col = NA_integer_,
        matrix_data = list(numeric(nrow(data))),
        design_spec = list(design$spec)
      ))
    } else {
      return(tibble::tibble(
        dpar = character(),
        segment = integer(),
        matrix_name = character(),
        display_name = character(),
        code_name = character(),
        term_key = character(),
        par_type = character(),
        order = integer(),
        x_factor = character(),
        design_id = character(),
        design_col = integer(),
        design_spec = list(),
        matrix_data = list()
      ))
    }
  }

  predictors = tibble::tibble(
    dpar = dpar,
    segment = segment,
    matrix_name = matrix_name,
    display_name,
    code_name = code_name,
    term_key = term_key,
    par_type = dplyr::case_when(
      is_intercept == TRUE ~ "Intercept",
      is_dummy == TRUE ~ "dummy",
      TRUE ~ "slope"
    ),
    order = ifelse(is.null(order), NA, order),
    x_factor = x_factor,
    design_id = design_id,
    design_col = seq_len(ncol(mat_without_x)),
    design_spec = rep(list(design$spec), ncol(mat_without_x)),
    matrix_data = lapply(seq_len(ncol(mat_without_x)), function(i) mat_without_x[, i])
  )

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
  collision_names = unique(predictors$code_name[duplicated(predictors$code_name)])
  if (length(collision_names) == 0)
    return(invisible(predictors))

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


# Detect terms that contain a particular variable
term_contains = function(par_x, terms) {
  vapply(terms, function(term) par_x %in% all.vars(str2lang(term)), logical(1))
}


# Normalize coefficient-sharing selectors before ordinary formula parsing. The
# parser below only needs to know which expanded source terms a selector names;
# all downstream code receives resolved definition/occurrence rows.

# Check whether an expression is a call to same()
is_same_call = function(expr) {
  is.call(expr) && identical(deparse1(expr[[1]]), "same")
}


# Recursively check whether an expression contains any call to same()
contains_same_call = function(expr) {
  if (is_same_call(expr))
    return(TRUE)
  if (!is.call(expr))
    return(FALSE)

  any(vapply(as.list(expr)[-1], contains_same_call, logical(1)))
}


# Extract and validate formula terms selected by a same() call
same_selector_terms = function(selector, env) {
  # Disallow unsupported calls inside same()
  if (is.call(selector) && deparse1(selector[[1]]) %in% c("offset", "stats::offset"))
    stop("`same(offset(...))` is not supported. Repeat `offset(...)` in each active segment.", call. = FALSE)
  if (is.call(selector) && deparse1(selector[[1]]) %in% c("ar", "ma", known_dpar_wrappers()))
    stop("Select terms inside their component, e.g. `sigma(same(1))`, rather than wrapping a component in `same(sigma(1))`.", call. = FALSE)
  if (contains_same_call(selector))
    stop("Nested `same()` calls are not supported.", call. = FALSE)

  # Convert selector to a formula and extract terms
  selector_form = stats::as.formula(call("~", call("+", 0, selector)), env = env)
  attrs = attributes(stats::terms(selector_form))
  terms = attrs$term.labels

  # Validate that group terms and offsets are not included
  if (any(vapply(terms, is_group_term, logical(1))) || grepl("\\|", deparse1(selector)))
    stop("Sharing group-effect blocks is not supported yet. It will be added in a later mcp v0.4 step.", call. = FALSE)
  if (!is.null(attrs$offset))
    stop("`same(offset(...))` is not supported. Repeat `offset(...)` in each active segment.", call. = FALSE)

  # Include intercept if present and ensure at least one term was selected
  terms = c(if (attrs$intercept == 1) "(Intercept)", terms)
  if (length(terms) == 0)
    stop("`same()` must select at least one coefficient.", call. = FALSE)

  terms
}


# Parse and validate a same() call into a selector specification
parse_same_call = function(expr, segment, env) {
  # Check argument count and naming
  args = as.list(expr)[-1]
  arg_names = names(args)
  if (is.null(arg_names))
    arg_names = rep("", length(args))
  if (length(args) == 0 || length(args) > 2 || arg_names[1] != "" ||
      (length(args) == 2 && arg_names[2] != "as"))
    stop("`same()` accepts one selection and an optional named `as` segment, e.g. `same(x, as = 1)`.", call. = FALSE)
  if (segment == 1)
    stop("`same()` cannot be used in segment 1 because there is no earlier coefficient to share.", call. = FALSE)

  # Resolve and validate the source segment
  source = segment - 1L
  if (length(args) == 2) {
    source_value = args[[2]]
    if (!is.numeric(source_value) || length(source_value) != 1 ||
        !is.finite(source_value) || source_value <= 0 || source_value != floor(source_value))
      stop("`same(..., as = )` must name an earlier positive integer segment.", call. = FALSE)
    source = source_value
  }
  if (source >= segment)
    stop("`same(..., as = )` must name an earlier segment.", call. = FALSE)

  # Return parsed selector specification
  tibble::tibble(
    term_key = same_selector_terms(args[[1]], env),
    source = as.integer(source)
  )
}


# Build a component design and replace selected terms with fitted source columns
get_shared_predictors = function(data, form_rhs, segment, dpar, par_x, order = NULL,
                                  check_rank = TRUE, design_id = NULL, previous = NULL) {
  # Expand additive selectors once, retaining ordinary formula coding context
  leaves = unpack_additive(form_rhs[[2]])
  shared = vapply(leaves, is_same_call, logical(1))
  if (any(vapply(leaves[!shared], contains_same_call, logical(1))))
    stop("`same()` must be an additive selector of a complete term, such as `same(x:z)`.", call. = FALSE)
  selected = dplyr::bind_rows(lapply(leaves[shared], parse_same_call,
    segment = segment, env = environment(form_rhs)))
  if (nrow(selected) > 0) {
    selected = dplyr::distinct(selected)
    if (anyDuplicated(selected$term_key))
      stop("Overlapping `same()` selections cannot use competing sources.", call. = FALSE)

    # Shared intercepts replace implicit intercepts, but conflict with explicit 1
    bare = leaves[!shared]
    shared_intercept = "(Intercept)" %in% selected$term_key
    if (shared_intercept && any(vapply(bare, function(x) identical(x, 1) || identical(x, 1L), logical(1))))
      stop("An explicit bare `1` and `same(1)` are competing declarations.", call. = FALSE)
    bare_expr = Reduce(function(a, b) call("+", a, b), c(list(0), bare))
    bare_terms = attr(stats::terms(stats::as.formula(call("~", bare_expr), env = environment(form_rhs))), "term.labels")
    if (any(selected$term_key %in% bare_terms))
      stop("A term cannot be both bare and shared in segment ", segment, ".", call. = FALSE)

    # Unwrap selectors for coding; 0 + same(1) still declares an intercept
    leaves[shared] = lapply(leaves[shared], function(x) x[[2]])
    expr = Reduce(function(a, b) call("+", a, b), leaves)
    if (shared_intercept)
      expr = call("+", expr, 1)
    form_rhs = stats::as.formula(call("~", expr), env = environment(form_rhs))
  }

  # Construct bare columns in the full unwrapped formula, then borrow source columns
  form_rhs = remove_terms(form_rhs, "varying")
  current = get_predictors_dpar(data, form_rhs, segment, dpar, par_x, order,
    check_rank = check_rank && nrow(selected) == 0, design_id = design_id) %>%
    dplyr::mutate(definition_name = .data$code_name, definition_segment = .data$segment)
  if (nrow(selected) == 0)
    return(current)

  borrowed = lapply(seq_len(nrow(selected)), function(i) {
    source = previous[previous$dpar == dpar & previous$segment == selected$source[i] &
      previous$order %in% (if (is.null(order)) NA_integer_ else order) &
      previous$term_key == selected$term_key[i], , drop = FALSE]
    if (is.null(source) || nrow(source) == 0)
      stop("`same()` cannot find `", selected$term_key[i], "` in segment ",
        selected$source[i], " for ", dpar, if (!is.null(order)) paste0(" lag ", order), ".", call. = FALSE)
    source$segment = segment
    source
  })
  bare = current[!current$term_key %in% selected$term_key, , drop = FALSE]

  # Keep destination offsets even when all its coefficients are shared
  if (nrow(bare) == 0 && nrow(current) > 0 && isTRUE(current$design_spec[[1]]$has_offset)) {
    bare = current[1, ]
    bare$par_type = "offset"
    bare$code_name = NA_character_
  }
  result = dplyr::bind_rows(bare, dplyr::bind_rows(borrowed))

  # Rank checking applies to the assembled design, including borrowed contrasts
  if (check_rank) {
    coefficients = result[result$par_type != "offset", , drop = FALSE]
    matrix = do.call(cbind, lapply(seq_len(nrow(coefficients)), function(i) {
      degree = if (coefficients$x_factor[i] == "1") 0 else
        if (coefficients$x_factor[i] == "x") 1 else as.numeric(sub("x^", "", coefficients$x_factor[i], fixed = TRUE))
      coefficients$matrix_data[[i]] * data[[par_x]]^degree
    }))
    colnames(matrix) = coefficients$matrix_name
    assert_rank(matrix, segment, dpar)
  }
  result
}


# Does an explicit dpar formula provide no initial predictor at all? This is
# deliberately syntactic: a nonzero offset or a group-only formula counts as
# a declaration even though it has no population coefficient.
is_empty_initial_predictor = function(form) {
  attrs = attributes(stats::terms(form))
  group_terms = attrs$term.labels[vapply(attrs$term.labels, is_group_term, logical(1))]
  population_terms = setdiff(attrs$term.labels, group_terms)
  group_is_active = vapply(group_terms, function(term) {
    coefficient_form = stats::as.formula(call("~", str2lang(term)[[2]]), env = environment(form))
    coefficient_attrs = attributes(stats::terms(coefficient_form))
    coefficient_attrs$intercept == 1 || length(coefficient_attrs$term.labels) > 0
  }, logical(1))

  # Check regular coefficient-terms; if no issue is found return TRUE
  if (attrs$intercept == 1 || length(population_terms) > 0 || any(group_is_active))
    return(FALSE)
  if (is.null(attrs$offset))
    return(TRUE)

  # Check offset terms
  offset_terms = vapply(
    attrs$offset,
    function(i) deparse1(attrs$variables[[i + 1]]),
    character(1)
  )
  all(offset_terms %in% c("offset(0)", "stats::offset(0)"))
}


#' @aliases get_predictors_segment
#' @keywords internal
#' @noRd
#' @describeIn get_predictors_dpar Apply `get_predictors_dpar` to
#'   each formula in a segment
get_predictors_segment = function(form_rhs, segment, family, data, par_x, check_rank = TRUE, previous = NULL) {
  checkmate::assert_formula(form_rhs)
  checkmate::assert_int(segment, lower = 1)
  checkmate::assert_true(is.mcpfamily(family), .var.name = "family")
  checkmate::assert_data_frame(data)
  checkmate::assert_string(par_x)

  # Components have already been canonicalized by get_predictor_tables()
  form_env = environment(form_rhs)
  attrs = attributes(stats::terms(form_rhs))
  term_labels = attrs$term.labels

  # Formula wrappers belonging to distributional parameters and AR/MA
  model_dpars = family$dpar_specs$dpar
  arma_components = c("ar", "ma")

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
      dpar_form = stats::as.formula("~1", env = form_env)
      dpar_pars[[dpar]] = get_shared_predictors(
        data, dpar_form, segment, dpar = dpar, par_x, NULL, check_rank,
        design_id = paste("population", dpar, segment, sep = ":"), previous = previous
      ) %>%
        dplyr::mutate(explicit = FALSE)
    } else if (length(dpar_term) > 0) {
      dpar_form = get_term_content(dpar_term, form_env)
      if (segment == 1 && spec$require_initial_predictor && is_empty_initial_predictor(dpar_form))
        stop(
          "`", dpar, "(0)` cannot be the initial predictor for family = ", family$family,
          "(). Declare an initial predictor, such as `", dpar, "(1)` or `", dpar, "(0 + x)`.",
          call. = FALSE
        )
      dpar_pars[[dpar]] = get_shared_predictors(
        data, dpar_form, segment, dpar = dpar, par_x, NULL, check_rank,
        design_id = paste("population", dpar, segment, sep = ":"), previous = previous
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
      component_form = get_term_content(component_stuff$form_str, form_env)
      if (length(get_group_terms(component_form)) > 0)
        stop(
          "Group-level effects inside ", component,
          "() are not currently supported. Found one in segment ", segment, "."
        )
      # Expand one formula into a separate regression parameter for each lag.
      arma_pars[[component]] = if (component_stuff$order == 0) {
        get_shared_predictors(
          data, component_form, segment, component, par_x, order = 1L, check_rank,
          design_id = paste("population", component, 0, segment, sep = ":"), previous = previous
        ) %>%
          dplyr::mutate(boundary = component_stuff$boundary, explicit = TRUE)
      } else {
        lapply(
          seq_len(component_stuff$order),
          function(order) get_shared_predictors(
            data, component_form, segment, component, par_x, order, check_rank,
            design_id = paste("population", component, order, segment, sep = ":"), previous = previous
          )
        ) %>%
          dplyr::bind_rows() %>%
          dplyr::mutate(boundary = component_stuff$boundary, explicit = TRUE)
      }
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
    dplyr::bind_rows(dpar_pars),
    dplyr::bind_rows(arma_pars)
  )
}


#' @aliases get_predictor_tables
#' @keywords internal
#' @describeIn get_predictors_dpar Apply `get_predictors_segment`
#'   to all segments of a model.
get_predictor_tables = function(model, data, family, par_x, check_rank = TRUE) {
  # Canonicalize segment formulas to wrap bare mu terms into mu(...)
  rhs = lapply(model, get_rhs)
  rhs = lapply(rhs, canonicalize_rhs, family = family)

  # Resolve each component against earlier active occurrences, including chains
  parsed_predictors = NULL
  for (segment in seq_along(rhs)) {
    parsed_predictors = dplyr::bind_rows(parsed_predictors,
      get_predictors_segment(rhs[[segment]], segment, family, data, par_x,
        check_rank, previous = parsed_predictors))
  }
  parsed_predictors = dplyr::arrange(parsed_predictors, .data$dpar, .data$segment)
  predictors = parsed_predictors %>%
    dplyr::filter(.data$par_type != "offset") %>%
    dplyr::mutate(matrix_col = dplyr::row_number())

  # Only the defining occurrence introduces a parameter and prior
  predictor_definitions = predictors %>%
    dplyr::filter(.data$segment == .data$definition_segment) %>%
    dplyr::select(-"design_spec")
  assert_unique_predictor_names(predictor_definitions)
  if ("boundary" %notin% names(predictors))
    predictors$boundary = rep(NA_real_, nrow(predictors))

  # Population intercepts reset the current population predictor. First find
  # the next intercept for each distributional parameter and AR/MA order.
  # Strategy: (1) select one row for segments with intercepts for each dpar (filter)
  #           (2) save this segment number in the last segment that had an intercept (lag)
  #           (3) left-join this into predictors
  #           (4) fill downwards into intermittent segments without intercepts
  df_next_intercept = predictors %>%
    dplyr::filter(.data$par_type == "Intercept") %>%
    get_definition_lifetimes(c("dpar", "order")) %>%
    dplyr::rename(next_intercept = "next_segment")

  # Population predictors: left-join and fill-down. NA means "there is no next
  # intercept segment".
  predictors = predictors %>%
    dplyr::left_join(df_next_intercept, by = c("dpar", "segment", "order")) %>%
    dplyr::group_by(.data$dpar, .data$order) %>%
    tidyr::fill("next_intercept", .direction = "down") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(next_intercept = dplyr::if_else(.data$segment >= .data$next_intercept, NA_integer_, .data$next_intercept))

  # Population non-local terms are active only in the segment where they are
  # declared. Local par_x terms retain their endpoint through joined segments
  # and therefore still end at the next population intercept.
  predictors = predictors %>%
    dplyr::mutate(
      next_segment = dplyr::if_else(
        .data$par_type != "Intercept" & .data$x_factor == "1" &
          .data$dpar %notin% c("ar", "ma"),
        dplyr::if_else(.data$segment < length(rhs), .data$segment + 1L, NA_integer_),
        .data$next_intercept
      )
    ) %>%
    dplyr::select(-"next_intercept")

  # AR/MA retain endpoints through declared, joined lags until an intercept reset
  arma_declarations = get_arma_declarations(rhs)
  arma_next_segment = vapply(seq_len(nrow(predictors)), function(i) {
    if (predictors$dpar[i] %notin% c("ar", "ma"))
      return(NA_integer_)

    next_segment = predictors$segment[i] + 1L
    joins = predictors$par_type[i] == "Intercept" || predictors$x_factor[i] != "1"
    while (joins && next_segment <= length(rhs)) {
      declaration = arma_declarations[
        arma_declarations$dpar == predictors$dpar[i] &
          arma_declarations$segment == next_segment,
        , drop = FALSE
      ]
      resets = any(predictors$dpar == predictors$dpar[i] &
        predictors$order %in% predictors$order[i] & predictors$segment == next_segment &
        predictors$par_type == "Intercept")
      if (nrow(declaration) != 1 || declaration$order < predictors$order[i] || resets)
        break
      next_segment = next_segment + 1L
    }
    if (next_segment > length(rhs)) NA_integer_ else next_segment
  }, integer(1))
  predictors = predictors %>%
    dplyr::mutate(
      next_segment = dplyr::if_else(
        .data$dpar %in% c("ar", "ma"),
        arma_next_segment,
        .data$next_segment
      )
    )

  # Group intercepts and non-local terms are active only where declared.
  # Local group-x terms retain their endpoint through joined blocks until a
  # later group intercept or explicit `(0 | group)` resets the block.
  definitions = lapply(
    seq_along(rhs),
    function(segment) get_predictor_group_definitions_segment(
      rhs[[segment]], segment, family, data, par_x, check_rank
    )
  ) %>%
    dplyr::bind_rows()

  predictor_group_effects = definitions
  if (nrow(definitions) > 0) {
    predictor_group_effects = definitions %>%
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
        "correlated", "design_spec"
      )

    resets = definitions %>%
      dplyr::group_by(.data$dpar, .data$group_col, .data$segment) %>%
      dplyr::summarise(
        reset = any(!.data$active | .data$par_type == "Intercept"),
        .groups = "drop"
      ) %>%
      dplyr::filter(.data$reset)

    next_reset = vapply(seq_len(nrow(predictor_group_effects)), function(i) {
      candidates = resets$segment[
        resets$dpar == predictor_group_effects$dpar[i] &
          resets$group_col == predictor_group_effects$group_col[i] &
          resets$segment > predictor_group_effects$segment[i]
      ]
      if (length(candidates) == 0) NA_integer_ else min(candidates)
    }, integer(1))

    predictor_group_effects = predictor_group_effects %>%
      dplyr::mutate(next_segment = dplyr::if_else(
        .data$x_factor == "1",
        dplyr::if_else(.data$segment < length(rhs), .data$segment + 1L, NA_integer_),
        next_reset
      ))
  }

  # Store each fitted component specification once, not on every coefficient.
  design_specs = collect_design_specs(parsed_predictors, definitions)
  predictors = predictors %>%
    dplyr::select(-"design_spec")
  if ("design_spec" %in% names(predictor_group_effects))
    predictor_group_effects = dplyr::select(predictor_group_effects, -"design_spec")
  if ("name" %notin% names(predictor_group_effects))
    predictor_group_effects = tibble::tibble(name = character(), segment = integer())

  # Attach occurrence segment and definition metadata
  predictors = predictors %>%
    dplyr::mutate(occurrence_segment = .data$segment)
  group_definitions = predictor_group_effects %>%
    dplyr::mutate(definition_name = .data$name, definition_segment = .data$segment)
  predictor_group_effects = predictor_group_effects %>%
    dplyr::mutate(definition_name = .data$name, definition_segment = .data$segment, occurrence_segment = .data$segment)

  list(
    predictors = predictors,
    group_effects = predictor_group_effects,
    predictor_definitions = predictor_definitions,
    group_definitions = group_definitions,
    design_specs = design_specs
  )
}


#' @aliases get_predictors
#' @keywords internal
#' @describeIn get_predictors_dpar Return only the population predictor table.
get_predictors = function(model, data, family, par_x, check_rank = TRUE) {
  get_predictor_tables(model, data, family, par_x, check_rank)$predictors
}
