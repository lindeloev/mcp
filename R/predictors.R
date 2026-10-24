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
#' @param form A one-sided predictor formula.
#' @param data A data frame used to fit the design.
#' @param spec An optional fitted specification returned by this function.
#' @return A list with `matrix`, `offset` (`NULL` without `offset()`), and `spec`.
#' @noRd
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
#' @param ... Predictor tables that may contain `design_id` and `design_spec`.
#' @return A named list of fitted design specifications.
#' @noRd
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


# First segment where each term no longer applies (`next_segment`; NA = through the last segment).
# Joined terms apply until the next segment that resets their `key`; other terms apply in their own segment only.
# - key, segment, joins: One value per term.
# - resets: A data frame with `key` and `segment` columns.
get_next_segment = function(key, segment, joins, resets, n_segments) {
  vapply(seq_along(segment), function(i) {
    later = if (joins[i]) resets$segment[resets$key == key[i] & resets$segment > segment[i]] else segment[i] + 1L
    later = later[later <= n_segments]
    if (length(later) == 0) NA_integer_ else as.integer(min(later))
  }, integer(1))
}


#' Rewrite x-terms to be measured from the change point
#'
#' Only bare `x` and polynomial powers expressed as `I(x^k)` are x-terms, which
#' are measured from the change point. Other functions, transformations, and basis
#' expansions (e.g., `poly(x, 2, raw = TRUE)`) use the original values of `x`
#' and do not join.
#'
#' @param form A one-sided predictor formula.
#' @param par_x Name of the change-point axis.
#' @return The rewritten formula and the local degree of each formula term.
#' @noRd
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
make_parameter_name = function(x) {
  out = gsub("[^A-Za-z0-9_.]", "_", x)
  out = gsub("_+", "_", out)
  needs_prefix = !grepl("^[A-Za-z]", out)
  out[needs_prefix] = paste0("b_", out[needs_prefix])
  out
}

#' Detects par_x and verifies model-data fit
#'
#' @inheritParams mcp
#' @param vars Model variables by role, returned by `get_formula_vars()`.
#' @return The column name of par_x.
#' @noRd
get_par_x = function(vars, data, par_x = NULL) {
  checkmate::assert_data_frame(data)
  checkmate::assert_string(par_x, null.ok = TRUE)

  if (!is.null(par_x)) {
    # Check user-supplied par_x
    if (par_x %notin% colnames(data))
      stop("par_x = '", par_x, "' not found in data.")
    if (inherits(data[[par_x]], c("Date", "POSIXt")))
      stop("par_x = '", par_x, "' is a date. Convert it to a number first, e.g., days since the first date: `as.numeric(difftime(data$", par_x, ", min(data$", par_x, "), units = \"days\"))`.")
    if (!is_continuous(data[[par_x]]))
      stop("par_x = '", par_x, "' has to be continuous. Is it binary or categorical?")
  } else {
    # Detect exactly one continuous predictor; exclude grouping, offset, and series variables
    rhs_vars = setdiff(vars$predictor, vars$group)
    continuous_cols = vapply(dplyr::select(data, dplyr::all_of(rhs_vars)), is_continuous, logical(1))
    par_x_candidates = names(continuous_cols)[continuous_cols]
    if (length(par_x_candidates) == 0)
      stop("No continuous column for change points found in the formulas. Either provide mcp(..., par_x = 'my_col') or update the model.")
    if (length(par_x_candidates) > 1)
      stop("Could not automatically determine the change-point variable (multiple candidates: ", and_collapse(par_x_candidates), "). Set it explicitly using mcp(..., par_x = 'my_col').")
    par_x = par_x_candidates
  }

  # Warn if par_x is used in offset()
  if (par_x %in% vars$offset) {
    warning(
      "The change-point predictor '", par_x, "' is also used in offset().\n",
      "If '", par_x, "' represents cumulative exposure/duration from time 0, note that mcp applies ",
      "the segment's rate to the entire duration rather than integrating piecewise rates over time. ",
      "Offsets are valid when rows represent discrete intervals with local exposure.",
      call. = FALSE
    )
  }

  par_x
}


# Zero-row predictor table with the columns returned by get_predictors_dpar()
empty_predictors = function() {
  tibble::tibble(
    dpar = character(),
    segment = integer(),
    matrix_name = character(),
    name = character(),
    term_key = character(),
    par_type = character(),
    order = integer(),
    x_degree = integer(),
    design_id = character(),
    design_col = integer(),
    design_spec = list(),
    matrix_data = list()
  )
}


#' Get predictors for one distributional parameter
#'
#' This function extracts a `par_x`-less design matrix.
#' x-terms are measured from the change point, so `par_x` is multiplied in the formula
#' (`jags_code` and `fit$simulate()`).
#'
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
#'   - `name`: parameter name, e.g., `x_2` or `x.stateB_1`, used in priors, JAGS, and draws.
#'   - `term_key`: identifier for the formula term which generated the
#'     coefficient. Multi-column terms share one key.
#'   - `par_type`: One of "Intercept", "dummy", or "slope". Used for setting priors and for change point indicator func.
#'   - `order`: positive integer or NA. Only relevant for `ar` and `ma`.
#'   - `x_degree`: The power of x (measured from the change point) that multiplies the coefficient, e.g., 0 for none and 2 for `x^2`.
#'   - `design_id`: key of the fitted component formula that produced the row.
#'   - `design_col`: column occupied by the row in that component's model matrix.
#'   - `matrix_data`: column of the design matrix less the `par_x` term.
#' @noRd
get_predictors_dpar = function(data, form_rhs, segment, dpar, par_x, order = NULL,
                               check_rank = TRUE, design_id = NULL) {
  # Empty segments return no rows
  if (all(as.character(form_rhs) == c("~", "0")))
    return(empty_predictors())

  checkmate::assert_data_frame(data)
  checkmate::assert_formula(form_rhs)
  checkmate::assert_int(length(form_rhs), lower = 2, upper = 2, .var.name = "length(form_rhs)")
  checkmate::assert_int(segment, lower = 1)
  checkmate::assert_string(dpar)
  checkmate::assert_string(par_x)
  checkmate::assert_string(design_id)
  checkmate::assert_int(order, lower = 1, null.ok = TRUE)
  order_value = if (is.null(order)) NA_integer_ else as.integer(order)  # Stored in the returned table

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

  # Rewrite formula with a placeholder for x measured from the change point
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
  design$spec$order = order_value
  if (!is.null(source_design$offset)) {
    design$spec$has_offset = TRUE
    design$spec$offset_name = paste0("offset_", dpar, order, "_", segment, "_")
    design$spec$offset_data = as.numeric(source_design$offset)
  }

  # Verify that rewriting local_x preserved matrix values and numerical rank
  mat = design$matrix
  if (!identical(dim(mat), dim(source_design$matrix)) ||
      !isTRUE(all.equal(unname(mat), unname(source_design$matrix))))
    stop_github("Rewriting x-terms to be measured from the change point changed the model matrix.")
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
  intercept_name = ifelse(dpar %in% c("ar", "ma"), "", "Intercept")
  pars[is_intercept] = intercept_name

  # Parameter names, e.g., xE2_1 for I(x^2) and x.stateB_1 for x:stateB. Symbols are renamed as in brms,
  # except that interactions use "." because JAGS names cannot contain ":".
  symbols = c("^" = "E", "-" = "M", "+" = "P", "*" = "MU", "/" = "D", ":" = ".")
  name = gsub("\\(|\\)", "", pars)
  for (symbol in names(symbols))
    name = gsub(symbol, symbols[[symbol]], name, fixed = TRUE)
  name = paste0(dpar_prefix, name, "_", segment)

  # Multi-column bases often include arguments and namespaces in their matrix
  # names. Name both the position within the basis and the segment explicitly.
  is_basis = grepl("[,=]|::", matrix_name)
  term_index = attr(mat, "assign")[is_basis]
  basis_name = formula_terms[term_index]
  basis_name = sub("^.*::", "", basis_name)
  basis_name = gsub("[^A-Za-z0-9]+", "_", basis_name)
  basis_name = gsub("^_|_$", "", basis_name)
  basis_col = ave(term_index, term_index, FUN = seq_along)
  name[is_basis] = paste0(
    dpar_prefix, basis_name, "_basis", basis_col, "_", segment
  )
  name = make_parameter_name(gsub(" ", "", name, fixed = TRUE))

  # is_dummy
  is_dummy = apply(mat, 2, function(x) all(x %in% c(0, 1)))


  ################
  # GET X_DEGREE #
  ################

  # x-terms (bare par_x and supported powers) are measured from the change point. The
  # model-matrix assign vector maps expanded factor columns back to terms.
  x_degree = c(0L, local$degree)[attr(mat, "assign") + 1L]
  checkmate::assert_integerish(x_degree, lower = 0, .var.name = "exponents in formula")



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
        name = NA_character_,
        term_key = "(offset)",
        par_type = "offset",
        order = order_value,
        x_degree = 0L,
        design_id = design_id,
        design_col = NA_integer_,
        matrix_data = list(numeric(nrow(data))),
        design_spec = list(design$spec)
      ))
    } else {
      return(empty_predictors())
    }
  }

  predictors = tibble::tibble(
    dpar = dpar,
    segment = segment,
    matrix_name = matrix_name,
    name = name,
    term_key = term_key,
    par_type = dplyr::case_when(
      is_intercept == TRUE ~ "Intercept",
      is_dummy == TRUE ~ "dummy",
      TRUE ~ "slope"
    ),
    order = order_value,
    x_degree = x_degree,
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
#' @param predictors A data frame returned by `get_predictors_dpar()`.
#' @return `predictors`, invisibly. Stops with an informative error on
#'   collision.
#' @noRd
assert_unique_predictor_names = function(predictors) {
  collision_names = unique(predictors$name[duplicated(predictors$name)])
  if (length(collision_names) == 0)
    return(invisible(predictors))

  collision_lines = vapply(collision_names, function(name) {
    rows = predictors[predictors$name == name, , drop = FALSE]
    sources = paste0(
      "`", rows$matrix_name, "` (", rows$dpar,
      ", segment ", rows$segment, ")"
    )
    paste0("  `", name, "`: ", and_collapse(sources))
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


# Standardize coefficient-sharing selectors before ordinary formula parsing. The
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
    stop("`same(offset(...))` is not supported. Repeat `offset(...)` in each segment where it applies.", call. = FALSE)
  if (is.call(selector) && deparse1(selector[[1]]) %in% c("ar", "ma", known_dpar_wrappers()))
    stop("Select terms inside their component, e.g. `sigma(same(1))`, rather than wrapping a component in `same(sigma(1))`.", call. = FALSE)
  if (contains_same_call(selector))
    stop("Nested `same()` calls are not supported.", call. = FALSE)

  # Convert selector to a formula and extract terms
  selector_form = stats::as.formula(call("~", call("+", 0, selector)), env = env)
  attrs = attributes(stats::terms(selector_form))
  terms = attrs$term.labels

  # Offsets have no coefficient to share
  if (!is.null(attrs$offset))
    stop("`same(offset(...))` is not supported. Repeat `offset(...)` in each segment where it applies.", call. = FALSE)

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
    stop("`same()` cannot be used in segment 1 because there is no earlier coefficient to reuse.", call. = FALSE)

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


# Standardize selectors once for both population terms and complete group blocks
standardize_shared_component = function(form_rhs, segment) {
  # Expand additive selectors once, retaining ordinary formula coding context
  leaves = unpack_additive(form_rhs[[2]])
  shared = vapply(leaves, is_same_call, logical(1))
  if (any(vapply(leaves[!shared], contains_same_call, logical(1))))
    stop("`same()` must be an additive selector of a complete term, such as `same(x:z)`.", call. = FALSE)
  selected = dplyr::bind_rows(lapply(leaves[shared], parse_same_call,
    segment = segment, env = environment(form_rhs)))
  if (nrow(selected) == 0)
    selected = tibble::tibble(term_key = character(), source = integer())
  if (nrow(selected) > 0) {
    selected = dplyr::distinct(selected)
    if (anyDuplicated(selected$term_key))
      stop("Overlapping `same()` selections cannot use competing sources.", call. = FALSE)

    # Shared intercepts replace implicit intercepts, but conflict with explicit 1
    bare = leaves[!shared]
    shared_intercept = "(Intercept)" %in% selected$term_key
    if (shared_intercept && any(vapply(bare, function(x) identical(x, 1) || identical(x, 1L), logical(1))))
      stop("An explicit bare `1` and `same(1)` are competing intercepts; include only one.", call. = FALSE)
    bare_expr = Reduce(function(a, b) call("+", a, b), c(list(0), bare))
    bare_terms = attr(stats::terms(stats::as.formula(call("~", bare_expr), env = environment(form_rhs))), "term.labels")
    if (any(selected$term_key %in% bare_terms))
      stop("A term cannot be both included directly and reused with `same()` in segment ", segment, ".", call. = FALSE)

    # Unwrap selectors for coding; 0 + same(1) still includes an intercept
    leaves[shared] = lapply(leaves[shared], function(x) x[[2]])
    expr = Reduce(function(a, b) call("+", a, b), leaves)
    if (shared_intercept)
      expr = call("+", expr, 1)
    form_rhs = stats::as.formula(call("~", expr), env = environment(form_rhs))
  }
  list(form = form_rhs, selected = selected)
}


# Build a component design and replace selected terms with fitted source columns
get_shared_predictors = function(data, form_rhs, segment, dpar, par_x, order = NULL,
                                  check_rank = TRUE, design_id = NULL, previous = NULL) {
  standardized = standardize_shared_component(form_rhs, segment)
  form_rhs = standardized$form
  selected = standardized$selected
  selected = selected[!vapply(selected$term_key, is_group_term, logical(1)), , drop = FALSE]

  # Construct bare columns in the full unwrapped formula, then borrow source columns
  form_rhs = remove_terms(form_rhs, "group")
  current = get_predictors_dpar(data, form_rhs, segment, dpar, par_x, order,
    check_rank = check_rank && nrow(selected) == 0, design_id = design_id) %>%
    dplyr::mutate(definition_segment = .data$segment)
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

  # Keep current segment offsets even when all its coefficients are shared
  if (nrow(bare) == 0 && nrow(current) > 0 && isTRUE(current$design_spec[[1]]$has_offset)) {
    bare = current[1, ]
    bare$par_type = "offset"
    bare$name = NA_character_
  }
  result = dplyr::bind_rows(bare, dplyr::bind_rows(borrowed))

  # Rank checking applies to the assembled design, including borrowed contrasts
  if (check_rank) {
    coefficients = result[result$par_type != "offset", , drop = FALSE]
    matrix = do.call(cbind, lapply(seq_len(nrow(coefficients)), function(i) {
      coefficients$matrix_data[[i]] * data[[par_x]]^coefficients$x_degree[i]
    }))
    colnames(matrix) = coefficients$matrix_name
    assert_rank(matrix, segment, dpar)
  }
  result
}


# Does an explicit dpar formula provide no initial predictor at all? This is
# deliberately syntactic: a nonzero offset or a group-only formula counts as
# an initial predictor even though it has no population coefficient.
is_empty_initial_predictor = function(form) {
  attrs = attributes(stats::terms(form))
  group_terms = attrs$term.labels[vapply(attrs$term.labels, is_group_term, logical(1))]
  population_terms = setdiff(attrs$term.labels, group_terms)
  group_has_coefficients = vapply(group_terms, function(term) {
    coefficient_form = stats::as.formula(call("~", str2lang(term)[[2]]), env = environment(form))
    coefficient_attrs = attributes(stats::terms(coefficient_form))
    coefficient_attrs$intercept == 1 || length(coefficient_attrs$term.labels) > 0
  }, logical(1))

  # Check regular coefficient-terms; if no issue is found return TRUE
  if (attrs$intercept == 1 || length(population_terms) > 0 || any(group_has_coefficients))
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


#' Apply `get_predictors_dpar` to each formula in a segment
#'
#' @param components The segment's components from `get_segment_components()`.
#' @noRd
get_predictors_segment = function(components, segment, family, data, par_x, check_rank = TRUE, previous = NULL) {
  checkmate::assert_list(components)
  checkmate::assert_int(segment, lower = 1)
  checkmate::assert_class(family, "mcpfamily")
  checkmate::assert_data_frame(data)
  checkmate::assert_string(par_x)

  #############################
  # DISTRIBUTIONAL PARAMETERS #
  #############################
  dpar_pars = list()
  for (dpar in family$dpar_specs$dpar) {
    spec = get_dpar_spec(family, dpar)
    dpar_form = components[[dpar]]

    # An implicit dpar includes an intercept in segment 1 and is joined in
    # later segments until the user includes another dpar intercept.
    if (is.null(dpar_form) && spec$implicit && segment == 1) {
      dpar_form = stats::as.formula("~1", env = environment(components$mu))
    } else if (!is.null(dpar_form) && segment == 1 && spec$require_initial_predictor && is_empty_initial_predictor(dpar_form)) {
      stop(
        "`", dpar, "(0)` cannot be the initial predictor for family = ", family$family,
        "(). Include an initial predictor, such as `", dpar, "(1)` or `", dpar, "(0 + x)`.",
        call. = FALSE
      )
    }

    if (!is.null(dpar_form)) {
      dpar_pars[[dpar]] = get_shared_predictors(
        data, dpar_form, segment, dpar = dpar, par_x, NULL, check_rank,
        design_id = paste("population", dpar, segment, sep = ":"), previous = previous
      )
    }
  }

  #########
  # AR/MA #
  #########
  arma_pars = list()
  for (component in intersect(c("ar", "ma"), names(components))) {
    arma = components[[component]]
    if (length(get_group_terms(arma$form)) > 0)
      stop(
        "Group-level effects inside ", component,
        "() are not currently supported. Found one in segment ", segment, "."
      )

    # Expand one formula into a separate regression parameter for each lag.
    # Order 0 turns off all lags, coded as lag 1 with formula ~ 0.
    arma_pars[[component]] = lapply(seq_len(max(1, arma$order)), function(order) {
      get_shared_predictors(
        data, arma$form, segment, component, par_x, order, check_rank,
        design_id = paste("population", component, if (arma$order == 0) 0 else order, segment, sep = ":"),
        previous = previous
      )
    }) %>%
      dplyr::bind_rows()
  }

  ##########
  # RETURN #
  ##########
  dplyr::bind_rows(
    dplyr::bind_rows(dpar_pars),
    dplyr::bind_rows(arma_pars)
  )
}


#' Apply `get_predictors_segment` to all segments of a model.
#'
#' @noRd
get_predictor_tables = function(model, data, family, par_x, check_rank = TRUE) {
  # Parse the components of each segment once (mu, other dpars, ar, ma)
  components = lapply(model, function(form) get_segment_components(get_rhs(form), family))
  n_segments = length(components)
  assert_offsets(components, family)

  # Resolve each component against earlier occurrences, including chains
  parsed_predictors = NULL
  for (segment in seq_len(n_segments)) {
    parsed_predictors = dplyr::bind_rows(parsed_predictors,
      get_predictors_segment(components[[segment]], segment, family, data, par_x,
        check_rank, previous = parsed_predictors))
  }
  parsed_predictors = dplyr::arrange(parsed_predictors, .data$dpar, .data$segment)

  # An implicit dpar without a formula anywhere is the parameter itself, e.g., sigma_1
  # rather than the link-scale intercept sigma_Intercept_1 of a modeled sigma().
  unmodeled = parsed_predictors$dpar %in% setdiff(family$dpar_specs$dpar, c("mu", unlist(lapply(components, names))))
  parsed_predictors$name[unmodeled] = paste0(parsed_predictors$dpar[unmodeled], "_1")
  predictors = parsed_predictors %>%
    dplyr::filter(.data$par_type != "offset") %>%
    dplyr::mutate(matrix_col = dplyr::row_number())

  assert_unique_predictor_names(definition_rows(predictors))

  # AR and MA transform the same response history, so one threshold applies to all segments.
  thresholds = unique(unlist(lapply(components, function(x) stats::na.omit(c(x$ar$threshold, x$ma$threshold)))))
  if (length(thresholds) > 1)
    stop("Only one `threshold` value can be specified across all ar() and ma() terms.")
  threshold = if (length(thresholds) == 1) thresholds else 0.1  # Default
  predictors$threshold = ifelse(predictors$dpar %in% c("ar", "ma"), threshold, NA_real_)

  # How long each population term applies. x-terms and intercepts join later segments until the
  # next intercept for the same dpar and AR/MA lag, which starts a disjoined segment. AR/MA terms
  # also end in a segment that does not include their lag. Other terms apply in their own segment only.
  population_key = paste(predictors$dpar, predictors$order)
  is_intercept = predictors$par_type == "Intercept"
  included_arma = get_included_arma(components)
  arma_lags = dplyr::distinct(predictors[predictors$dpar %in% c("ar", "ma"), c("dpar", "order")])
  arma_resets = lapply(seq_len(nrow(arma_lags)), function(i) {
    included = included_arma$segment[included_arma$dpar == arma_lags$dpar[i] & included_arma$order >= arma_lags$order[i]]
    tibble::tibble(key = paste(arma_lags$dpar[i], arma_lags$order[i]), segment = setdiff(seq_len(n_segments), included))
  })
  population_resets = dplyr::bind_rows(
    tibble::tibble(key = population_key[is_intercept], segment = predictors$segment[is_intercept]),
    arma_resets
  )
  predictors$next_segment = get_next_segment(
    population_key, predictors$segment, joins = is_intercept | predictors$x_degree > 0,
    population_resets, n_segments
  )

  # Parse group blocks segment by segment so `same()` can find its source
  definitions = NULL
  for (segment in seq_len(n_segments)) {
    definitions = dplyr::bind_rows(definitions,
      get_predictor_group_definitions_segment(components[[segment]], segment, family,
        data, par_x, check_rank, previous = definitions))
  }

  group_predictors = definitions
  if (nrow(definitions) > 0) {
    # Rows without coefficients are `(0 | group)` turn-offs, which only mark resets below.
    # Link each coefficient to its population counterpart, if there is one.
    group_predictors = definitions %>%
      dplyr::filter(.data$has_coefficients) %>%
      dplyr::mutate(
        population_name = dplyr::if_else(
          .data$population_name %in% predictors$name,
          .data$population_name,
          NA_character_
        ),
        matrix_col = nrow(predictors) + dplyr::row_number()
      ) %>%
      dplyr::select(
        "population_name", "name", "group_col", "segment", "dpar",
        "sd_name", "par_type", "matrix_name", "order",
        "x_degree", "design_id", "design_col", "matrix_col", "matrix_data",
        "design_spec", "definition_segment"
      )

    # How long each group term applies. Group x-terms join later segments until a group intercept
    # or `(0 | group)` for the same dpar and grouping column, or a population intercept for that dpar.
    # Group intercepts and other terms apply in their own segment only.
    group_key = function(x) paste(x$dpar, x$group_col)
    is_reset = !definitions$has_coefficients | definitions$par_type %in% "Intercept"
    disjoined = dplyr::distinct(definitions, .data$dpar, .data$group_col) %>%
      dplyr::inner_join(dplyr::distinct(predictors[is_intercept, c("dpar", "segment")]), by = "dpar", relationship = "many-to-many")
    group_resets = dplyr::bind_rows(
      tibble::tibble(key = group_key(definitions)[is_reset], segment = definitions$segment[is_reset]),
      tibble::tibble(key = group_key(disjoined), segment = disjoined$segment)
    )
    group_predictors$next_segment = get_next_segment(
      group_key(group_predictors), group_predictors$segment,
      joins = group_predictors$x_degree > 0, group_resets, n_segments
    )
  }

  # Store each fitted component specification once, not on every coefficient.
  design_specs = collect_design_specs(parsed_predictors, definitions)
  predictors = predictors %>%
    dplyr::select(-"design_spec")
  if ("design_spec" %in% names(group_predictors))
    group_predictors = dplyr::select(group_predictors, -"design_spec")
  if ("name" %notin% names(group_predictors))
    group_predictors = tibble::tibble(
      population_name = character(), name = character(), group_col = character(), segment = integer(),
      dpar = character(), sd_name = character(), par_type = character(), matrix_name = character(),
      order = integer(), x_degree = integer(), design_id = character(), design_col = integer(),
      matrix_col = integer(), matrix_data = list(), definition_segment = integer(), next_segment = integer()
    )

  list(
    components = components,
    predictors = predictors,
    group_predictors = group_predictors,
    design_specs = design_specs
  )
}
