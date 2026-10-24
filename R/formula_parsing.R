# ABOUT: Formula-structure helpers. These parse mcp's segment-formula syntax
# (tildes, response, change point, and right-hand-side sections) without
# building any design matrices.
# ------------------------------------------------------------

# Takes any formula-like input (formula or string) and returns a formula
to_formula = function(form) {
  checkmate::assert(
    checkmate::check_character(form, min.len = 1, max.len = 3),
    checkmate::check_formula(form),
    .var.name = "form"
  )
  if (is.character(form)) {
    # Add tilde
    if (!stringr::str_detect(form, "^(\\s|)~")) {
      form = paste0("~", form)
    }
    form = stats::as.formula(form)
  }

  form
}


# Convert formula to string using deparse1() to guarantee a single string.
# One-sided formulas get a space after the tilde ("~ 0 + x"), like two-sided ones.
formula_to_char = function(form) {
  checkmate::assert_formula(form)
  sub("^~", "~ ", deparse1(form))
}


# Returns the right-hand-side of a formula, excluding the change point in cp()
get_rhs = function(form) {
  checkmate::assert_formula(form)
  rhs = split_cp(form[[length(form)]])$rhs
  stats::as.formula(call("~", rhs), env = environment(form))
}


# Check whether an expression is a formula call, e.g., `y ~ x`
is_tilde_call = function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("~"))
}


# Check whether an expression is a call to cp()
is_cp_call = function(expr) {
  is.call(expr) && identical(deparse1(expr[[1]]), "cp")
}


# Names of functions called in an expression, except namespaced calls like splines::ns(x)
called_functions = function(expr) {
  if (!is.call(expr))
    return(character())
  head = expr[[1]]
  is_namespaced = is.call(head) && deparse1(head[[1]]) %in% c("::", ":::")
  c(if (is.symbol(head) && !is_namespaced) as.character(head), unlist(lapply(as.list(expr)[-1], called_functions)))
}


# Recursively check whether an expression contains any call to cp()
contains_cp_call = function(expr) {
  if (is_cp_call(expr))
    return(TRUE)
  if (!is.call(expr))
    return(FALSE)

  any(vapply(as.list(expr)[-1], contains_cp_call, logical(1)))
}


# Split a segment RHS expression into the cp() content (NULL if absent) and the remaining terms
# E.g., `cp(1 + (1 | id)) + 0 + x` --> list(cp = `1 + (1 | id)`, rhs = `0 + x`).
split_cp = function(rhs) {
  leaves = unpack_additive(rhs)
  is_cp = vapply(leaves, is_cp_call, logical(1))
  if (!any(is_cp))
    return(list(cp = NULL, rhs = rhs))
  if (sum(is_cp) > 1)
    stop("Only one `cp()` allowed in each segment formula. Got `", deparse1(rhs), "`.", call. = FALSE)

  cp = leaves[is_cp][[1]]
  if (length(cp) != 2)
    stop("`cp()` takes exactly one formula, e.g., `cp(1 + (1 | id))`. Got `", deparse1(cp), "`.", call. = FALSE)

  rest = leaves[!is_cp]
  list(
    cp = cp[[2]],
    rhs = if (length(rest) == 0) 1 else Reduce(function(a, b) call("+", a, b), rest)
  )
}




# Classify the variables in segment formulas by their role in the model. Returns a list with
# `response` and its `auxiliary` columns (e.g., `trials(N)`), and the right-hand-side
# `predictor`, `group` (grouping factors, including in `cp()`), `offset`, and `series` variables.
# E.g., `y ~ 1 + x + offset(log(t)) + (1 | id)` gives response = "y", predictor = "x",
# group = "id", and offset = "t".
get_formula_vars = function(model) {
  checkmate::assert_true(is_mcpmodel(model), .var.name = "model")
  vars = list(response = character(), auxiliary = character(), predictor = character(),
              group = character(), offset = character(), series = character())
  add = function(role, x) vars[[role]] <<- union(vars[[role]], x)

  # Add variables to `role`, except in grouping, offset, and AR/MA-argument positions.
  # `role = NULL` collects only those special positions.
  walk = function(expr, role) {
    if (!is.call(expr)) {
      if (!is.null(role))
        add(role, all.vars(expr))
      return(invisible(NULL))
    }
    fun = deparse1(expr[[1]])
    if (fun %in% c("offset", "stats::offset")) {
      add("offset", all.vars(expr))
    } else if (fun %in% c("|", "||")) {
      walk(expr[[2]], role)
      add("group", all.vars(expr[[3]]))
    } else if (fun %in% c("ar", "ma")) {
      arma = unpack_arma(expr)
      if (!is.na(arma$series))
        add("series", arma$series)
      walk(arma$formula, role)
    } else {
      for (arg in as.list(expr)[-1])
        walk(arg, role)
    }
  }

  for (i in seq_along(model)) {
    form = model[[i]]
    lhs = if (length(form) == 3) form[[2]] else NULL
    if (i == 1 && is.call(lhs) && identical(lhs[[1]], as.name("|"))) {
      add("response", all.vars(lhs[[2]]))
      add("auxiliary", all.vars(lhs[[3]]))
    } else if (i == 1) {
      add("response", all.vars(lhs))
    } else {
      walk(lhs, role = NULL)  # Deprecated change point on the left-hand side, e.g., `1 + (1 | id) ~ 0 + x`
    }
    walk(form[[length(form)]], role = "predictor")
  }

  if (length(vars$series) > 1)
    stop("Only one `series` column can be specified across all ar() and ma() terms.")

  vars
}

#' Remove group-level or population-level terms from a formula
#'
#' WARNING: removes response side from the formula
#'
#' @param form A formula
#' @param remove Either "group" or "population". These are removed.
#' @return A formula
#' @noRd
remove_terms = function(form, remove) {
  checkmate::assert_formula(form)
  remove = rlang::arg_match0(remove, c("group", "population"))

  # Find terms with "|"
  attrs = attributes(stats::terms(form))
  term.labels = attrs$term.labels
  is_group = stringr::str_detect(term.labels, "\\|")

  # Add parenthesis back to them
  term.labels[is_group] = paste0("(", term.labels[is_group], ")")

  # Remove non-matching types
  if (remove == "group") {
    term.labels = term.labels[!is_group]
    # base::terms() omits offset() from term.labels; re-attach it for population formulas
    if (!is.null(attrs$offset)) {
      offset_terms = vapply(attrs$offset, function(i) deparse1(attrs$variables[[i + 1]]), character(1))
      term.labels = c(term.labels, offset_terms)
    }
    term.labels = c(attrs$intercept, term.labels)  # Add intercept indicator
  } else if (remove == "population") {
    term.labels = term.labels[is_group]
  }

  # Build formula from terms and return
  if (length(term.labels) == 0) {
    return(NULL)
  } else {
    formula_terms = paste0(term.labels, collapse = " + ")
    formula_str = paste0("~", formula_terms)
    # Rebuild in the caller's environment so local transformations still resolve.
    return(stats::as.formula(formula_str, env = environment(form)))
  }
}


#' Get formula inside a wrapper
#'
#' @param term E.g., "mu(1 + x)", "sigma(0 + I(x^2))", etc.
#' @param env Environment in which the term was originally defined.
#' @return char formula with the content inside the brackets.
#' @noRd
get_term_content = function(term, env = parent.frame()) {
  if (length(term) > 1)
    stop("Only one ", sub("\\(.*", "()", term[1]), " allowed in each formula. Got ", and_collapse(term), ".")

  # Get formula inside wrapper
  content = sub("^[^(]*\\((.*)\\)$", "\\1", term)
  if (content == "")
    stop("Empty terms not allowed in the formulas. Found '", term, "'.")

  # Wrapper terms are strings, so carry their source formula environment explicitly.
  stats::as.formula(paste0("~", content), env = env)
}


# Unpack an additive expression into a list of leaf terms
unpack_additive = function(expr) {
  if (is.call(expr) && identical(deparse1(expr[[1]]), "+")) {
    c(unpack_additive(expr[[2]]), unpack_additive(expr[[3]]))
  } else {
    list(expr)
  }
}


# Split a segment's right-hand side into its components: a one-sided formula for each included
# distributional parameter, where bare terms belong to mu, and the unpacked ar()/ma() terms
# (see `unpack_arma()`, with the formula as `form`).
# E.g., `~ 1 + x + sigma(1 + x) + ar(2)` --> list(mu = ~1 + x, sigma = ~1 + x, ar = list(order = 2, form = ~1, ...)).
# As in other R formulas, a segment without mu terms includes an intercept: `~ sigma(1)` has mu = ~1.
get_segment_components = function(form_rhs, family) {
  checkmate::assert_formula(form_rhs)
  checkmate::assert_class(family, "mcpfamily")
  as_form = function(expr) stats::as.formula(call("~", expr), env = environment(form_rhs))
  sum_exprs = function(exprs) Reduce(function(a, b) call("+", a, b), exprs)

  # Unpack top-level additive terms and their wrapper names
  leaves = unpack_additive(form_rhs[[2]])
  leaf_names = vapply(leaves, function(x) if (is.call(x)) deparse1(x[[1]]) else "", character(1))
  wrappers = c(family$dpar_specs$dpar, "ar", "ma")

  unsupported = setdiff(intersect(leaf_names, known_dpar_wrappers()), family$dpar_specs$dpar)
  if (length(unsupported) > 0) {
    family_call = paste0(family$family, "()")
    stop(
      paste0("`", unsupported, "()`", collapse = " and "), " is not a distributional parameter for family = ", family_call, ". ",
      "Available parameters: ", and_collapse(paste0("`", family$dpar_specs$dpar, "`")), "."
    )
  }
  repeated = unique(leaf_names[duplicated(leaf_names) & leaf_names %in% wrappers])
  if (length(repeated) > 0)
    stop("Only one ", and_collapse(paste0("`", repeated, "()`")), " allowed in each segment.", call. = FALSE)

  # Bare terms are added to the contents of an explicit mu()
  is_bare = leaf_names %notin% wrappers
  components = list()
  for (i in which(!is_bare)) {
    if (leaf_names[i] %in% c("ar", "ma")) {
      arma = unpack_arma(leaves[[i]])
      arma$form = as_form(arma$formula)
      components[[leaf_names[i]]] = arma[c("order", "form", "threshold", "series")]
    } else if (length(leaves[[i]]) != 2) {
      stop("`", leaf_names[i], "()` takes exactly one formula, e.g., `", leaf_names[i], "(1 + x)`. Got `", deparse1(leaves[[i]]), "`.", call. = FALSE)
    } else {
      components[[leaf_names[i]]] = as_form(leaves[[i]][[2]])
    }
  }
  mu_leaves = c(if (!is.null(components$mu)) unpack_additive(components$mu[[2]]), leaves[is_bare])
  components$mu = as_form(if (length(mu_leaves) == 0) 1 else sum_exprs(mu_leaves))
  components
}


# Canonicalize a segment formula into a uniform AST structure
# Returns a list: response (expr), cp (expr or NULL for segment 1), rhs (expr), and form (canonical formula)
# E.g., `~ cp((1 | id)) + 0 + x` gives cp = `(1 | id)` and rhs = `0 + x`.
canonicalize_segment = function(form, i, default_response = NULL) {
  form = to_formula(form)
  if (i > 1 && length(form) == 3) {
    if (identical(form[[2]], default_response))
      stop("Error in segment ", i, ": segments 2+ use the response from segment 1. Remove `", deparse1(form[[2]]), " ~`, e.g., `~ 0 + x`.", call. = FALSE)
    form = legacy_cp_lhs(form, i, default_response)
  }

  split = split_cp(form[[length(form)]])
  if (contains_cp_call(split$rhs))
    stop("Error in segment ", i, ": `cp()` must be a separate term in the formula, e.g., `~ cp((1 | id)) + 0 + x`. Got `", deparse1(form), "`.", call. = FALSE)

  # Segment 1 must define the response and cannot contain a change point
  if (i == 1) {
    if (length(form) == 2)
      stop("No response variable in segment 1.", call. = FALSE)
    if (is_tilde_call(form[[2]]))
      stop("The first segment must have exactly one tilde. Got two.", call. = FALSE)
    if (!is.null(split$cp))
      stop("Segment 1 has no change point, so `cp()` is only allowed in segments 2+.", call. = FALSE)

    return(list(response = form[[2]], cp = NULL, rhs = split$rhs, form = form))
  }

  # Segment > 1: the response is inherited from segment 1 and the change point is in cp(), defaulting to cp(1)
  cp = if (is.null(split$cp)) 1 else split$cp
  list(response = default_response, cp = cp, rhs = split$rhs, form = form)
}


#' Unpacks y variable name
#'
#' @inheritParams mcp
#' @param form_y Character representation of formula
#' @param i Segment number
#' @return A one-row tibble with the response and auxiliary-data columns.
#' @noRd
unpack_y = function(form_y, i, family, env = parent.frame()) {
  if (is.language(form_y))
    form_y = deparse1(form_y)

  family_aux = names(family$response$auxiliary)
  aux_names = unique(c("trials", "weights", family_aux))
  response = stats::setNames(as.list(rep(NA_character_, length(aux_names) + 1)), c("y", aux_names))

  # If NA and not segment 1, just return empty
  if (is.na(form_y)) {
    if (i == 1)
      stop("A response must be defined in segment 1, e.g., 'y ~ 1'")

    return(tibble::as_tibble(response))
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
  if (!is.symbol(str2lang(y_col)))
    stop("The response must be a data column name. Got '", y_col, "' in segment ", i, ". Transform it in the data before fitting, e.g., `data$y_log = log(data$y)`.", call. = FALSE)
  response$y = y_col

  term_labels = character()
  if (length(y_split) == 2) {
    rhs = y_split[2]
    term_labels = attr(stats::terms(to_formula(rhs)), "term.labels")
    ok_terms = if (length(family_aux) == 0) rep(FALSE, length(term_labels)) else
      vapply(term_labels, function(term) any(stringr::str_detect(term, paste0("^", family_aux, "\\("))), logical(1))
    if (!all(ok_terms))
      stop(
        "Only ", if (length(family_aux) == 0) "no terms are" else and_collapse(paste0("`", family_aux, "()`")),
        " allowed after the pipe for family = ", family$family, "(). Got '", rhs, "'."
      )
  }

  for (name in family_aux) {
    term_index = stringr::str_detect(term_labels, paste0("^", name, "\\("))
    got_term = any(term_index)
    if (family$response$auxiliary[[name]]$required && !got_term)
      stop("Error in response of segment ", i, ": need a valid ", name, "() specification.")
    if (!got_term)
      next
    if (sum(term_index) > 1)
      stop("Only one ", name, "() term is allowed in segment ", i, ".")

    term = term_labels[term_index]
    content = get_term_content(term, env)
    column = attr(stats::terms(content), "term.labels")
    if (length(column) != 1)
      stop("There must be exactly one term inside ", name, "(). Got ", term, " in segment ", i)
    response[[name]] = column
  }

  tibble::as_tibble(response)
}


#' Takes the content of cp() and returns its properties
#'
#' @param form_cp The expression inside `cp()`, e.g., `1 + (1 | id)`, or `NULL` for segment 1.
#' @param i segment number
#' @return A one-row tibble with columns:
#'   * `cp_group_col`: char or NA. Which data column defines the grouping factor?
#' @noRd
unpack_cp = function(form_cp, i, env = parent.frame()) {
  # Segment 1 has no change point
  if (is.null(form_cp)) {
    return(tibble::tibble(
      cp_group_col = NA
    ))
  }

  form_cp = if (rlang::is_formula(form_cp)) form_cp else stats::as.formula(call("~", form_cp), env = env)

  # Group-level effects
  form_group = remove_terms(form_cp, "population")

  if (!is.null(form_group)) {
    group_terms = attr(stats::terms(form_group), "term.labels")
    if (length(group_terms) > 1)
      stop("Error in segment ", i, " (change point): only one group-level effect is allowed. Found ", deparse1(form_cp))

    group_parts = strsplit(gsub(" ", "", group_terms), "\\|\\|?")[[1]]  # (1 || id) is (1 | id) for intercepts
    if (!group_parts[1] == "1")
      stop("Error in segment ", i, " (change point): Only plain intercepts are allowed in group-level effects, e.g., `cp((1 | id))`.")

    if (!grepl("^[A-Za-z._0-9]+$", group_parts[2]))
      stop("Error in segment ", i, " (change point): invalid grouping-variable format in group-level effect. Got: ", group_parts[2])
  }

  # Population-level effects
  attrs = attributes(stats::terms(remove_terms(form_cp, "group")))
  if (length(attrs$term.labels) > 0)
    stop("Error in segment ", i, " (change point): Only intercepts are allowed in population-level effects, e.g., `cp(1 + (1 | id))`.")

  if (attrs$intercept == 0)
    stop("Error in segment ", i, " (change point): change points always include a population-level intercept. Remove `0 +`, e.g., `cp((1 | id))`.")

  tibble::tibble(
    cp_group_col = if (is.null(form_group)) NA else group_parts[2]
  )
}
