#' RowSums of element-wise products.
#' Like an inner product, just vectorized for many y.
#' This ensures R <--> JAGS code compatibility
#'
#' @aliases inprod
#' @param x A matrix
#' @param y A matrix
#' @return A vector
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
inprod = function(x, y) {
  checkmate::assert_matrix(x)
  checkmate::assert_matrix(y)
  rowSums(x * y)
}


# Converts logical(0) to null. Returns x otherwise
logical0_to_null = function(x) {
  if (length(x) > 0)
    return(x)
  else
    return(NULL)
}


# if((a %in% b) == FALSE) --> if(a %notin% b)
`%notin%` = Negate(`%in%`)


# List of categorical column names and their unique levels
get_categorical_levels = function(df) {
  checkmate::assert_data_frame(df)
  categorical_cols = colnames(df)[sapply(df, class) %in% c("factor", "logical", "character")]
  lapply(df[, categorical_cols, drop = FALSE], unique)
}


# Ask reminder questions for CRAN export
release_questions = function() {
  c(
    #"TEST: Have you run the extensive tests? options(test_mcp_allmodels = TRUE)",
    "TEST: Have you run the test of fits? (uncomment skip() in test-fits-examples.R and helper-fits.R)",
    "TEST: Have you reviewed all example plots? Before running test-fits-examples.R, use `Sys.setenv(MCP_MAKE_PLOT_TEST_FILE = '~/mcp-example-plots.pdf')`.",
    "TEST: Have you run `revdepcheck::revdep_check()`?",


    "DATA: Have you generated a new `demo_fit` if needed?",

    "DOC: Have you built the README plots and checked them? source('vignettes/figures/make_README_plots.R')",
    "DOC: Have you re-built the site using pkgdown::build_site() AFTER deleting caches of articles in 'vignettes/*_cache/'?",
    "DOC: Have you checked all articles and plots after re-building the site?",
    "DOC: Have you run the script to insert the correct logo.png in the HTML meta?"
  )
}


#' Remove varying or population terms from a formula
#'
#' WARNING: removes response side from the formula
#'
#' @aliases remove_terms
#' @keywords internal
#' @noRd
#' @param form A formula
#' @param remove Either "varying" or "population". These are removed.
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
remove_terms = function(form, remove) {
  checkmate::assert_formula(form)
  remove = rlang::arg_match0(remove, c("varying", "population"))

  # Find terms with "|"
  attrs = attributes(stats::terms(form))
  term.labels = attrs$term.labels
  varying_bool = stringr::str_detect(term.labels, "\\|")

  # Add parenthesis back to them
  term.labels[varying_bool] = paste0("(", term.labels[varying_bool], ")")

  # Remove non-matching types
  if (remove == "varying") {
    term.labels = term.labels[!varying_bool]
    term.labels = c(attrs$intercept, term.labels)  # Add intercept indicator
  } else if (remove == "population") {
    term.labels = term.labels[varying_bool]
  }

  # Build formula from terms and return
  if (length(term.labels) == 0) {
    return(NULL)
  } else {
    formula_terms = paste0(term.labels, collapse = " + ")
    formula_str = paste0("~", formula_terms)
    return(stats::as.formula(formula_str, env=globalenv()))
  }
}


#' Takes any formula-like input and returns a formula
#' @aliases to_formula
#' @keywords internal
#' @noRd
#' @param form Formula or character (with or without initial tilde/"~")
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
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


#' Homogonize enumerating strings in mcp
#'
#' Nice for error messages.
#'
#' @aliases collapse_quote
#' @keywords internal
#' @noRd
#' @param x A character vector
#' @return Character
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
and_collapse = function(x) {
  checkmate::assert_character(x)
  paste0(x, collapse = " and ")
}


#' Warn about poorly mixed posterior chains
#'
#' Thresholds (Rhat > 1.01, bulk/tail ESS < 400) follow the recommendations in
#' Vehtari, Gelman, Simpson, Carpenter, & Bürkner (2021). "Rank-normalization,
#' folding, and localization: An improved Rhat for assessing convergence of
#' MCMC". Bayesian Analysis, 16(2), 667-718. \doi{10.1214/20-BA1221}. The same
#' thresholds are used by Stan/`cmdstanr`/`brms`.
#'
#' @aliases warn_nonconvergence
#' @keywords internal
#' @noRd
#' @param mcmc_post An `mcmc.list` of posterior draws.
#' @return `NULL`, invisibly. Called for the warning side-effect.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
warn_nonconvergence = function(mcmc_post) {
  if (coda::nchain(mcmc_post) < 2)
    return(invisible(NULL))  # Rhat/ESS need >= 2 chains

  diagnostics = posterior::summarise_draws(
    posterior::as_draws_df(mcmc_post),
    rhat = posterior::rhat,
    ess_bulk = function(x) suppressWarnings(posterior::ess_bulk(x)),
    ess_tail = function(x) suppressWarnings(posterior::ess_tail(x))
  )

  bad_rhat = diagnostics$variable[!is.na(diagnostics$rhat) & diagnostics$rhat > 1.01]
  bad_ess = diagnostics$variable[
    (!is.na(diagnostics$ess_bulk) & diagnostics$ess_bulk < 400) |
    (!is.na(diagnostics$ess_tail) & diagnostics$ess_tail < 400)
  ]

  if (length(bad_rhat) == 0 && length(bad_ess) == 0)
    return(invisible(NULL))

  warning(
    "Some parameters may not have converged well:\n",
    if (length(bad_rhat) > 0) paste0("  * Rhat > 1.01: ", and_collapse(bad_rhat), "\n"),
    if (length(bad_ess) > 0) paste0("  * ess_bulk or ess_tail < 400: ", and_collapse(bad_ess), "\n"),
    "Inspect `summary(fit)` and `plot_pars(fit)`, and consider increasing ",
    "`iter`/`adapt` or simplifying the model before trusting these results.",
    call. = FALSE
  )
  invisible(NULL)
}


#' Converts formula to string
#'
#' Note: this uses base R and circumvents the length-limitation of `deparse()`
#' and `format()`.
#'
#' @aliases formula_to_char
#' @keywords internal
#' @noRd
#' @param form Any valid formula with any number of tildes.
#' @return A character.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
formula_to_char = function(form) {
  checkmate::assert_formula(form)
  form_char = as.character(form)
  if (length(form_char) == 2 & form_char[1] == "~") {
    return(paste0(form_char, collapse = " "))
  } else if (length(form_char == 3) & form_char[1] == "~") {
    return(paste0(form_char[c(2, 3)], collapse = " ~ "))
  } else {
    stop_github("Could not decode formula ", deparse(form, width.cutoff = 500))
  }
}


#' Returns the right-hand-side of a formula
#'
#' @aliases get_rhs
#' @keywords internal
#' @noRd
#' @param form Formula, e.g. `~x`, `y ~ x` or `y ~ z ~ x`
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs = function(form) {
  checkmate::assert_formula(form)
  if (length(form) == 2) {
    return(form)
  } else if (length(form) == 3) {
    return(form[-2])
  }
}


#' Returns all variables in the predictor parts of an mcpmodel
#'
#' @aliases get_rhs_vars
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return Character vector with unique term names
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs_vars = function(model) {
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")

  model %>%
    lapply(get_rhs) %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()
}


#' Returns grouping-factor variables in predictor group-level terms
#'
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return Character vector of grouping-factor column names.
get_rhs_group_vars = function(model) {
  find_groups = function(expr) {
    if (!is.call(expr))
      return(character())
    if (as.character(expr[[1]]) %in% c("|", "||"))
      return(all.vars(expr[[3]]))
    nested_groups = lapply(rlang::call_args(expr), find_groups)
    unique(unlist(nested_groups))
  }

  model %>%
    lapply(get_rhs) %>%
    lapply(function(form) find_groups(form[[2]])) %>%
    unlist() %>%
    unique()
}

#' Returns all variables in the predictor parts of an mcpmodel
#'
#' @aliases get_model_vars
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return Character vector with unique term names
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_model_vars = function(model) {
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")

  model %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()
}


#' Get names of columns in the predictor design matrix
#'
#' @keywords internal
#' @noRd
#' @param predictors The output of `get_predictors()`.
#' @param group_effects The output of `get_group_effects()`.
#' @return Character vector ordered by design-matrix column.
get_predictor_design_names = function(predictors, group_effects = NULL) {
  group_predictors = if (is.null(group_effects) || "matrix_col" %notin% names(group_effects)) {
    tibble::tibble(matrix_col = integer(), name = character())
  } else {
    group_effects[group_effects$part == "predictor", c("matrix_col", "name"), drop = FALSE]
  }
  design = dplyr::bind_rows(
    tibble::tibble(matrix_col = predictors$matrix_col, name = predictors$code_name),
    group_predictors
  )
  design$name[order(design$matrix_col)]
}


#' Create a model matrix from the predictor metadata tables
#'
#' Combines population and group-level predictor design columns.
#' @aliases get_predictor_matrix
#' @keywords internal
#' @noRd
#' @param predictors The output of `get_predictors()`.
#' @param group_effects The output of `get_group_effects()`.
#' @return A matrix with one column for each population or group-level
#'   predictor coefficient.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_predictor_matrix = function(predictors, group_effects = NULL) {
  checkmate::assert_data_frame(predictors)
  group_predictors = if (is.null(group_effects) || "matrix_col" %notin% names(group_effects)) {
    tibble::tibble(
      matrix_col = integer(), name = character(), matrix_data = list()
    )
  } else {
    group_effects[group_effects$part == "predictor", , drop = FALSE]
  }
  design = dplyr::bind_rows(
    dplyr::transmute(
      predictors,
      matrix_col = .data$matrix_col,
      name = .data$code_name,
      matrix_data = .data$matrix_data
    ),
    dplyr::transmute(
      group_predictors,
      matrix_col = .data$matrix_col,
      name = .data$name,
      matrix_data = .data$matrix_data
    )
  ) %>%
    dplyr::arrange(.data$matrix_col)

  suppressMessages(dplyr::bind_cols(design$matrix_data, .name_repair = "unique")) %>% # Suppress message about lacking column names
    as.matrix() %>%
    magrittr::set_colnames(design$name)
}


#' Evaluated-draw identity helpers
#'
#' Within one call to `pp_eval()`, `data_row` identifies a row in the evaluated data.
#' Together, `.draw` and `data_row` must identify exactly one evaluated value.
#' Plotting groups, predictor values, and varying-effect columns are metadata for
#' that row and never determine which draws may be summarised together.
#'
#' @keywords internal
#' @noRd
#' @param samples Evaluated draws in tidy format.
#' @param type Name of the evaluated-value column, matching `pp_eval(type)`.
#' @return `samples`, invisibly.
validate_eval_draws = function(samples, type) {
  checkmate::assert_data_frame(samples)
  checkmate::assert_string(type)
  assert_data_cols(samples, c(".draw", "data_row", type))

  if (anyNA(samples$.draw) || anyNA(samples$data_row))
    stop_github("Evaluated draws contain missing `.draw` or `data_row` keys.")

  draw_ids = unique(samples$.draw)
  data_rows = unique(samples$data_row)
  draw_index = match(samples$.draw, draw_ids)
  row_index = match(samples$data_row, data_rows)
  keys = draw_index + length(draw_ids) * (row_index - 1L)
  if (anyDuplicated(keys))
    stop_github("Evaluated draws must contain one `", type, "` value per `.draw` and `data_row`.")

  if (nrow(samples) != length(draw_ids) * length(data_rows))
    stop_github("Every `data_row` must contain the same complete set of posterior draws.")

  invisible(samples)
}


#' Convert evaluated draws from tidy to matrix format
#'
#' Converts the output of `pp_eval(fit, samples_format = "tidy")` to an `N_draws`
#' by `N_evaluation_rows` matrix. Explicit keys prevent duplicate or incomplete
#' draw/row combinations from being silently reshaped.
#'
#' @aliases tidy_to_matrix
#' @keywords internal
#' @noRd
#' @inheritParams validate_eval_draws
#' @param data_rows Optional `data_row` values to use as matrix columns.
#' @return A numeric matrix.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
tidy_to_matrix = function(samples, type, data_rows = NULL) {
  checkmate::assert_string(type)
  assert_data_cols(samples, c(".draw", "data_row", type))

  if (is.null(data_rows))
    data_rows = sort(unique(samples$data_row))
  if (anyDuplicated(data_rows))
    stop_github("Requested `data_row` values must be unique.")

  missing_rows = setdiff(data_rows, unique(samples$data_row))
  if (length(missing_rows) > 0)
    stop_github("Requested evaluation rows are absent: ", paste(missing_rows, collapse = ", "), ".")
  samples = dplyr::filter(samples, .data$data_row %in% data_rows)

  draw_ids = sort(unique(samples$.draw))
  result = matrix(NA_real_, nrow = length(draw_ids), ncol = length(data_rows), dimnames = list(NULL, as.character(data_rows)))
  matrix_rows = match(samples$.draw, draw_ids)
  matrix_cols = match(samples$data_row, data_rows)
  result[cbind(matrix_rows, matrix_cols)] = samples[[type]]
  result
}


#' Compute evaluated-draw quantiles without pooling evaluation rows
#'
#' @aliases get_quantiles
#' @keywords internal
#' @noRd
#' @inheritParams validate_eval_draws
#' @param quantiles Vector of quantiles between zero and one.
#' @param keep Evaluation-row metadata to rejoin after summarising.
#' @return A tibble with `data_row`, `quantile`, the `type` column, and requested metadata.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_quantiles = function(samples, quantiles, type, keep = NULL) {
  keep = unique(keep)
  assert_data_cols(samples, c("data_row", type, keep))
  grid = samples %>% dplyr::select("data_row", dplyr::all_of(keep)) %>% dplyr::distinct()
  if (anyDuplicated(grid$data_row))
    stop_github("Evaluation-row metadata differs across draws for the same `data_row`.")

  result = samples %>%
    dplyr::group_by(.data$data_row) %>%
    dplyr::reframe(quantile = quantiles,
                   !!type := stats::quantile(.data[[type]], probs = quantiles, names = FALSE))

  dplyr::left_join(result, grid, by = "data_row", relationship = "many-to-one")
}


#' Print mcplist
#'
#' Shows a list in a more condensed format using `str(list)`.
#' @aliases print.mcplist
#' @inheritParams print.mcpfit
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
print.mcplist = function(x, ...) {
  checkmate::assert_list(x)
  rlang::check_dots_empty()

  # For all-formula list (typically fit$model)
  if (all(sapply(x, rlang::is_formula)) == TRUE) {
    cat(paste0("List of ", length(x), "\n"))
    for (i in x) {
      cat(" $ ")
      print(i)
    }
  } else {
    # Other lists
    utils::str(x, vec.len = Inf, give.head = FALSE, give.attr = FALSE)
  }
}


#' Nice Printing of Multiline Texts
#'
#' Useful for `print(fit$jags_code)`, `print(mcp_demo$call)`, etc.
#'
#' @aliases print.mcptext
#' @param x Character, often with newlines.
#' @param ... Currently ignored.
#' @return NULL
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' mytext = "line1 = 2\n line2 = 'horse'"
#' class(mytext) = "mcptext"
#' print(mytext)
print.mcptext = function(x, ...) {
  checkmate::assert_string(x)
  rlang::check_dots_empty()
  cat(x)
}


# Set model environment to parent.frame() for prettier printing
# and because it was created in a different environment than inteded for use.
fix_model_environment = function(model) {
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")
  for (i in seq_along(model))
    environment(model[[i]]) = globalenv()
  model
}
