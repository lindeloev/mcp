# RowSums of element-wise products of matrices. Like an inner product, just vectorized for many y.
# This ensures R <--> JAGS code compatibility.
inprod = function(x, y) {
  checkmate::assert_matrix(x)
  checkmate::assert_matrix(y)
  rowSums(x * y)
}


# Converts zero-length vectors to NULL. Returns x otherwise
logical0_to_null = function(x) {
  if (length(x) > 0) x
}


# if((a %in% b) == FALSE) --> if(a %notin% b)
`%notin%` = Negate(`%in%`)


# Whether user-supplied JAGS code may have diverged from mcp's R-side model.
has_custom_jags_code = function(fit) {
  isTRUE(fit$.internal$custom_jags_code)
}


# Warn before evaluating an R-side model that was generated from formulas rather
# than the custom JAGS code used for fitting.
warn_custom_jags_code = function(fit) {
  if (has_custom_jags_code(fit)) {
    warning(
      "Custom `jags_code` was supplied. mcp's R-side simulations and predictions may not match the fitted JAGS model.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


# List of categorical column names and their unique levels
get_categorical_levels = function(df) {
  checkmate::assert_data_frame(df)
  categorical_cols = colnames(df)[vapply(df, function(x) is.factor(x) || is.logical(x) || is.character(x), logical(1))]
  lapply(df[, categorical_cols, drop = FALSE], unique)
}


# Collapse vector into a string joined by " and "
and_collapse = function(x) {
  checkmate::assert_character(x)
  paste0(x, collapse = " and ")
}


# Merge partial diagnostics with the defaults.
resolve_diagnostics = function(diagnostics = list()) {
  defaults = list(
    rhat = 1.01,
    ess_bulk = 400,
    ess_tail = 400,
    ar = 0.10,
    ma = 0.10
  )
  if (identical(diagnostics, FALSE))
    return(lapply(defaults, function(x) NULL))
  if (!is.list(diagnostics))
    stop("`diagnostics` must be a named list or FALSE.")
  if (length(diagnostics) > 0 && (is.null(names(diagnostics)) || any(names(diagnostics) == "")))
    stop("Every element of `diagnostics` must be named.")

  unknown = setdiff(names(diagnostics), names(defaults))
  if (length(unknown) > 0)
    stop("Unknown diagnostic(s): ", and_collapse(unknown), ".")
  defaults[names(diagnostics)] = diagnostics
  defaults
}


# Check MCMC diagnostic criteria for a summary data.frame.
# - results: A data.frame with a `parameter` column and optionally `rhat`, `ess_bulk`, `ess_tail`.
# - diagnostics: Resolved or partial diagnostics list.
# - fixed_pars: Parameter names for constant params (e.g. "x_2" for x_2 = 5)
# - check_rhat: Logical. Whether Rhat should be evaluated (FALSE if nchain < 2).
check_convergence = function(results, diagnostics = list(), fixed_pars = character(), check_rhat = TRUE) {
  diagnostics = resolve_diagnostics(diagnostics)
  if (nrow(results) == 0)
    return(list(bad = logical(0), bad_names = character(0), n_bad = 0L, thresholds = ""))

  base_name = sub("\\[.*\\]$", "", results$parameter)
  is_fixed = results$parameter %in% fixed_pars | base_name %in% fixed_pars

  check_rhat = check_rhat && !is.null(diagnostics$rhat) && "rhat" %in% names(results)
  check_bulk = !is.null(diagnostics$ess_bulk) && "ess_bulk" %in% names(results)
  check_tail = !is.null(diagnostics$ess_tail) && "ess_tail" %in% names(results)

  bad = rep(FALSE, nrow(results))
  if (check_rhat)
    bad = bad | is.na(results$rhat) | results$rhat > diagnostics$rhat
  if (check_bulk)
    bad = bad | is.na(results$ess_bulk) | results$ess_bulk < diagnostics$ess_bulk
  if (check_tail)
    bad = bad | is.na(results$ess_tail) | results$ess_tail < diagnostics$ess_tail

  bad = bad & !is_fixed

  thresholds = c(
    if (check_rhat) paste0("rhat > ", diagnostics$rhat),
    if (check_bulk) paste0("ess_bulk < ", diagnostics$ess_bulk),
    if (check_tail) paste0("ess_tail < ", diagnostics$ess_tail)
  )

  list(
    bad = bad,
    bad_names = results$parameter[bad],
    n_bad = sum(bad),
    thresholds = paste(thresholds, collapse = " or ")
  )
}


# Warn about poorly mixed posterior chains.
# Thresholds (rhat > 1.01, bulk/tail ESS < 400) follow Vehtari et al. (2021).
# - mcmc_post: An mcmc.list of posterior draws.
# - diagnostics: A resolved diagnostics configuration.
# - fixed_pars: Parameter names intentionally fixed to constants.
warn_nonconvergence = function(mcmc_post, diagnostics = list(), fixed_pars = character()) {
  diagnostics = resolve_diagnostics(diagnostics)
  if (is.null(unlist(diagnostics[c("rhat", "ess_bulk", "ess_tail")])))
    return(invisible(NULL))

  check_rhat = !is.null(diagnostics$rhat) && coda::nchain(mcmc_post) >= 2

  results = posterior::summarise_draws(
    posterior::as_draws_df(mcmc_post),
    rhat = posterior::rhat,
    ess_bulk = function(x) suppressWarnings(posterior::ess_bulk(x)),
    ess_tail = function(x) suppressWarnings(posterior::ess_tail(x))
  ) %>%
    dplyr::rename(parameter = "variable")

  conv = check_convergence(results, diagnostics, fixed_pars = fixed_pars, check_rhat = check_rhat)
  if (conv$n_bad == 0)
    return(invisible(NULL))

  warning(
    "Some parameters may not have converged well:\n",
    "  * ", conv$thresholds, ": ",
    and_collapse(conv$bad_names), "\n",
    "Inspect `summary(fit)` and `plot_pars(fit)`, and consider increasing ",
    "`iter`/`warmup`, set inits manually, or simplifying the model.",
    call. = FALSE
  )
  invisible(NULL)
}


# Get column names of predictor designs ordered by design-matrix column
# - predictors: Output of get_predictors()
# - group_effects: Output of get_group_effects()
get_predictor_design_names = function(predictors, group_effects = NULL) {
  # Extract group-level predictor columns
  group_predictors = if (is.null(group_effects) || "matrix_col" %notin% names(group_effects)) {
    tibble::tibble(matrix_col = integer(), name = character(), segment = integer())
  } else {
    group_effects[group_effects$part == "predictor", c("matrix_col", "name", "segment"), drop = FALSE]
  }

  # Combine population and group-level predictors
  design = dplyr::bind_rows(
    tibble::tibble(
      matrix_col = predictors$matrix_col,
      name = predictors$name,
      segment = predictors$segment
    ),
    group_predictors
  )

  # Disambiguate duplicate names across occurrences and ensure uniqueness
  design = design[order(design$matrix_col), , drop = FALSE]
  duplicate = duplicated(design$name) | duplicated(design$name, fromLast = TRUE)
  design$name[duplicate] = paste0(design$name[duplicate], "_occurrence_", design$segment[duplicate])
  make.unique(design$name)
}


#' Create a model matrix from the predictor metadata tables
#'
#' Combines population and group-level predictor design columns.
#'
#' @param predictors The output of `get_predictors()`.
#' @param group_effects The output of `get_group_effects()`.
#' @return A matrix with one column for each population or group-level
#'   predictor coefficient.
#' @noRd
get_predictor_matrix = function(predictors, group_effects = NULL) {
  checkmate::assert_data_frame(predictors)
  group_predictors = if (is.null(group_effects) || "matrix_col" %notin% names(group_effects)) {
    tibble::tibble(
      matrix_col = integer(), name = character(), matrix_data = list()
    )
  } else {
    group_effects[group_effects$part == "predictor", , drop = FALSE]
  }
  # Each coefficient contributes one column; both sources share matrix indices.
  column_order = order(c(predictors$matrix_col, group_predictors$matrix_col))
  matrix_data = c(predictors$matrix_data, group_predictors$matrix_data)[column_order]
  column_names = get_predictor_design_names(predictors, group_effects)

  predictor_matrix = do.call(cbind, matrix_data)
  colnames(predictor_matrix) = column_names
  predictor_matrix
}


#' Evaluated-draw identity helpers
#'
#' Within one call to `pp_eval()`, `.row` identifies a row in the evaluated data.
#' Together, `.draw` and `.row` must identify exactly one evaluated value.
#' Plotting groups, predictor values, and group-level columns are metadata for
#' that row and never determine which draws may be summarised together.
#'
#' @param draws Evaluated draws in tidy format.
#' @param type Name of the evaluated-value column, matching `pp_eval(type)`.
#' @return `draws`, invisibly.
#' @noRd
validate_eval_draws = function(draws, type) {
  checkmate::assert_data_frame(draws)
  checkmate::assert_string(type)
  assert_data_cols(draws, c(".draw", ".row", type))

  if (anyNA(draws$.draw) || anyNA(draws$.row))
    stop_github("Evaluated draws contain missing `.draw` or `.row` keys.")

  draw_ids = unique(draws$.draw)
  data_rows = unique(draws$.row)
  draw_index = match(draws$.draw, draw_ids)
  row_index = match(draws$.row, data_rows)
  keys = draw_index + length(draw_ids) * (row_index - 1L)
  if (anyDuplicated(keys))
    stop_github("Evaluated draws must contain one `", type, "` value per `.draw` and `.row`.")

  if (nrow(draws) != length(draw_ids) * length(data_rows))
    stop_github("Every `.row` must contain the same complete set of posterior draws.")

  invisible(draws)
}


#' Convert evaluated draws from tidy to matrix format
#'
#' Converts the output of `pp_eval(fit, draws_format = "tidy")` to an `N_draws`
#' by `N_evaluation_rows` matrix. Explicit keys prevent duplicate or incomplete
#' draw/row combinations from being silently reshaped.
#'
#' @inheritParams validate_eval_draws
#' @param data_rows Optional `.row` values to use as matrix columns.
#' @return A numeric matrix.
#' @noRd
tidy_to_matrix = function(draws, type, data_rows = NULL) {
  checkmate::assert_string(type)
  assert_data_cols(draws, c(".draw", ".row", type))

  if (is.null(data_rows))
    data_rows = sort(unique(draws$.row))
  if (anyDuplicated(data_rows))
    stop_github("Requested `.row` values must be unique.")

  missing_rows = setdiff(data_rows, unique(draws$.row))
  if (length(missing_rows) > 0)
    stop_github("Requested evaluation rows are absent: ", paste(missing_rows, collapse = ", "), ".")
  draws = dplyr::filter(draws, .data$.row %in% data_rows)

  draw_ids = sort(unique(draws$.draw))
  result = matrix(NA_real_, nrow = length(draw_ids), ncol = length(data_rows), dimnames = list(NULL, as.character(data_rows)))
  matrix_rows = match(draws$.draw, draw_ids)
  matrix_cols = match(draws$.row, data_rows)
  result[cbind(matrix_rows, matrix_cols)] = draws[[type]]
  result
}


#' Compute evaluated-draw quantiles without pooling evaluation rows
#'
#' @inheritParams validate_eval_draws
#' @param quantiles Vector of quantiles between zero and one.
#' @param na.rm Whether to remove missing values before computing quantiles.
#' @param keep Evaluation-row metadata to rejoin after summarising.
#' @return A tibble with `.row`, `quantile`, the `type` column, and requested metadata.
#' @noRd
get_quantiles = function(draws, quantiles, type, keep = NULL, na.rm = FALSE) {
  keep = unique(keep)
  assert_data_cols(draws, c(".row", type, keep))
  grid = draws %>% dplyr::select(dplyr::all_of(c(".row", keep))) %>% dplyr::distinct()
  if (anyDuplicated(grid$.row))
    stop_github("Evaluation-row metadata differs across draws for the same `.row`.")

  result = draws %>%
    dplyr::group_by(.data$.row) %>%
    dplyr::reframe(quantile = quantiles,
                   !!type := stats::quantile(.data[[type]], probs = quantiles, names = FALSE, na.rm = na.rm))

  dplyr::left_join(result, grid, by = ".row", relationship = "many-to-one")
}


# Invert the posterior predictive mixture CDF for a single quantile probability `p`.
# The mixture quantile lies between the smallest and largest per-draw quantiles.
find_mixture_quantile = function(family, dpars, data, p) {
  mix_cdf = function(y) mean(family$r$cdf(y, dpars, data))
  bracket = range(family$r$quantile(p, dpars, data))
  if (bracket[1] == bracket[2])
    return(bracket[1])

  # Continuous: root-finding. extendInt guards against rounding at the bracket ends.
  if (!isTRUE(family$response$is_discrete))
    return(stats::uniroot(function(y) mix_cdf(y) - p, bracket, extendInt = "upX", tol = .Machine$double.eps)$root)

  # Discrete: binary search for the smallest integer with mix_cdf(y) >= p.
  # Above 2^53, doubles cannot represent all integers, so stop at a relative tolerance.
  lo = bracket[1]
  hi = bracket[2]
  while (hi - lo > max(1, hi * 1e-12)) {
    mid = floor((lo + hi) / 2)
    if (mix_cdf(mid) >= p) {
      hi = mid
    } else {
      lo = mid + 1
    }
  }
  if (lo < hi && mix_cdf(lo) < p) hi else lo
}


# Compute exact predictive quantiles via mixture CDF inversion
#
# Evaluates the exact posterior predictive mixture CDF at each data row
# and finds the target quantiles via root-finding (uniroot for continuous,
# integer search for discrete count families).
get_mixture_quantiles = function(draws, quantiles, family, keep = NULL, rate = FALSE, dpars = attr(draws, "dpars"), response_data = attr(draws, "response_data")) {
  keep = unique(keep)
  grid = draws %>% dplyr::select(dplyr::all_of(c(".row", keep))) %>% dplyr::distinct()

  # Fallback to empirical quantiles if no cdf method or dpars available
  if (is.null(family$r$cdf) || is.null(dpars))
    return(get_quantiles(draws, quantiles, ".prediction", keep = keep))

  # Split draws by evaluation row to process one data point across all posterior draws
  row_indices = split(seq_len(nrow(draws)), draws$.row)
  unique_rows = as.integer(names(row_indices))

  res_list = lapply(seq_along(row_indices), function(i) {
    idx = row_indices[[i]]
    data_row_id = unique_rows[i]

    # Slice distributional parameters for this data row across all MCMC draws
    dp_slice = lapply(dpars, function(v) v[idx])
    rd_slice = if (!is.null(response_data)) lapply(response_data, function(v) v[idx]) else NULL

    # Invert the mixture CDF for each requested quantile probability and convert to the rate scale
    q_vals = vapply(quantiles, function(q) find_mixture_quantile(family, dp_slice, rd_slice, q), numeric(1))
    q_vals = family$response$observed(q_vals, lapply(rd_slice, `[`, 1), rate)

    df_res = data.frame(
      quantile = quantiles,
      .prediction = q_vals,
      stringsAsFactors = FALSE
    )
    df_res$.row = data_row_id
    df_res
  })

  # Combine and re-attach predictor / grouping metadata
  result = dplyr::bind_rows(res_list)
  dplyr::left_join(result, grid, by = ".row", relationship = "many-to-one")
}


#' Print mcplist
#'
#' Shows a list in a more condensed format using `str(list)`.
#' @aliases print.mcplist
#' @inheritParams print.mcpfit
#' @return `x`, invisibly.
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples 
#' print(demo_fit$model)
#' print(demo_fit$prior)
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
  invisible(x)
}


#' Nice Printing of Multiline Texts
#'
#' Useful for `print(fit$jags_code)`, `print(mcp_demo$example_code)`, etc.
#'
#' @aliases print.mcptext
#' @param x Character, often with newlines.
#' @param ... Must be empty. Reserved for future use.
#' @return `x`, invisibly.
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
  cat(x, "\n", sep = "")
  invisible(x)
}
