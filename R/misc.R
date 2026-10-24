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
# - convergence: A data.frame with a `parameter` column and optionally `rhat`, `ess_bulk`, `ess_tail`.
# - diagnostics: Resolved or partial diagnostics list.
# - constant_pars: Names of parameters with constant priors, e.g., "x_2" for x_2 = 5.
# - check_rhat: Logical. Whether Rhat should be evaluated (FALSE if nchain < 2).
check_convergence = function(convergence, diagnostics = list(), constant_pars = character(), check_rhat = TRUE) {
  diagnostics = resolve_diagnostics(diagnostics)
  if (nrow(convergence) == 0)
    return(list(bad = logical(0), bad_names = character(0), n_bad = 0L, thresholds = ""))

  base_name = sub("\\[.*\\]$", "", convergence$parameter)
  is_constant = convergence$parameter %in% constant_pars | base_name %in% constant_pars

  check_rhat = check_rhat && !is.null(diagnostics$rhat) && "rhat" %in% names(convergence)
  check_bulk = !is.null(diagnostics$ess_bulk) && "ess_bulk" %in% names(convergence)
  check_tail = !is.null(diagnostics$ess_tail) && "ess_tail" %in% names(convergence)

  bad = rep(FALSE, nrow(convergence))
  if (check_rhat)
    bad = bad | is.na(convergence$rhat) | convergence$rhat > diagnostics$rhat
  if (check_bulk)
    bad = bad | is.na(convergence$ess_bulk) | convergence$ess_bulk < diagnostics$ess_bulk
  if (check_tail)
    bad = bad | is.na(convergence$ess_tail) | convergence$ess_tail < diagnostics$ess_tail

  bad = bad & !is_constant

  thresholds = c(
    if (check_rhat) paste0("rhat > ", diagnostics$rhat),
    if (check_bulk) paste0("ess_bulk < ", diagnostics$ess_bulk),
    if (check_tail) paste0("ess_tail < ", diagnostics$ess_tail)
  )

  list(
    bad = bad,
    bad_names = convergence$parameter[bad],
    n_bad = sum(bad),
    thresholds = paste(thresholds, collapse = " or ")
  )
}


# rhat and bulk/tail ESS of each variable in `draws` (an mcmc.list or posterior draws object).
# Unrounded. Shared by the fit-time warning, summary(), fixef(), and ranef().
convergence_table = function(draws) {
  posterior::summarise_draws(
    posterior::as_draws_array(draws),
    rhat = posterior::rhat,
    ess_bulk = function(x) suppressWarnings(posterior::ess_bulk(x)),
    ess_tail = function(x) suppressWarnings(posterior::ess_tail(x))
  ) %>%
    dplyr::rename(parameter = "variable")
}


# The convergence_table() of the posterior draws, cached by mcp(). NULL if the draws changed since then.
cached_convergence = function(fit) {
  convergence = fit$.internal$convergence
  if (!is.null(convergence) && identical(convergence$hash, rlang::hash(mcmclist_draws(fit, error = FALSE))))
    convergence$table
}


# Name the `parameter`s that are `bad`, counting group-level deviations by level instead of listing them,
# e.g., "cp_1 and Intercept_1_id (37 of 300 levels)".
format_bad_pars = function(parameter, bad) {
  base_name = sub("\\[.*\\]$", "", parameter)
  is_group = base_name != parameter
  labels = vapply(unique(base_name[bad]), function(name) {
    is_level = is_group & base_name == name
    if (!any(is_level))
      return(name)
    paste0(name, " (", sum(bad & is_level), " of ", sum(is_level), " levels)")
  }, character(1))
  and_collapse(labels)
}


# Warn about poorly mixed posterior chains.
# Thresholds (rhat > 1.01, bulk/tail ESS < 400) follow Vehtari et al. (2021).
# - convergence: A convergence_table() of the posterior draws.
# - diagnostics: A resolved diagnostics configuration.
# - constant_pars: Names of parameters with constant priors.
# - check_rhat: Logical. Whether Rhat should be evaluated (FALSE if nchain < 2).
warn_nonconvergence = function(convergence, diagnostics = list(), constant_pars = character(), check_rhat = TRUE) {
  diagnostics = resolve_diagnostics(diagnostics)
  if (is.null(unlist(diagnostics[c("rhat", "ess_bulk", "ess_tail")])))
    return(invisible(NULL))

  convergence_check = check_convergence(convergence, diagnostics, constant_pars = constant_pars, check_rhat = check_rhat)
  if (convergence_check$n_bad == 0)
    return(invisible(NULL))

  has_group = any(convergence_check$bad & grepl("[", convergence$parameter, fixed = TRUE))
  warning(
    "Some parameters may not have converged well:\n",
    "  * ", convergence_check$thresholds, ": ",
    format_bad_pars(convergence$parameter, convergence_check$bad), "\n",
    "Inspect `summary(fit)`", if (has_group) ", `ranef(fit)`,", " and `plot_pars(fit)`, and consider increasing ",
    "`iter`/`warmup`, set inits manually, or simplifying the model.",
    call. = FALSE
  )
  invisible(NULL)
}


# Design matrix with one column per coefficient in `predictors` and `group_predictors`,
# ordered by `matrix_col`. The columns are in their `matrix_data`.
get_predictor_matrix = function(predictors, group_predictors = NULL) {
  columns = c(predictors$matrix_data, group_predictors$matrix_data)
  do.call(cbind, columns[order(c(predictors$matrix_col, group_predictors$matrix_col))])
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


# Exact quantiles of the posterior predictive distribution of each row, i.e., the mixture of
# the response distributions of all draws. Returns a rows x probs matrix.
# - dpars: Distributional parameters on the response scale (rows x draws).
# - data: Response auxiliaries, e.g., `trials`, one value per row.
get_mixture_quantiles = function(dpars, data, probs, family, rate = FALSE) {
  quantiles = vapply(seq_len(nrow(dpars$mu)), function(row) {
    row_dpars = lapply(dpars, function(value) value[row, ])
    row_data = lapply(data, `[`, row)
    row_quantiles = vapply(probs, function(p) find_mixture_quantile(family, row_dpars, row_data, p), numeric(1))
    family$response$observed(row_quantiles, row_data, rate)
  }, numeric(length(probs)))
  matrix(quantiles, ncol = length(probs), byrow = TRUE)
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
