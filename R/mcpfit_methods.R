# ABOUT: These are non-plotting functions that take an mcpfit as the first argument
# -----------------

#' Class `mcpfit` of Models Fitted with the \pkg{mcp} Package
#'
#' Models fitted with the \code{\link[mcp:mcp]{mcp}} function are represented as
#' an `mcpfit` object which contains the user input (model, data, family),
#' derived model characteristics (prior, parameter names, and jags code), and
#' the fit (prior and/or posterior MCMC draws).
#'
#' @name mcpfit-class
#' @aliases mcpfit
#' @docType class
#'
#' @details
#' See `methods(class = "mcpfit")` for an overview of available methods.
#'
#' Components:
#' * `call`: The matched call to `mcp()`.
#' * `model`: A list of user-provided formulas.
#' * `data`: The user-provided data frame reduced to model-used columns.
#' * `family`: An `mcpfamily` object.
#' * `prior`: A named list of priors.
#' * `mcmc_post` and `mcmc_prior`: \code{\link[coda]{mcmc.list}} objects with
#'   posterior and prior draws, respectively. Do not access these directly; 
#'   use as_draws(fit) or similar.
#' * `jags_code`: A string with JAGS code; use `cat(fit$jags_code)` to show it.
#'   [Read more](https://lindeloev.github.io/mcp/articles/implementation.html).
#' * `simulate`: A function to simulate data from supplied parameter values.
#' * `version`: A named list of [package_version]s for `mcp`, `R`, `rjags`,
#'   and `jags`, e.g., `fit$version$mcp >= "0.4.0"`. `rjags` and `jags` are
#'   `NULL` for fits with `sample = "none"`. `NULL` for fits made with mcp < 0.4.0.
#' * `.internal`: Information used internally by mcp.
NULL


# Parameter values used to simulate the response with `fit$simulate()`, or NULL. Only for fits of the
# simulating model, i.e., with the same parameters, so that other models are not compared to them.
get_simulated = function(fit) {
  simulated = as.list(attr(fit$data[, mcp_columns(fit)$response], "simulated"))
  simulated = simulated[vapply(simulated, is.numeric, logical(1))]  # Remove non-numeric
  if (has_custom_jags_code(fit) || !setequal(names(simulated), mcp_pars(fit)$parameter))
    return(NULL)
  simulated
}


# Simulated values named like the draws, e.g., "cp_1" and "cp_1_id[a]", or NULL
get_simulated_values = function(fit) {
  simulated = get_simulated(fit)
  if (is.null(simulated))
    return(NULL)

  # Use one value for each level of group-level deviations
  group_pars = mcp_pars(fit, scope = "group")
  for (this_group_effect in group_pars$parameter) {
    if (!is.null(simulated[[this_group_effect]])) {
      # Find the needed values and labels
      value = simulated[[this_group_effect]]  # Extract simulation values
      label_col = group_pars$group_col[group_pars$parameter == this_group_effect]
      labs = fit$data[[label_col]]  # Find the labels. Same length as `value`
      if (length(value) != length(labs)) {
        warning("This is simulated data, but the labels for group-level effect '", label_col, "' in data do not have the same length as the numeric parameters used for simulation.")
        next
      }

      # Name like the MCMC columns and use one value for each group level.
      keep = !duplicated(labs)
      value = value[keep]
      names(value) = paste0(this_group_effect, "[", labs[keep], "]")

      # Delete the simulation vector and add the new label-value pairs to list
      simulated[[this_group_effect]] = NULL
      simulated = c(simulated, as.list(value))
    }
  }
  unlist(simulated)
}


# Add the `sim` values and whether they are within the central interval (`match`)
add_sim_match = function(estimates, sim, fit) {
  estimates$sim = sim
  is_match = estimates$sim >= estimates$lower & estimates$sim <= estimates$upper

  # A disjoined change point has the same likelihood anywhere between two adjacent observed x.
  # Match if `sim` is between the same observed x as the interval.
  cols = mcp_columns(fit)
  x_obs = sort(unique(fit$data[[cols$par_x]][!is.na(fit$data[[cols$response]])]))
  gap = function(value) findInterval(value, x_obs, left.open = TRUE)
  is_disjoined = estimates$parameter %in% get_disjoined_cps(fit)
  is_match[is_disjoined] = (gap(estimates$lower) <= gap(estimates$sim) & gap(estimates$sim) <= gap(estimates$upper))[is_disjoined]
  estimates$match = ifelse(is_match, yes = "OK", no = "")
  estimates
}


# Parameter estimates and diagnostics of each variable in `draws`, as in summary.mcpfit.
# - convergence: A convergence_table() covering `draws`. Computed if NULL.
summarise_pars = function(draws, width, convergence = NULL) {
  tail_prob = (1 - width) / 2
  if (is.null(convergence))
    convergence = convergence_table(draws)

  posterior::summarise_draws(
    draws,
    mean = base::mean,
    sd = stats::sd,
    lower = function(x) stats::quantile(x, tail_prob, names = FALSE),
    upper = function(x) stats::quantile(x, 1 - tail_prob, names = FALSE)
  ) %>%
    dplyr::rename(parameter = "variable") %>%
    dplyr::left_join(convergence, by = "parameter") %>%
    dplyr::mutate(
      ess_bulk = round(.data$ess_bulk),
      ess_tail = round(.data$ess_tail)
  )
}


# Internal function for summary.mcpfit, fixef.mcpfit, and ranef.mcpfit
#
# - fit: An \code{\link{mcpfit}}` object.
# - scope: Which parameter scope to summarise: population-level parameters
#   or group-level deviations.
# - role: Optional parameter role to select within `scope`.
# - verbose: Logical. Include the `segment` and `dpar` columns.
# Returns: A tibble with summaries for each model parameter. With
#   `verbose = TRUE`, rows are labeled with `segment` and `dpar` columns (see
#   `summary.mcpfit`).
get_summary = function(fit, width, scope = c("population", "group"), role = NULL,
                       dpar = NULL, prior = FALSE, verbose = FALSE) {
  # Check arguments
  checkmate::assert_class(fit, "mcpfit")
  assert_width(width)
  scope = rlang::arg_match0(scope, c("population", "group"))
  checkmate::assert_character(role, null.ok = TRUE)
  checkmate::assert_character(dpar, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_flag(prior)
  checkmate::assert_flag(verbose)

  draws = posterior_draws(fit, prior = prior)
  if (scope == "group" && nrow(mcp_pars(fit, scope = "group")) == 0)
    return(NULL)

  # Select by the independent scope and role dimensions of the parameter table.
  all_cols = posterior::variables(draws)
  pars = mcp_pars(fit)
  selected = pars$scope == scope
  if (!is.null(role))
    selected = selected & pars$role %in% role
  if (!is.null(dpar))
    selected = selected & pars$dpar %in% dpar
  selected_names = pars$parameter[selected]

  if (scope == "population") {
    get_cols = all_cols[all_cols %in% selected_names]
  } else {
    get_cols = all_cols[vapply(
      all_cols,
      function(column) any(startsWith(column, paste0(selected_names, "["))),
      logical(1)
    )]
    if (length(get_cols) == 0)
      stop("There were no matching parameters in the model.")
  }

  # The selection can be empty, e.g., `fixef(fit, dpar = "mu")` when all mu segments are `~ 0`.
  # Return a regular empty summary instead of asking posterior to summarise no draws.
  if (length(get_cols) == 0) {
    estimates = tibble::tibble(
      parameter = character(), mean = numeric(), sd = numeric(), lower = numeric(),
      upper = numeric(), rhat = numeric(), ess_bulk = numeric(), ess_tail = numeric()
    )
    if (verbose) {
      estimates$segment = integer()
      estimates$dpar = character()
    }
    if (!is.null(get_simulated(fit))) {
      estimates$sim = numeric()
      estimates$match = character()
    }
    return(estimates)
  }

  draws = posterior::subset_draws(draws, variable = get_cols)
  estimates = summarise_pars(draws, width, convergence = if (!prior) cached_convergence(fit))

  # Order rows and add `segment`/`dpar` using the canonical parameter table
  # built in mcp(). Group-level columns (e.g. "cp_1_id[A]") are matched by
  # their base name; ties (levels of the same group-level effect) preserve their
  # incoming factor-level order from draws.
  base_name = sub("\\[.*\\]$", "", estimates$parameter)
  match_idx = match(base_name, pars$parameter)
  estimates$segment = pars$segment[match_idx]
  estimates$dpar = pars$dpar[match_idx]
  estimates = estimates[order(match_idx), ]

  # Add simulation parameters if the data is simulated
  simulated = get_simulated_values(fit)
  if (!is.null(simulated))
    estimates = add_sim_match(estimates, unname(simulated[estimates$parameter]), fit)

  # Return-columns and column-order
  if (!verbose)
    estimates = dplyr::select(estimates, -"segment", -"dpar")

  estimates = dplyr::select(
    estimates,
    dplyr::any_of(c(
      "parameter", "mean", "sd", "lower", "upper", "rhat", "ess_bulk", "ess_tail",
      "segment", "dpar", "sim", "match"
    ))
  )

  tibble::as_tibble(estimates)
}


#' Summarise mcpfit objects
#'
#' Summarise parameter estimates and model diagnostics.
#'
#' @aliases summary summary.mcpfit
#' @param object An \code{\link{mcpfit}} object.
#' @param width Float. The width of the central posterior interval (between 0
#'   and 1).
#' @param digits Positive integer. Number of significant digits of `sd` when
#'   printing the summary. Defaults to 2. `mean`, `lower`, `upper`, and `sim` are
#'   rounded to the same decimal places in each row. The returned tibble retains
#'   the unrounded values.
#' @param prior Logical. Summarise prior draws (`TRUE`) instead of posterior draws (`FALSE`, default)?
#' @param verbose Logical. Include the `segment` and `dpar` columns. Defaults
#'   to `FALSE` for a compact summary.
#' @inheritParams mcp
#' @param ... Must be empty. Reserved for future use.
#'
#' @return An `mcpsummary` tibble with parameter estimates and MCMC diagnostics.
#'   It prints with model information, separate tables for change points and
#'   population-level parameters, and convergence warnings. Rows
#'   are ordered by change point first, then `mu`, then the other
#'   distributional parameters, then AR/MA lags - each ascending by
#'   segment. Note: The change point distributions are often not unimodal and
#'   symmetric so the intervals can be deceiving. Plot them using
#'   `plot_pars(fit)`.
#'
#'   Columns:
#'   * `segment` (only with `verbose = TRUE`) is the segment the parameter belongs to.
#'     For change points, it is the segment that the change point starts, e.g., 2 for `cp_1`.
#'   * `dpar` (only with `verbose = TRUE`) is the distributional parameter (e.g., `"mu"`,
#'     `"sigma"`, or `"ar1"`; see [mcp-formula]) as accepted by the `dpar` argument of
#'     [fixef.mcpfit()], [fitted.mcpfit()], etc. It is `NA` for change points. See also [mcp_pars()].
#'   * `mean` is the posterior mean
#'   * `sd` is the posterior standard deviation.
#'   * `lower` and `upper` are the lower and upper ends of the central posterior interval
#'     given in `width`.
#'   * `rhat` is the rank-normalized split-Rhat convergence diagnostic.
#'   * `ess_bulk` and `ess_tail` are the bulk and tail effective sample sizes.
#'     Low effective sample sizes are also obvious as poor mixing in trace plots
#'     (see `plot_pars(fit)`). Read how to deal with such problems in
#'     [Diagnosing and fixing problems](https://lindeloev.github.io/mcp/articles/tips.html#convergence).
#'
#'  Group-level change-point deviations (`cp_i_id`) follow a hierarchical normal
#'  distribution around the population-level change point before truncation. The
#'  group-specific change points are truncated to the observed range of x and ordered, so the
#'  population-level change point is generally close to their mean, but not exactly.
#'  Predictor group-level deviations (such as `Intercept_1_id`) use standard
#'  hierarchical zero-mean priors, without these bounds.
#'  [Read more](https://lindeloev.github.io/mcp/articles/group_effects.html).
#'
#'  For simulated data, the summary contains two additional columns so that it
#'  is easy to inspect whether the model can recover the parameters. Run
#'  simulation and summary multiple times to get a sense of the robustness.
#'
#'   * `sim` is the value used to generate the data.
#'   * `match` is `"OK"` if `sim` is contained in the central posterior
#'     interval (`lower` to `upper`). For change points starting a disjoined
#'     segment, `sim` only needs to be between the same two adjacent observed
#'     `x` as the interval, since the likelihood is identical anywhere between them.
#'     Group-level deviations in `ranef()` are relative to their population-level parameter,
#'     so they share its error; `coef()` checks the group-specific values directly.
#'
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' # Typical usage
#' summary(demo_fit)
#' summary(demo_fit, width = 0.8, digits = 4)  # Set interval width
#'
#' # Get the results as a tibble
#' results = summary(demo_fit)
#' results$mean
#'
#' # Group-level deviations (random effects)
#' # ranef(my_fit)
#'
#' # Summarise prior
#' summary(demo_fit, prior = TRUE)
summary.mcpfit = function(object, width = 0.95, digits = 2, prior = FALSE, verbose = FALSE, diagnostics = NULL, ...) {
  mcmclist_draws(object, prior = prior)
  summarize_mcpfit(object, width, digits, prior, verbose, diagnostics, ...)
}


# Shared by summary() and print(), including unsampled models when printing.
# Returns the population-level summary as an `mcpsummary` tibble. Its "mcpsummary"
# attribute holds what print.mcpsummary() shows around the tables.
summarize_mcpfit = function(object, width = 0.95, digits = 2, prior = FALSE, verbose = FALSE, diagnostics = NULL, ...) {
  fit = object  # Standard name in mcp
  checkmate::assert_class(fit, "mcpfit")
  assert_width(width)
  checkmate::assert_int(digits, lower = 1)
  checkmate::assert_flag(prior)
  checkmate::assert_flag(verbose)
  rlang::check_dots_empty()

  if (is.null(diagnostics))
    diagnostics = fit$.internal$diagnostics
  diagnostics = resolve_diagnostics(diagnostics)

  draws = mcmclist_draws(fit, prior = prior, error = FALSE)

  # Model info
  forms = vapply(fit$model, formula_to_char, character(1))
  header = c(
    format(fit$family),
    if (!is.null(draws)) paste0("Draws: ", coda::niter(draws) * coda::nchain(draws), " (", coda::nchain(draws), " chains with iter = ", coda::niter(draws), ")"),
    "Segments:",
    paste0("  ", format(seq_along(forms)), ": ", forms)
  )

  # Data
  if (!is.null(draws)) {
    result = get_summary(fit, width, scope = "population", prior = prior, verbose = verbose)
    pars = mcp_pars(fit)
    cp_names = pars$parameter[pars$part == "cp" & pars$scope == "population"]

    # Convergence warning footer. Counts group-level deviations by level, which are summarised in ranef().
    footer = character()
    group_names = mcp_pars(fit, scope = "group")$parameter
    if (length(group_names) > 0) {
      footer = paste0(
        "Group-level effects: ", paste(group_names, collapse = ", "),
        ". Use `ranef(fit)` to inspect deviations by level."
      )
    }
    convergence = if (!prior) cached_convergence(fit)
    if (is.null(convergence))
      convergence = convergence_table(draws)
    check_rhat = !is.null(diagnostics$rhat) && coda::nchain(draws) >= 2
    prior_table = fit$.internal$prior_table
    constant_pars = prior_table$parameter[prior_table$kind == "constant"]
    convergence_check = check_convergence(convergence, diagnostics, constant_pars = constant_pars, check_rhat = check_rhat)
    is_group = grepl("[", convergence$parameter, fixed = TRUE)
    n_bad_population = sum(convergence_check$bad & !is_group)
    if (n_bad_population > 0) {
      param_str = if (n_bad_population == 1) "1 parameter shows" else paste0(n_bad_population, " parameters show")
      footer = c(footer, paste0("Warning: ", param_str, " poor convergence (", convergence_check$thresholds, ")."))
    }
    if (any(convergence_check$bad & is_group)) {
      footer = c(footer, paste0(
        "Warning: Group-level deviations show poor convergence (", convergence_check$thresholds, "): ",
        format_bad_pars(convergence$parameter, convergence_check$bad & is_group), "."
      ))
    }
  } else {
    result = tibble::tibble()
    cp_names = character()
    footer = if (!prior && !is.null(mcmclist_draws(fit, prior = TRUE, error = FALSE))) {
      "No posterior draws. Use `summary(fit, prior = TRUE)` to summarise the prior draws."
    } else {
      paste0("No ", if (prior) "prior" else "posterior", " draws. Nothing to summarise.")
    }
  }

  attr(result, "mcpsummary") = list(header = header, cp_names = cp_names, footer = footer, digits = digits)
  class(result) = c("mcpsummary", class(result))
  result
}


#' @aliases print.mcpsummary
#' @describeIn summary.mcpfit Print an `mcpsummary` returned by `summary()`: model
#'   information, change point and population-level parameters, and diagnostics.
#' @export
print.mcpsummary = function(x, digits = attr(x, "mcpsummary")$digits, ...) {
  info = attr(x, "mcpsummary")
  if (is.null(info) || (nrow(x) > 0 && "parameter" %notin% names(x)))
    return(NextMethod())  # e.g., after dplyr::select() dropped `parameter`

  cat(info$header, sep = "\n")
  if (nrow(x) > 0) {
    # Format before splitting, so both printed tables share column widths.
    display = data.frame(lapply(x, format, digits = digits), check.names = FALSE)

    # Round each row to the decimals that show `digits` significant digits of its sd.
    # Pad right to align decimal points (and add the "." missing for 0 decimals).
    decimals = pmax(0, digits - 1 - floor(log10(abs(x$sd))))
    decimals[!is.finite(decimals)] = digits  # E.g., sd = 0 or NA
    pad = strrep(" ", max(decimals) - decimals + (decimals == 0 & max(decimals) > 0))
    for (col in intersect(c("mean", "sd", "lower", "upper", "sim"), names(x)))
      display[[col]] = format(paste0(sprintf("%.*f", decimals, x[[col]]), pad), justify = "right")
    
    if ("rhat" %in% names(display))
      display$rhat = sprintf("%.2f", x$rhat)
    is_cp = x$parameter %in% info$cp_names

    # Start printing
    if (any(is_cp)) {
      cat("\nChange point parameters:\n")
      print(display[is_cp, , drop = FALSE], digits = digits, row.names = FALSE)
    }
    if (any(!is_cp)) {
      cat("\nPopulation-level parameters:\n")
      print(display[!is_cp, , drop = FALSE], digits = digits, row.names = FALSE)
    }
  }
  if (length(info$footer) > 0)
    cat(paste0("\n", info$footer), sep = "\n")
  invisible(x)
}



#' @aliases fixef fixef.mcpfit
#' @describeIn summary.mcpfit Population-level coefficients (fixed effects) of `mcpfit`.
#' @param dpar Distributional parameter(s) whose regression coefficients to
#'   return, e.g., `"mu"`, `"sigma"`, or `"ar1"` (see [mcp-formula]). For modeled
#'   distributional parameters such as `sigma()`, these coefficients are on the link scale.
#' @export fixef
#' @exportS3Method nlme::fixef
fixef.mcpfit = function(object, width = 0.95, prior = FALSE, verbose = FALSE, dpar = "mu", ...) {
  rlang::check_dots_empty()
  checkmate::assert_subset(dpar, c(object$family$dpar_specs$dpar, stats::na.omit(mcp_pars(object)$dpar)))
  get_summary(
    object, width, scope = "population", role = "coefficient",
    dpar = dpar, prior = prior, verbose = verbose
  )
}

#' @aliases ranef ranef.mcpfit
#' @describeIn summary.mcpfit Group-level deviations (random effects) of `mcpfit`.
#'   Change-point deviations are relative to their population-level change point;
#'   `cp_i_<group>_sd` (e.g., `cp_1_id_sd`) is the SD of their normal distribution before truncation.
#'   Use `coef()` for group-specific values.
#' @export ranef
#' @exportS3Method nlme::ranef
ranef.mcpfit = function(object, width = 0.95, prior = FALSE, verbose = FALSE, ...) {
  rlang::check_dots_empty()
  get_summary(object, width, scope = "group", prior = prior, verbose = verbose)
}

#' @describeIn summary.mcpfit Group-specific coefficients and change points of `mcpfit`:
#'   the population-level parameter plus each group-level deviation, as in `lme4` and `brms`.
#'   For example, `cp_1` for `group_col = "id"` and `level = "a"` summarises the draws of
#'   `cp_1 + cp_1_id[a]`, which is the change point of group `a`. Deviations without a
#'   population-level parameter, e.g., from `(0 + x || id)`, are added to zero.
#'   Coefficients of distributional parameters such as `sigma()` are on the link scale.
#'   For simulated data, `sim` is the sum of the simulated values.
#' @export
coef.mcpfit = function(object, width = 0.95, prior = FALSE, verbose = FALSE, ...) {
  rlang::check_dots_empty()
  assert_width(width)
  group_pars = mcp_pars(object, scope = "group")
  if (nrow(group_pars) == 0)
    stop("The model has no group-level effects. Use `fixef(fit)` for population-level coefficients.", call. = FALSE)

  # Sum population-level and group-level draws for each group level
  draws = posterior::as_draws_df(posterior_draws(object, prior = prior))
  columns = posterior::variables(draws)
  sums = list()
  rows = list()
  for (i in seq_len(nrow(group_pars))) {
    population = if (is.na(group_pars$population_parameter[i])) 0 else draws[[group_pars$population_parameter[i]]]
    group_columns = columns[startsWith(columns, paste0(group_pars$parameter[i], "["))]
    sums[group_columns] = lapply(group_columns, function(column) population + draws[[column]])
    rows[[i]] = tibble::tibble(
      parameter = sub(paste0("_", group_pars$group_col[i], "$"), "", group_pars$parameter[i]),
      group_col = group_pars$group_col[i],
      level = sub("^.*?\\[(.*)\\]$", "\\1", group_columns),
      segment = group_pars$segment[i],
      dpar = group_pars$dpar[i],
      .population = group_pars$population_parameter[i],
      .column = group_columns
    )
  }
  rows = dplyr::bind_rows(rows)
  sums = tibble::as_tibble(sums)
  sums$.chain = draws$.chain
  sums$.iteration = draws$.iteration

  estimates = summarise_pars(posterior::as_draws_df(sums), width)
  estimates = dplyr::bind_cols(dplyr::select(rows, -".population", -".column"), dplyr::select(estimates, -"parameter"))

  # Add simulation parameters if the data is simulated
  simulated = get_simulated_values(object)
  if (!is.null(simulated)) {
    population_sim = ifelse(is.na(rows$.population), 0, simulated[rows$.population])
    estimates = add_sim_match(estimates, unname(population_sim + simulated[rows$.column]), object)
  }
  if (!verbose)
    estimates = dplyr::select(estimates, -"segment", -"dpar")
  estimates
}


#' Extract Model Information from an `mcpfit`
#'
#' Standard R accessors for the model formulas, family, fitting data, and number
#' of observations stored in an `mcpfit`.
#'
#' @param object,x,formula An `mcpfit` object.
#' @param segment `NULL` to return all segment formulas, or a positive integer
#'   selecting one segment.
#' @param ... Must be empty. Reserved for future use.
#' @return `formula()` returns the complete list of segment formulas, or one
#'   formula when `segment` is supplied. `family()` returns an `mcpfamily`.
#'   `model.frame()` returns the data retained in the fit. `nobs()` returns the
#'   number of observed response values (excluding missing responses).
#' @name model-accessors-mcpfit
#' @examples
#' formula(demo_fit)  # Show all segment formulas
#' formula(demo_fit, segment = 2)  # Show the formula for segment 2
#' family(demo_fit)  # Show the response family and link
#' head(model.frame(demo_fit))  # Show the top rows of fitting data
#' nobs(demo_fit)  # Count observed response rows
NULL


#' @rdname model-accessors-mcpfit
#' @export
family.mcpfit = function(object, ...) {
  rlang::check_dots_empty()
  object$family
}


#' @rdname model-accessors-mcpfit
#' @export
nobs.mcpfit = function(object, ...) {
  # No dots check: stats::sigma() etc. call nobs(object, use.fallback = ...)
  y_col = mcp_columns(object)$response
  sum(!is.na(object$data[[y_col]]))
}


#' @rdname model-accessors-mcpfit
#' @export
model.frame.mcpfit = function(formula, ...) {
  rlang::check_dots_empty()
  formula$data
}


#' @rdname model-accessors-mcpfit
#' @export
formula.mcpfit = function(x, segment = NULL, ...) {
  rlang::check_dots_empty()
  checkmate::assert_int(segment, lower = 1, upper = length(x$model), null.ok = TRUE)
  if (is.null(segment))
    return(x$model)

  x$model[[segment]]
}


#' Prior and Posterior Covariance and Central Intervals for `mcpfit` Objects
#'
#' Summarise the joint and marginal uncertainty of population-level model
#' parameters using posterior or prior draws.
#'
#' @param object An `mcpfit` object.
#' @param correlation Return the correlation matrix instead of the covariance
#'   matrix?
#' @param pars Optional names of population-level parameters to extract, or
#'   `"population"` for all population-level parameters.
#' @param dpar Distributional parameter(s) to select when `pars = NULL`, e.g., `"mu"`, `"sigma"`, or `"ar1"`.
#' @param parm Optional names or positions of population-level parameters to
#'   include in the intervals.
#' @param level Width of the central interval.
#' @param prior Logical. Use prior draws instead of posterior draws?
#' @param ... Must be empty. Reserved for future use.
#' @return `vcov()` returns a covariance or correlation matrix. `confint()`
#'   returns a two-column matrix of central intervals.
#' @details As in `lme4`, `vcov()` covers the coefficients of [fixef()][fixef.mcpfit()] (by default those of `mu`),
#'   while `confint()` covers all population-level parameters, including change points.
#' @name posterior-uncertainty-mcpfit
#' @examples
#' # Posterior covariance of the primary-response coefficients, matching fixef().
#' vcov(demo_fit)
#'
#' # Central posterior intervals for all population-level parameters, or a selection.
#' confint(demo_fit)
#' confint(demo_fit, parm = "cp_1")
#' confint(demo_fit, parm = c("Intercept_1", "time_2"), level = 0.8)
#' confint(demo_fit, prior = TRUE)
#'
#' # Include change points, residual SDs, group SDs, and AR/MA parameters.
#' vcov(demo_fit, pars = "population")
#'
#' # Inspect posterior parameter correlations across the full population model.
#' # Useful to quickly check identifiability (high correlation). Inspecting
#' # `bayesplot::mcmc_pairs(as_draws(demo_fit))` is better, though.
#' vcov(demo_fit, pars = "population", correlation = TRUE)
NULL


#' @rdname posterior-uncertainty-mcpfit
#' @export
vcov.mcpfit = function(object, correlation = FALSE, pars = NULL, dpar = "mu",
                       prior = FALSE, ...) {
  rlang::check_dots_empty()
  checkmate::assert_flag(correlation)
  checkmate::assert_flag(prior)
  parameters = mcp_pars(object, scope = "population")
  if (is.null(pars)) {
    checkmate::assert_subset(dpar, c(object$family$dpar_specs$dpar, stats::na.omit(parameters$dpar)))
    pars = parameters$parameter[
      parameters$role == "coefficient" &
        parameters$dpar %in% dpar
    ]
  } else if (identical(pars, "population")) {
    pars = parameters$parameter
  } else {
    checkmate::assert_character(pars, any.missing = FALSE)
    checkmate::assert_subset(pars, parameters$parameter)
  }
  if (length(pars) == 0)
    return(NULL)

  draws = posterior::as_draws_matrix(posterior_draws(
    object, prior = prior
  ))
  if (correlation)
    return(stats::cor(draws[, pars, drop = FALSE]))

  stats::cov(draws[, pars, drop = FALSE])
}


#' @rdname posterior-uncertainty-mcpfit
#' @export
confint.mcpfit = function(object, parm, level = 0.95, prior = FALSE, ...) {
  rlang::check_dots_empty()
  assert_width(level, name = "level")
  checkmate::assert_flag(prior)

  population = mcp_pars(object, scope = "population")$parameter
  if (missing(parm)) {
    parm = population
  } else if (is.numeric(parm)) {
    parm = population[parm]
  } else {
    checkmate::assert_character(parm)
  }
  if (!all(parm %in% population))
    stop("`parm` must name population-level parameters.", call. = FALSE)

  # Compute credible interval
  probs = c((1 - level) / 2, 1 - (1 - level) / 2)
  draws = posterior::as_draws_matrix(posterior_draws(
    object, prior = prior
  ))
  intervals = parm |>
    vapply(
      function(parameter) stats::quantile(draws[, parameter], probs = probs, names = FALSE),
      numeric(2)
    ) |>
    t()

  colnames(intervals) = paste0(format(100 * probs, trim = TRUE), " %")
  intervals
}


# The default method would silently re-evaluate `fit$call` with whatever `model` and `data` are in the caller's environment.
#' @export
update.mcpfit = function(object, ...) {
  stop("`update()` is not implemented for mcpfit objects. Call `mcp()` again instead, e.g., `mcp(fit$model, fit$data, family = family(fit), ...)`.", call. = FALSE)
}


#' @aliases print print.mcpfit
#' @describeIn summary.mcpfit Print the posterior summary of an \code{\link{mcpfit}} object,
#'   or the prior summary with `prior = TRUE`. Without the requested draws, only the model is shown.
#' @param x An \code{\link{mcpfit}} object, or an `mcpsummary` returned by `summary()`.
#' @export
print.mcpfit = function(x, ...) {
  print(summarize_mcpfit(x, ...))
  invisible(x)
}


#' Checks if the Argument is an `mcpfit` Object
#'
#' @aliases is_mcpfit
#' @param x An `R` object.
#' @return Logical scalar (`TRUE` if `x` is an `mcpfit` object, `FALSE` otherwise).
#' @export
is_mcpfit = function(x) {
  inherits(x, "mcpfit")
}
