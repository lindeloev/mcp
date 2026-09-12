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
#' * `simulate`: A function to simulate data from supplied parameter values.
#' * `.internal`: Information used internally by mcp.
NULL


# Internal function for summary.mcpfit, fixef.mcpfit, and ranef.mcpfit
#
# - fit: An \code{\link{mcpfit}}` object.
# - scope: Which parameter scope to summarise: population-level parameters
#   or group-level deviations.
# - role: Optional parameter role to select within `scope`.
# - verbose: Logical. Include the `segment` and `dpar` columns.
# Returns: A data.frame with summaries for each model parameter. With
#   `verbose = TRUE`, rows are labeled with `segment` and `dpar` columns (see
#   `summary.mcpfit`).
get_summary = function(fit, width, scope = c("population", "group"), role = NULL,
                       dpar = NULL, prior = FALSE, verbose = FALSE) {
  # Check arguments
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_number(width, lower = 0, upper = 1)
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
  selected_names = pars$name[selected]

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

  # A model such as `y ~ 0` has no primary-response fixed effects. Return a
  # regular empty summary instead of asking posterior to summarise no draws.
  if (length(get_cols) == 0) {
    estimates = data.frame(
      variable = character(), mean = numeric(), sd = numeric(), lower = numeric(),
      upper = numeric(), rhat = numeric(), ess_bulk = numeric(), ess_tail = numeric()
    )
    if (verbose) {
      estimates$segment = integer()
      estimates$dpar = character()
    }
    if (!is.null(attr(fit$data[, mcp_columns(fit)$response], "simulated")) &&
        !has_custom_jags_code(fit)) {
      estimates$sim = numeric()
      estimates$match = character()
    }
    return(estimates)
  }

  draws = posterior::subset_draws(draws, variable = get_cols)

  # Get parameter estimates and diagnostics
  tail_prob = (1 - width) / 2
  estimates = posterior::summarise_draws(
    draws,
    mean = base::mean,
    sd = stats::sd,
    lower = function(x) stats::quantile(x, tail_prob, names = FALSE),
    upper = function(x) stats::quantile(x, 1 - tail_prob, names = FALSE),
    rhat = posterior::rhat,
    ess_bulk = function(x) suppressWarnings(posterior::ess_bulk(x)),
    ess_tail = function(x) suppressWarnings(posterior::ess_tail(x))
  ) %>%
    dplyr::mutate(
      ess_bulk = round(.data$ess_bulk),
      ess_tail = round(.data$ess_tail)
  )

  # Order rows and add `segment`/`dpar` using the canonical parameter table
  # built in mcp(). Group-level columns (e.g. "cp_1_id[A]") are matched by
  # their base name; ties (i.e., levels of the same group-level effect) are
  # broken alphabetically by the full column name.
  base_name = sub("\\[.*\\]$", "", estimates$variable)
  match_idx = match(base_name, pars$name)
  estimates$segment = pars$segment[match_idx]
  estimates$dpar = pars$dpar[match_idx]
  estimates = estimates[order(match_idx, estimates$variable), ]

  # Add simulation parameters if the data is simulated
  sim_list = attr(fit$data[, mcp_columns(fit)$response], "simulated")
  if (has_custom_jags_code(fit))
    sim_list = NULL
  if(!is.null(sim_list)) {
    simulated = as.list(sim_list)  # Get as oroper list
    simulated = simulated[sapply(simulated, is.numeric)]  # Remove non-numeric

    # Handle group-level deviations. Find the matching labels.
    for (this_group_effect in mcp_pars(fit, scope = "group")$name) {
      if (!is.null(simulated[[this_group_effect]])) {
        # Find the needed values and labels
        value = simulated[[this_group_effect]]  # Extract simulation values
        group_effects = get_fit_model_tables(fit)$group_effects
        label_col = group_effects$group_col[group_effects$name == this_group_effect]
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

    # Now unpack the whole bunch to a left_join() friendly data.frame.
    simulated = unlist(simulated)  # as named vector
    simulated = data.frame(
      variable = names(simulated),
      sim = as.numeric(simulated),  # without row names
      stringsAsFactors = FALSE
    )

    # Add simulation values for comparison with the fitted parameters.
    estimates = estimates %>%
      dplyr::left_join(simulated, by = "variable", relationship = "one-to-one") %>%
      dplyr::mutate(
        cp_width = ifelse(stringr::str_detect(.data$variable, "^cp_[0-9]+"), .data$upper - .data$lower, 0),
        match = ifelse(
          .data$sim >= (.data$lower - 0.05 * .data$cp_width) &
          .data$sim <= (.data$upper + 0.05 * .data$cp_width),
          yes = "OK", no = ""
        )
      ) %>%
      dplyr::select(-"cp_width")
  }

  # Return-columns and column-order
  if (!verbose)
    estimates = dplyr::select(estimates, -"segment", -"dpar")

  estimates = dplyr::select(
    estimates,
    dplyr::any_of(c(
      "variable", "mean", "sd", "lower", "upper", "rhat", "ess_bulk", "ess_tail",
      "segment", "dpar", "sim", "match"
    ))
  )

  data.frame(estimates, row.names = NULL)
}


#' Summarise mcpfit objects
#'
#' Summarise parameter estimates and model diagnostics.
#'
#' @aliases summary summary.mcpfit
#' @param object An \code{\link{mcpfit}} object.
#' @param width Float. The width of the central posterior interval (between 0
#'   and 1).
#' @param digits Non-negative integer. Number of significant digits used when
#'   printing the summary. Defaults to 2. The invisibly returned data frame
#'   retains the unrounded values.
#' @param prior Logical. Summarise prior draws (`TRUE`) instead of posterior draws (`FALSE`, default)?
#' @param verbose Logical. Include the `segment` and `dpar` columns. Defaults
#'   to `FALSE` for a compact, v0.3.4-compatible summary.
#' @inheritParams mcp
#' @param ... Must be empty. Reserved for future use.
#'
#' @return A data frame with parameter estimates and MCMC diagnostics. Rows
#'   are ordered by change point first, then `mu`, then the other
#'   distributional parameters, then `ar`/`ma` components - each ascending by
#'   segment. OBS: The change point distributions are often not unimodal and
#'   symmetric so the intervals can be deceiving. Plot them using
#'   `plot_pars(fit)`.
#'
#'   With `verbose = TRUE`:
#'   * `segment` is the segment the parameter belongs to.
#'   * `dpar` is the distributional parameter (`"cp"`, `"mu"`, `"sigma"`,
#'     `"ar"`, `"ma"`, etc.) the parameter belongs to. For AR/MA terms, the
#'     lag order is encoded in `variable`, e.g. `ar2_1`.
#'   * `mean` is the posterior mean
#'   * `sd` is the posterior standard deviation.
#'   * `lower` and `upper` are the bounds of the central posterior interval
#'     given in `width`.
#'   * `rhat` is the rank-normalized split-Rhat convergence diagnostic.
#'   * `ess_bulk` and `ess_tail` are the bulk and tail effective sample sizes.
#'     Low effective sample sizes are also obvious as poor mixing in trace plots
#'     (see `plot_pars(fit)`). Read how to deal with such problems [here](https://lindeloev.github.io/mcp/articles/tips.html)
#'
#'  Group-level change-point deviations (`cp_i_id`) follow a standard hierarchical
#'  normal distribution around the population change point. Their realized
#'  locations are truncated to remain in range and ordered.
#'  Predictor group-level effects (such as `Intercept_1_id`) also use standard
#'  hierarchical zero-mean priors, without change-point constraints.
#'
#'  For simulated data, the summary contains two additional columns so that it
#'  is easy to inspect whether the model can recover the parameters. Run
#'  simulation and summary multiple times to get a sense of the robustness.
#'
#'   * `sim` is the value used to generate the data.
#'   * `match` is `"OK"` if `sim` is contained in the central posterior
#'     interval (`lower` to `upper`).
#'
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' # Typical usage
#' summary(demo_fit)
#' summary(demo_fit, width = 0.8, digits = 4)  # Set interval width
#'
#' # Get the results as a data frame
#' results = summary(demo_fit)
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


# Shared display, including unsampled models when printing.
summarize_mcpfit = function(object, width = 0.95, digits = 2, prior = FALSE, verbose = FALSE, diagnostics = NULL, ...) {
  fit = object  # Standard name in mcp
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_number(width, lower = 0, upper = 1)
  checkmate::assert_int(digits, lower = 0)
  checkmate::assert_flag(prior)
  checkmate::assert_flag(verbose)
  rlang::check_dots_empty()

  if (is.null(diagnostics)) {
    diagnostics = .subset2(fit, ".internal")[["diagnostics"]]
    if (is.null(diagnostics))
      diagnostics = list()
  }
  diagnostics = resolve_diagnostics(diagnostics)

  draws = mcmclist_draws(fit, prior = prior, error = FALSE)

  # Model info
  cat(format(fit$family), "\n", sep = "")
  if (!is.null(draws))
    cat("Iterations: ", coda::niter(draws), " from ", coda::nchain(draws), " chains.\n", sep="")
  cat("Segments:\n")
  for (i in seq_along(fit$model)) {
    cat("  ", i, ": ", formula_to_char(fit$model[[i]]), "\n", sep = "")
  }

  # Data
  if (!is.null(draws)) {
    # Print and return population-level summaries invisibly.
    result = get_summary(fit, width, scope = "population", prior = prior, verbose = verbose)
    pars = mcp_pars(fit)
    cp_names = pars$name[pars$part == "cp" & pars$scope == "population"]
    is_cp = result$variable %in% cp_names

    # Format before splitting, so both printed tables share column widths.
    display = data.frame(lapply(result, format, digits = digits), check.names = FALSE)
    if ("rhat" %in% names(display))
      display$rhat = format(result$rhat, digits = digits, nsmall = digits)
    result_cp = display[is_cp, , drop = FALSE]
    result_population = display[!is_cp, , drop = FALSE]

    if (nrow(result_cp) > 0) {
      cat("\nChange point parameters:\n")
      print(data.frame(result_cp), digits = digits, row.names = FALSE)
    }
    if (nrow(result_population) > 0) {
      cat("\nPopulation-level parameters:\n")
      print(data.frame(result_population), digits = digits, row.names = FALSE)
    }

    # Convergence warning footer
    all_res = result
    group_names = mcp_pars(fit, scope = "group")$name
    if (length(group_names) > 0) {
      ran_res = get_summary(fit, width, scope = "group", prior = prior, verbose = verbose)
      cat(
        "\nGroup-level effects: ", paste(group_names, collapse = ", "),
        ". Use `ranef(fit)` to inspect deviations by level.\n", sep = ""
      )
      all_res = dplyr::bind_rows(all_res, ran_res)
    }
    bad_mask = rep(FALSE, nrow(all_res))
    if (!is.null(diagnostics$rhat))
      bad_mask = bad_mask | (!is.na(all_res$rhat) & all_res$rhat > diagnostics$rhat)
    if (!is.null(diagnostics$ess_bulk))
      bad_mask = bad_mask | (!is.na(all_res$ess_bulk) & all_res$ess_bulk < diagnostics$ess_bulk)
    if (!is.null(diagnostics$ess_tail))
      bad_mask = bad_mask | (!is.na(all_res$ess_tail) & all_res$ess_tail < diagnostics$ess_tail)
    n_bad = sum(bad_mask)
    if (n_bad > 0) {
      param_str = if (n_bad == 1) "1 parameter shows" else paste0(n_bad, " parameters show")
      thresholds = c(
        if (!is.null(diagnostics$rhat)) paste0("rhat > ", diagnostics$rhat),
        if (!is.null(diagnostics$ess_bulk)) paste0("ess_bulk < ", diagnostics$ess_bulk),
        if (!is.null(diagnostics$ess_tail)) paste0("ess_tail < ", diagnostics$ess_tail)
      )
      thresholds = paste(thresholds, collapse = " or ")
      cat("\nWarning: ", param_str, " poor convergence (", thresholds, ").\n", sep = "")
    }

    return(invisible(result))
  }
  else {
    cat("\nNo draws. Nothing to summarise.")
    return(invisible(NULL))
  }
}



#' @aliases fixef fixef.mcpfit
#' @describeIn summary.mcpfit Population-level fixed effects (regression coefficients) of `mcpfit`.
#' @param dpar Distributional parameter(s) whose regression coefficients to
#'   return. For modeled distributional parameters such as `sigma()`, these
#'   coefficients are on the link scale.
#' @export fixef
#' @exportS3Method nlme::fixef
fixef.mcpfit = function(object, width = 0.95, prior = FALSE, verbose = FALSE, dpar = "mu", ...) {
  rlang::check_dots_empty()
  checkmate::assert_subset(dpar, object$family$dpar_specs$dpar)
  get_summary(
    object, width, scope = "population", role = c("fixed_effect", "dpar_effect"),
    dpar = dpar, prior = prior, verbose = verbose
  )
}

#' @aliases ranef ranef.mcpfit
#' @describeIn summary.mcpfit Group-level deviations (random effects) of `mcpfit`.
#'   Change-point deviations are relative to their population change point;
#'   `cp_i_sd` is the scale of their latent normal distribution.
#' @export ranef
#' @exportS3Method nlme::ranef
ranef.mcpfit = function(object, width = 0.95, prior = FALSE, verbose = FALSE, ...) {
  rlang::check_dots_empty()
  get_summary(object, width, scope = "group", prior = prior, verbose = verbose)
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
  rlang::check_dots_empty()
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
#'   `"all"` for all population-level parameters.
#' @param dpar Distributional parameter(s) to select when `pars = NULL`.
#' @param parm Optional names or positions of population-level parameters to
#'   include in the intervals.
#' @param level Width of the central interval.
#' @param prior Logical. Use prior draws instead of posterior draws?
#' @param ... Must be empty. Reserved for future use.
#' @return `vcov()` returns a covariance or correlation matrix. `confint()`
#'   returns a two-column matrix of central intervals.
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
#' vcov(demo_fit, pars = "all")
#'
#' # Inspect posterior parameter correlations across the full population model.
#' # Useful to quickly check identifiability (high correlation). Inspecting
#' # `bayesplot::mcmc_pairs(as_draws(demo_fit))` is better, though.
#' vcov(demo_fit, pars = "all", correlation = TRUE)
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
    checkmate::assert_subset(dpar, object$family$dpar_specs$dpar)
    pars = parameters$name[
      parameters$role %in% c("fixed_effect", "dpar_effect") &
        parameters$dpar %in% dpar
    ]
  } else if (identical(pars, "all")) {
    pars = parameters$name
  } else {
    checkmate::assert_character(pars, any.missing = FALSE)
    checkmate::assert_subset(pars, parameters$name)
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
  checkmate::assert_number(level, lower = 0, upper = 1)
  checkmate::assert_true(level > 0 && level < 1, .var.name = "level")
  checkmate::assert_flag(prior)

  population = mcp_pars(object, scope = "population")$name
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


#' @aliases print print.mcpfit
#' @describeIn summary.mcpfit Print the posterior summary of an \code{\link{mcpfit}} object.
#' @param x An \code{\link{mcpfit}} object.
#' @export
print.mcpfit = function(x, ...) {
  args = list(...)
  if (isTRUE(args$prior))
    mcmclist_draws(x, prior = TRUE)
  if (!isTRUE(args$prior) && !coda::is.mcmc.list(.subset2(x, "mcmc_post")) &&
      coda::is.mcmc.list(.subset2(x, "mcmc_prior"))) {
    message("Posterior was not drawn. Using prior draws. Set `prior = TRUE` to mute this message.")
    args$prior = TRUE
  }
  do.call(summarize_mcpfit, c(list(object = x), args))
}


#' Checks if the Argument is an `mcpfit` Object
#'
#' @aliases is.mcpfit
#' @param x An `R` object.
#' @return Logical scalar (`TRUE` if `x` is an `mcpfit` object, `FALSE` otherwise).
#' @export
is.mcpfit = function(x) {
  inherits(x, "mcpfit")
}


#' @export
`$.mcpfit` = function(x, name) {
  if (name == "fit" && is.null(.subset2(x, "fit"))) {
    warning("`mcp_example()` now returns an `mcpfit` object directly instead of a list with `$fit`. Returning the object itself.", call. = FALSE)
    return(x)
  }
  if (name %in% c("mcmc_post", "mcmc_prior")) {
    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = I(paste0("fit$", name)),
      with = I("as_draws(fit) or coda::as.mcmc(fit)")
    )
  }
  if (name == "pars") {
    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = I("fit$pars"),
      with = I("mcp_pars(fit) and mcp_columns(fit)")
    )
  }
  if (name == "log_lik" && !name %in% names(x)) {
    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = I("fit$log_lik"),
      with = I("log_lik(fit)")
    )
  }
  .subset2(x, name)
}

#' @export
`[[.mcpfit` = function(x, i, exact = TRUE) {
  if (is.character(i) && i %in% c("mcmc_post", "mcmc_prior")) {
    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = I(paste0("fit$", i)),
      with = I("as_draws(fit) or coda::as.mcmc(fit)")
    )
  }
  if (is.character(i) && identical(i, "pars")) {
    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = I("fit$pars"),
      with = I("mcp_pars(fit) and mcp_columns(fit)")
    )
  }
  .subset2(x, i)
}
