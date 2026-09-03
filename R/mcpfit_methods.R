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

  if (scope == "group" && nrow(mcp_pars(fit, scope = "group")) == 0)
    return(NULL)

  draws = posterior_draws(fit, prior = prior)

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
    object, prior = prior, fallback_to_prior = FALSE
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
    object, prior = prior, fallback_to_prior = FALSE
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
  summary(x, ...)
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


# Internal function to get draws.
#
# Returns posterior draws, if available. If not, then prior draws. If not,
# then throw an informative error. This is useful for summary and plotting, that
# works on both.
#
# - fit: An \code{\link{mcpfit}} object
# - message: TRUE: gives a message if returning prior draws. FALSE = no message
# - error: TRUE: err if there are no draws. FALSE: return NULL
# - fallback_to_prior: TRUE: use prior draws when posterior draws are unavailable
mcmclist_draws = function(fit, prior = FALSE, message = TRUE, error = TRUE,
                          fallback_to_prior = TRUE) {
  check_mcpfit_version(fit)
  mcmc_prior = .subset2(fit, "mcmc_prior")
  mcmc_post = .subset2(fit, "mcmc_post")
  if (prior == TRUE) {
    if (coda::is.mcmc.list(mcmc_prior)) {
      return(mcmc_prior)
    } else {
      stop("Prior requested but the prior was not drawn.")
    }
  }

  if (coda::is.mcmc.list(mcmc_post)) {
    return(mcmc_post)
  } else if (fallback_to_prior && coda::is.mcmc.list(mcmc_prior)) {
    if (message)
      message("Posterior was not drawn. Using prior draws. Set `prior = TRUE` to mute this message.")
    return(mcmc_prior)
  } else if (error == TRUE) {
    if (fallback_to_prior)
      stop("This mcpfit contains no posterior or prior draws.")
    stop("Posterior requested but the posterior was not drawn.")
  }

  NULL
}


# Get draws as a posterior draws array
#
# This is the single internal conversion from the stored
# \code{\link[coda]{mcmc.list}} representation to a posterior draws object.
posterior_draws = function(fit, prior = FALSE, message = TRUE, error = TRUE,
                           fallback_to_prior = TRUE) {
  draws = mcmclist_draws(
    fit,
    prior = prior,
    message = message,
    error = error,
    fallback_to_prior = fallback_to_prior
  )
  if (is.null(draws))
    return(NULL)

  posterior::as_draws_array(draws)
}


#' Extract MCMC Draws from `mcpfit` Objects
#'
#' Extract posterior or prior draws using \pkg{posterior}, \pkg{tidybayes}, or \pkg{coda} S3 generics.
#'
#' @aliases as_draws as_draws.mcpfit as_draws_df.mcpfit as_draws_array.mcpfit as_draws_matrix.mcpfit as_draws_rvars.mcpfit as.mcmc.mcpfit tidy_draws.mcpfit
#' @param x An \code{\link{mcpfit}} object.
#' @param prior Logical. Extract prior draws (`TRUE`) instead of posterior draws
#'   (`FALSE`)? Errors if the requested draws are unavailable.
#' @param ... Passed to \pkg{posterior} or \pkg{tidybayes} format conversion functions.
#' @return A \pkg{posterior} `draws` object or a \pkg{coda} `mcmc.list` object.
#' @examples
#' # Default posterior draws, with one row per iteration and chain
#' draws = as_draws(demo_fit)  # Return a posterior::draws object
#' head(as_draws_df(demo_fit))  # Convert draws to a data frame
#'
#' # Other posterior formats are useful in different downstream packages
#' as_draws_matrix(demo_fit)[1:3, 1:3]  # Matrix of draws by parameter
#' as_draws_array(demo_fit)[1:2, , 1:2]  # Iteration-by-chain-by-parameter array
#' as_draws_rvars(demo_fit)[c("cp_1", "cp_2")]  # Random-variable representation
#'
#' # mcp also supports the coda and tidybayes conventions
#' head(coda::as.mcmc(demo_fit)[[1]])  # First chain as a coda mcmc object
#' head(tidybayes::tidy_draws(demo_fit))  # Tidybayes-compatible draw data
#' @exportS3Method posterior::as_draws
as_draws.mcpfit = function(x, prior = FALSE, ...) {
  posterior_draws(x, prior = prior, fallback_to_prior = FALSE)
}

#' @exportS3Method posterior::as_draws_df
as_draws_df.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_df(posterior_draws(x, prior = prior, fallback_to_prior = FALSE), ...)
}

#' @exportS3Method posterior::as_draws_array
as_draws_array.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_array(posterior_draws(x, prior = prior, fallback_to_prior = FALSE), ...)
}

#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_matrix(posterior_draws(x, prior = prior, fallback_to_prior = FALSE), ...)
}

#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_rvars(posterior_draws(x, prior = prior, fallback_to_prior = FALSE), ...)
}

#' @rdname as_draws.mcpfit
#' @export
as_draws = posterior::as_draws

#' @rdname as_draws.mcpfit
#' @export
as_draws_df = posterior::as_draws_df

#' @rdname as_draws.mcpfit
#' @export
as_draws_array = posterior::as_draws_array

#' @rdname as_draws.mcpfit
#' @export
as_draws_matrix = posterior::as_draws_matrix

#' @rdname as_draws.mcpfit
#' @export
as_draws_rvars = posterior::as_draws_rvars

#' @exportS3Method coda::as.mcmc
as.mcmc.mcpfit = function(x, prior = FALSE, ...) {
  mcmclist_draws(x, prior = prior)
}

#' @importFrom tidybayes tidy_draws
#' @exportS3Method tidybayes::tidy_draws
tidy_draws.mcpfit = function(model, ...) {
  posterior::as_draws_df(model, ...)
}


.onLoad = function(libname, pkgname) {
  if (requireNamespace("posterior", quietly = TRUE)) {
    registerS3method("as_draws", "mcpfit", as_draws.mcpfit, envir = asNamespace("posterior"))
    registerS3method("as_draws_df", "mcpfit", as_draws_df.mcpfit, envir = asNamespace("posterior"))
    registerS3method("as_draws_array", "mcpfit", as_draws_array.mcpfit, envir = asNamespace("posterior"))
    registerS3method("as_draws_matrix", "mcpfit", as_draws_matrix.mcpfit, envir = asNamespace("posterior"))
    registerS3method("as_draws_rvars", "mcpfit", as_draws_rvars.mcpfit, envir = asNamespace("posterior"))
    registerS3method("nchains", "mcpfit", nchains.mcpfit, envir = asNamespace("posterior"))
    registerS3method("ndraws", "mcpfit", ndraws.mcpfit, envir = asNamespace("posterior"))
    registerS3method("niterations", "mcpfit", niterations.mcpfit, envir = asNamespace("posterior"))
  }
  if (requireNamespace("coda", quietly = TRUE)) {
    registerS3method("as.mcmc", "mcpfit", as.mcmc.mcpfit, envir = asNamespace("coda"))
  }
  if (requireNamespace("tidybayes", quietly = TRUE)) {
    registerS3method("tidy_draws", "mcpfit", tidy_draws.mcpfit, envir = asNamespace("tidybayes"))
  }
  if (requireNamespace("rstantools", quietly = TRUE)) {
    registerS3method("posterior_epred", "mcpfit", posterior_epred.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("posterior_predict", "mcpfit", posterior_predict.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("posterior_linpred", "mcpfit", posterior_linpred.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("log_lik", "mcpfit", log_lik.mcpfit, envir = asNamespace("rstantools"))
  }
  setHook(packageEvent("rstantools", "onLoad"), function(...) {
    registerS3method("posterior_epred", "mcpfit", posterior_epred.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("posterior_predict", "mcpfit", posterior_predict.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("posterior_linpred", "mcpfit", posterior_linpred.mcpfit, envir = asNamespace("rstantools"))
    registerS3method("log_lik", "mcpfit", log_lik.mcpfit, envir = asNamespace("rstantools"))
  })
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


#' Index \code{mcpfit} objects
#'
#' Index variables, iterations, chains, and draws.
#'
#' @inheritParams fitted.mcpfit
#' @param x An `mcpfit` object or a posterior draws object.
#' @return An integer count of iterations, chains, or draws.
#' @name draws-index-mcp
#' @examples
#' niterations(demo_fit)
#' nchains(as_draws(demo_fit, prior = TRUE))
NULL


#' @aliases niterations.mcpfit
#' @describeIn draws-index-mcp Number of iterations per chain of an `mcpfit` object.
#' @exportS3Method posterior::niterations
niterations.mcpfit = function(x, ...) {
  coda::niter(mcmclist_draws(x))
}

#' @aliases nchains.mcpfit
#' @describeIn draws-index-mcp Number of chains of an `mcpfit` object.
#' @exportS3Method posterior::nchains
nchains.mcpfit = function(x, ...) {
  coda::nchain(mcmclist_draws(x))
}

#' @rdname draws-index-mcp
#' @exportS3Method posterior::ndraws
ndraws.mcpfit = function(x, ...) {
  draws = mcmclist_draws(x)
  sum(vapply(draws, nrow, integer(1)))
}

#' @rdname draws-index-mcp
#' @export
ndraws = posterior::ndraws

#' @rdname draws-index-mcp
#' @export
nchains = posterior::nchains

#' @rdname draws-index-mcp
#' @export
niterations = posterior::niterations

# Get information about group-level parameters
#
# Returns parameters, data columns, and effect metadata given parameter
# name(s), model part(s), or column(s).
#
# - pars: `NULL`/`FALSE` for nothing. `TRUE` for all. A character vector
#   containing `"cp"`, `"predictor"`, or exact group-level parameter names.
# - cols: `NULL`/`FALSE` for nothing. `TRUE` for all. A vector of grouping
#   column names for specifics. Usually provided via `facet_by` elsewhere.
# Returns: A list. See details.
#   
# Returns a list with
# * `pars`: Character vector of parameter names, or `NULL` if empty.
# * `cols`: Character vector of data column names, or `NULL` if empty.
# * `indices`: Logical vector indexing the group-effects table.
# * `effects`: The selected rows of the group-effects table.
unpack_group_effects = function(fit, pars = NULL, cols = NULL) {
  checkmate::assert_multi_class(pars, c("logical", "character"), null.ok = TRUE)
  checkmate::assert_multi_class(cols, c("logical", "character"), null.ok = TRUE)
  if (is.logical(pars))
    checkmate::assert_flag(pars)
  if (is.logical(cols))
    checkmate::assert_flag(cols)
  group_effects = get_fit_model_tables(fit)$group_effects
  use_group = rep(FALSE, nrow(group_effects))

  # If everything is NULL, just return NULLs
  if ((is.null(pars) && is.null(cols))) {
    return(list(
      pars = NULL,
      cols = NULL,
      indices = use_group,
      effects = group_effects[use_group, , drop = FALSE]
    ))
  } else if (!is.null(pars) && !is.null(cols)) {
    stop("One of `pars` and `cols` must be NULL.")
  }


  if (!is.null(pars)) {
    if (all(pars == FALSE)) {
      # Select no group-level effects
      use_group[] = FALSE
    } else if (all(pars == TRUE)) {
      # Select all group-level effects
      use_group[] = TRUE
    } else if (is.character(pars)) {
      allowed = c("cp", "predictor", group_effects$name)
      unknown = pars[pars %notin% allowed]
      if (length(unknown) > 0)
        stop(
          "Unknown group-effect selection: ", and_collapse(unknown), ". ",
          "Use TRUE, FALSE, \"cp\", \"predictor\", or a group-effect name."
        )
      use_group = group_effects$part %in% pars | group_effects$name %in% pars
    }
  } else if (!is.null(cols)) {
    if (all(cols == TRUE)) {
      use_group[] = TRUE
    } else if (!all(cols == FALSE)) {
      use_group = group_effects$group_col %in% cols
    }
  }

  # Return
  list(
    pars = logical0_to_null(group_effects$name[use_group]),
    cols = logical0_to_null(group_effects$group_col[use_group]),
    indices = use_group,
    effects = group_effects[use_group, , drop = FALSE]
  )
}



# Resolve the deprecated `nsamples` argument
resolve_ndraws = function(ndraws, nsamples, ndraws_missing, what,
                         samples = lifecycle::deprecated(),
                         env = rlang::caller_env(),
                         user_env = rlang::caller_env(2)) {
  if (lifecycle::is_present(samples)) {
    lifecycle::deprecate_soft(
      "0.4.0",
      paste0(what, "(samples)"),
      env = env,
      user_env = user_env
    )
  }
  if (lifecycle::is_present(nsamples)) {
    lifecycle::deprecate_soft(
      "0.4.0",
      paste0(what, "(nsamples)"),
      paste0(what, "(ndraws)"),
      env = env,
      user_env = user_env
    )
    if (!ndraws_missing)
      stop("Use only one of `ndraws` and deprecated `nsamples`.")
    ndraws = nsamples
  }
  ndraws
}


# Resolve the deprecated `samples_format` argument
resolve_draws_format = function(draws_format, samples_format, draws_format_missing, what,
                                env = rlang::caller_env(),
                                user_env = rlang::caller_env(2)) {
  if (lifecycle::is_present(samples_format)) {
    lifecycle::deprecate_soft(
      "0.4.0",
      paste0(what, "(samples_format)"),
      paste0(what, "(draws_format)"),
      env = env,
      user_env = user_env
    )
    if (!draws_format_missing)
      stop("Use only one of `draws_format` and deprecated `samples_format`.")
    draws_format = samples_format
  }
  rlang::arg_match0(draws_format, c("tidy", "matrix"))
}


#' Get tidy draws with or without group-level effects
#'
#' Extract posterior or prior draws formatted as tidy data frames
#'
#' Returns in a format useful for `fit$simulate()` with population-level parameters in wide format
#' and group-level deviations in long format (the number of rows is multiplied
#' by the number of selected group levels).
#'
#' @aliases mcp_draws
#' @keywords internal
#' @noRd
#' @inheritParams mcmclist_draws
#' @inheritParams pp_eval
#' @param population
#'   * `TRUE` All population-level model parameters.
#'   * `FALSE` No population-level effects. Same as `c()`.
#'   * Character vector: Only include specified population-level parameters.
#' @param varying Group-level effects. One of:
#'   * `TRUE` All group-level deviations.
#'   * `FALSE` No group-level deviations (`c()`).
#'   * `"cp"` or `"predictor"`: All group-level deviations belonging to that part of
#'     the model.
#'   * Character vector: Only include specified group-level parameters.
#' @param absolute
#'   * `TRUE` Returns the absolute location of all group-specific change points.
#'   * `FALSE` Return the group-level deviations.
#'   * Character vector: Apply the absolute transform only to these group-level parameters.
#'
#' @return `tibble` of posterior draws in `tidybayes` format.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
mcp_draws = function(
  fit,
  population = TRUE,
  varying = TRUE,
  absolute = FALSE,
  prior = FALSE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated()
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "mcp_draws")

  # General argument checks
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_multi_class(population, c("logical", "character"))
  checkmate::assert_multi_class(varying, c("logical", "character"), null.ok = TRUE)
  checkmate::assert_multi_class(absolute, c("logical", "character"), null.ok = TRUE)
  checkmate::assert_flag(prior)
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)

  if (all(population == FALSE) && all(varying == FALSE))
    stop("At least one TRUE or one parameter must be provided through either the `varying` or the `population` arguments.")


  # ----- IDENTIFY PARAMETERS -----
  # Group-level parameters formatted for tidybayes.
  group_info = unpack_group_effects(fit, pars = varying)
  group_terms = paste0(group_info$pars, "[", group_info$cols, "]")
  if (all(group_terms == "[]")) group_terms = ""  # quick fix

  # Population-level parameters. Result is `pars_population`.
  if (all(population == FALSE)) {
    pars_population = c()  # Empty if no absolute group-level change points
  } else if (all(population == TRUE)) {
    pars_population = mcp_pars(fit, scope = "population")$name
  } else if (is.character(population)) {
    if (!all(population %in% mcp_pars(fit, scope = "population")$name))
      stop("Not all `population` selections are population-level parameters.")

    pars_population = population
  }

  # Absolute effects. Results are `absolute_cps` and `absolute`.
  if (all(absolute == TRUE)) {
    cp_effects = group_info$effects[group_info$effects$part == "cp", , drop = FALSE]
    absolute = cp_effects$name
    absolute_cps = cp_effects$population_name
  } else if (all(absolute == FALSE)) {
    absolute_cps = NULL
  } else if (is.character(absolute)) {
    # Check
    is_group_selection = absolute %in% group_info$pars
    if (any(!is_group_selection))
      stop("The following parameter names in `absolute` are not in `varying`: ", and_collapse(absolute[!is_group_selection]))
    absolute_effects = group_info$effects[
      group_info$effects$name %in% absolute, , drop = FALSE
    ]
    if (any(absolute_effects$part != "cp"))
      stop("`absolute` can select change-point group-level effects only.")

    absolute_cps = absolute_effects$population_name
  }

  # ----- GET THESE PARAMETERS AS TIDY DRAWS -----
  # Select posterior/prior draws
  draws = mcmclist_draws(fit, prior = prior)

  # Build code for tidybayes::spread_draws() and execute it
  all_terms = unique(c(pars_population, group_terms, absolute_cps))
  code = paste0("tidybayes::spread_draws(draws, ", paste0(all_terms, collapse = ", "), ", ndraws = ndraws)")
  draws = eval(str2lang(code))

  # Preserve factor grouping columns from fit$data.
  if (length(group_info$cols) > 0) {
    is_factor = lapply(fit$data, is.factor)[group_info$cols]
    cols_to_factorize = group_info$cols[as.logical(is_factor)]
    draws = dplyr::mutate_at(draws, cols_to_factorize, as.factor)
  }

  # Add population-level change points to deviations, then remove helper columns.
  if (length(absolute_cps) > 0) {
    draws[, absolute_cps] = draws[, absolute_cps] + draws[, absolute]
    draws = dplyr::select(draws, -dplyr::all_of(absolute))
  }

  # Unassigned group-level deviations are simulated as zero (the population-level mean).
  remaining_group_cols = dplyr::setdiff(mcp_pars(fit, scope = "group")$name, colnames(draws))
  draws[, remaining_group_cols] = 0

  # Return with chain etc. first
  draws %>%
    dplyr::relocate(".chain", ".iteration", ".draw")
}


# Deprecated internal helper for MCMC draw extraction
tidy_samples = function(...) {
  lifecycle::deprecate_soft(
    when = "0.4.0",
    what = "tidy_samples()",
    with = "as_draws_df() or tidybayes::spread_draws()"
  )
  mcp_draws(...)
}


#' Fits and Predictions given Draws and data
#'
#' @aliases pp_eval pp_eval.mcpfit
#' @keywords internal
#' @param varying Group-level effects. One of:
#'   * `TRUE` All group-level deviations.
#'   * `FALSE` No group-level deviations (`c()`).
#'   * `"cp"` or `"predictor"`: All group-level deviations belonging to that part of
#'     the model.
#'   * Character vector: Only include specified group-level parameters.
#' @param object An `mcpfit` object.
#' @param newdata A `tibble` or a `data.frame` containing predictors in the model. 
#' - If `NULL` (default), the original data is used.
#' - For models with `ar()` or `ma()`: `fitted()`, `predict()`, `residuals()`, or `log_lik()` 
#'   conditions on the response history,
#'   so `newdata` must include the response. For `fitted()`, `predict()`, and
#'   `residuals()`, missing response histories are supported only in the original
#'   fitted data, using retained posterior imputations. `log_lik()` is unavailable
#'   when a missing response enters a later observed history. Use [`posterior_predict()`][rstantools::posterior_predict]
#'   to generate fresh response series recursively from predictor-only `newdata`.
#' - For models with `y | weights()`: Require the weights column except for `fitted()` and `predict()`.
#' @param summary Summarise at each x-value
#' @param type One of:
#'   - `"fitted"`: return the expected response. When `dpar` names a
#'     distributional parameter (e.g., `"mu"` or `"sigma"`), that parameter is returned instead.
#'     See also `fitted()`.
#'   - `"predict"`: return predicted values (e.g., `y_predict = rnorm(N, y_fitted, sigma_fitted)` for `family = gaussian()`).
#'     See also `predict()`.
#'   - `"residuals"`: observed y-values minus the fitted values. See also `residuals()`.
#'   - `"loglik"`: return the log-likelihood for each draw for each data point. See also `log_lik()`.
#'     Requires `scale = "response"`.
#' @param probs Vector of quantiles. Only in effect when `summary == TRUE`.
#' @param rate Logical scalar. For binomial models, return counts (`rate = FALSE`) or
#'   the observed or expected success proportion (`rate = TRUE`). Predictions and
#'   count-scale fitted values require a trials column in `newdata`.
#'   Distributional parameters such as `dpar = "mu"` evaluate the parameter itself (e.g., success probability)
#'   and are unaffected by `rate`.
#' @param prior Logical. Evaluate prior draws (`TRUE`) instead of posterior draws (`FALSE`, default)? Useful for `mcp(..., sample = "both")`.
#' @param dpar What distributional parameter to evaluate. This is only relevant when `type == "fitted"`. E.g.,
#'
#'   * `"epred"` (default): Expected response from the full model (or `NULL` for compatibility with brms etc.).
#'   * `"mu"`: The conditional mean (or success probability per trial for binomial/bernoulli models), on the link or response scale.
#'   * `"sigma"`: The standard deviation of the residuals.
#'   * `"ar1"`, `"ar2"`, `"ma1"`, `"ma2"`, etc. depending on which AR or MA
#'     coefficient you want to evaluate.
#' @param arma Whether to include AR and MA effects.
#'   * `TRUE` Compute the GARMA residual recurrence. Requires the response variable in `newdata`.
#'   * `FALSE` Disregard AR and MA effects. For `family = gaussian()`, `predict()` uses only `sigma` for residuals.
#'   For posterior evaluation of the original data, retained JAGS imputations
#'   supply missing GARMA histories. In models with group-level effects, this
#'   currently requires including all such effects (`varying = TRUE`).
#' @param ndraws Integer or `NULL`. Number of posterior draws to return/summarise.
#'   If there are group-level effects, this is the number of draws from each group.
#'   `NULL` means "all". More draws trade speed for accuracy.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @param draws_format One of "tidy" or "matrix". Controls the output format when `summary == FALSE` (for `fitted()`, `predict()`, and `log_lik()`). `residuals()` always returns tidy output.
#' @param samples_format Deprecated. Use `draws_format` instead.
#'   See more under "value"
#' @param scale One of
#'   * `"response"`: return on the observed scale, i.e., after applying the inverse link function.
#'   * `"linear"`: return on the linear-predictor (link) scale, where the linear
#'     trends are modeled.
#'     A linear scale is only applicable when `type == "fitted"` and `dpar` is not `NULL`.
#' @param .include_fitted Internal. Include fitted values with unsummarised predictions.
#' @param .include_dpars Internal. Include distributional parameters and response data as attributes with unsummarised predictions.
#' @param .garma_replicate Internal. For GARMA predictions, generate each
#'   response history recursively instead of conditioning on observed responses.
#' @return
#'   * If `summary = TRUE`: A data frame with the draw mean and SD (`sd`) for
#'     each row in `newdata`. With posterior draws (the default), `sd` is the
#'     posterior predictive SD for `type = "predict"` and the posterior SD of the
#'     evaluated quantity otherwise. With `prior = TRUE`, these are the analogous
#'     prior summaries. If `newdata` is `NULL`, the data in `fit$data` is used.
#'
#'   * If `summary = FALSE` and `draws_format = "tidy"`: A `tidybayes` `tibble` with all the posterior
#'     draws (`Nd`) evaluated at each row in `newdata` (`Nn`), i.e., with `Nd x Nn` rows. If there are
#'     group-level effects, the returned data is expanded with the relevant levels for each row.
#'
#'     The return columns are:
#'
#'      - Predictors from `newdata`, plus its response column when supplied.
#'      - Draw descriptors: ".chain", ".iteration", ".draw" (see the `posterior` and `tidybayes` packages), and `data_row`, the row number in the evaluated `newdata`.
#'      - Draw values: one column for each parameter in the model.
#'      - The estimate. Either ".epred", ".prediction", ".residual", or ".loglik" (matching tidybayes/ggdist conventions).
#'
#'   * If `summary = FALSE` and `draws_format = "matrix"`: An `N_draws` X `nrows(newdata)` matrix with fitted/predicted
#'       values (depending on `type`). This format is used by `brms` and it's useful as `yrep` in
#'      `bayesplot::ppc_*` functions.
#' @seealso \code{\link{fitted.mcpfit}} \code{\link{predict.mcpfit}} \code{\link{residuals.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
pp_eval = function(
  object,
  newdata = NULL,
  summary = TRUE,
  type = "fitted",
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  dpar = "epred",
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = 'response',
  .include_fitted = FALSE,
  .include_dpars = FALSE,
  .garma_replicate = FALSE,
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated()
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "pp_eval")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "pp_eval")

  # Recode
  fit = object
  checkmate::assert_class(fit, "mcpfit")
  warn_custom_jags_code(fit)
  if (!is.mcpfamily(fit$family))
    fit$family = mcpfamily(fit$family)
  if (is.null(fit$family$r$cdf))
    fit$family$r$cdf = mcpfamily(fit$family)$r$cdf
  dpar = assert_dpar(dpar, fit = fit, type = type)

  # What data to use
  using_original_data = is.null(newdata) || identical(data.frame(newdata), fit$data)
  if (using_original_data)
    newdata = fit$data

  data_columns = mcp_columns(fit)
  checkmate::assert_flag(.garma_replicate)
  replicate_garma = .garma_replicate && arma && is_arma(fit)
  assert_arma_series(newdata, data_columns$series)
  if (type == "loglik")
    assert_loglik_garma_history(fit, newdata, arma)

  conditional_garma = arma && is_arma(fit) && !replicate_garma &&
    (type %in% c("predict", "residuals") ||
       (type == "fitted" && dpar %in% c("epred", "mu")))
  if (conditional_garma && !using_original_data &&
      data_columns$response %in% names(newdata) && anyNA(newdata[[data_columns$response]]))
    stop(
      "Conditional GARMA evaluation with missing responses is supported only for ",
      "the original fitted data, where retained posterior imputations are available. ",
      "Use `posterior_predict()` to generate fresh replicated response series from ",
      "predictor-only `newdata`.",
      call. = FALSE
    )

  response_return = if (data_columns$response %in% colnames(newdata))
    newdata[, data_columns$response, drop = FALSE] else NULL


  ###############
  # FIX NEWDATA #
  ###############
  # Identify grouping columns to exclude based on the `varying` argument
  group_info = unpack_group_effects(fit, pars = varying)
  model_tables = get_fit_model_tables(fit)
  group_cols = unique(stats::na.omit(model_tables$group_effects$group_col))
  exclude_group_cols = setdiff(group_cols, c(group_info$cols, data_columns$series))

  # Determine which auxiliary columns are needed for this operation
  operation = switch(type, predict = "rng", loglik = "log_lik", fitted = "epred", residuals = "epred")
  aux_operations = c(operation, if (arma && is_arma(fit)) "garma")
  if (type == "fitted" && (rate || dpar != "epred"))
    aux_operations = setdiff(aux_operations, "epred")
  aux_columns = get_family_aux_columns(fit$family, model_tables$segments)
  aux_used = names(get_family_aux_columns(fit$family, model_tables$segments, aux_operations))
  unused_aux_columns = unname(aux_columns[names(aux_columns) %notin% aux_used])

  # Build list of required columns and validate presence in newdata
  required_cols = colnames(fit$data)  # Only predictive columns were saved in fit$data
  required_cols = required_cols[required_cols %notin% unused_aux_columns]
  required_cols = required_cols[required_cols %notin% exclude_group_cols]
  response_not_required = type %in% c("fitted", "predict") &&
    (arma == FALSE || is_arma(fit) == FALSE || replicate_garma)
  if (response_not_required) {
    required_cols = required_cols[required_cols != data_columns$response]
  } else if (data_columns$response %notin% colnames(newdata)) {
    stop("`newdata` must contain a response column named '", data_columns$response, "' for when `arma == TRUE` and/or `type == 'residuals'`")
  }
  assert_data_cols(newdata, required_cols)

  # Validate against reserved output namespace
  assert_reserved_output_namespace(colnames(newdata), context = "newdata")

  # Filter newdata columns and attach unique evaluation row index
  kept_cols = colnames(newdata)[colnames(newdata) %notin% exclude_group_cols]
  if (replicate_garma)
    kept_cols = kept_cols[kept_cols != data_columns$response]
  newdata = data.frame(newdata[, kept_cols, drop = FALSE])
  newdata$.mcp_data_row = seq_len(nrow(newdata))  # Evaluation key throughout summaries, matrices, plots, and metrics
  newdata_return = newdata
  if (!is.null(response_return) && data_columns$response %notin% colnames(newdata_return))
    newdata_return[[data_columns$response]] = response_return[[data_columns$response]]

  ########################
  # ASSERTS AND RECODING #
  ########################
  checkmate::assert_flag(summary)
  assert_typescale(type, scale)
  checkmate::assert(
    checkmate::check_flag(probs),
    checkmate::check_numeric(probs, any.missing = FALSE),
    .var.name = "probs"
  )
  if (is.numeric(probs))
    checkmate::assert_numeric(probs, lower = 0, upper = 1, any.missing = FALSE)
  if (is.logical(probs) && all(probs == TRUE))
    probs = c(0.025, 0.975)
  checkmate::assert_flag(rate)
  checkmate::assert_flag(prior)
  checkmate::assert_flag(arma)
  checkmate::assert_flag(.include_fitted)
  checkmate::assert_flag(.include_dpars)
  if (.include_fitted && (type != "predict" || summary))
    stop_github("`.include_fitted` requires `type = 'predict'` and `summary = FALSE`.")
  if (.include_dpars && (type != "predict" || summary))
    stop_github("`.include_dpars` requires `type = 'predict'` and `summary = FALSE`.")
  if (.garma_replicate && type != "predict")
    stop_github("`.garma_replicate` requires `type = 'predict'`.")
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)


  ########################
  # GET FITS/PREDICTIONS #
  ########################
  simulate_type = ifelse(type == "residuals", yes = "fitted", no = type)
  if (length(group_info$cols) > 0) {
    # Match group-level draws to each row of data.
    draws = dplyr::left_join(
      add_rhs_predictors(newdata, fit),
      mcp_draws(fit, population = TRUE, varying = varying, prior = prior, ndraws = ndraws),
      by = unique(group_info$cols),
      relationship = "many-to-many"
    )
  } else {
    # Without group-level effects, use all draws for each row of data.
    mcmc_draws = tibble::as_tibble(mcp_draws(fit, population = TRUE, varying = varying, prior = prior, ndraws = ndraws))
    predictors = tibble::as_tibble(add_rhs_predictors(newdata, fit))
    draws = dplyr::bind_cols(
      mcmc_draws[rep(seq_len(nrow(mcmc_draws)), each = nrow(predictors)), , drop = FALSE],
      predictors[rep(seq_len(nrow(predictors)), times = nrow(mcmc_draws)), , drop = FALSE]
    )
  }

  # Use imputed response draws for missing responses.
  # Requires special handling for varying and GARMA.
  imputed_response = rep(NA_real_, nrow(draws))
  has_posterior_draws = coda::is.mcmc.list(.subset2(fit, "mcmc_post"))
  needs_garma_history = arma && is_arma(fit) &&
    (type %in% c("predict", "residuals") ||
       (type == "fitted" && dpar %in% c("epred", "mu")))
  needs_imputed_response = !replicate_garma && (type == "predict" || needs_garma_history)
  if (needs_imputed_response && using_original_data && !prior && has_posterior_draws &&
      anyNA(fit$data[[data_columns$response]])) {
    full_varying = nrow(model_tables$group_effects) == 0 ||
      (is.logical(varying) && length(varying) == 1 && isTRUE(varying))
    if (arma && is_arma(fit) && !full_varying)
      stop(
        "This model has group-level effects, and its retained missing-response ",
        "histories are conditional on all of them. GARMA evaluation with missing ",
        "responses therefore currently requires `varying = TRUE`.",
        call. = FALSE
      )
    if (arma && is_arma(fit) && is.null(fit$.internal$imputed_response))
      stop(
        "This fit does not retain the missing response draws needed for coherent ",
        "GARMA evaluation. ",
        if (has_custom_jags_code(fit)) {
          "Automatic response imputation is unavailable with custom `jags_code`."
        } else {
          "Refit the model with the current version of mcp."
        },
        call. = FALSE
      )
    if (full_varying && (!is_arma(fit) || arma))
      imputed_response = get_imputed_response_draws(fit, draws)
    use_imputed = !is.na(imputed_response)
    if (arma && is_arma(fit) && any(use_imputed))
      draws[[data_columns$response]][use_imputed] = imputed_response[use_imputed]
  }

  # This is the important step!Evaluate the mcp model on newdata and draws.
  # Group-level joins are row-major, while GARMA recurrences require each
  # draw's data rows to be contiguous. Evaluate in draw/data order, then
  # restore the public row order below.
  evaluation_order = if (arma && is_arma(fit)) {
    order(draws$.draw, draws$.mcp_data_row)
  } else {
    seq_len(nrow(draws))
  }
  evaluation_data = draws[evaluation_order, , drop = FALSE]
  evaluate = function() rlang::exec(simulate_vectorized, fit, !!!evaluation_data, .type = simulate_type, .rate = rate, .dpar = dpar, .arma = arma, .scale = scale, .include_fitted = .include_fitted)
  evaluated = if (replicate_garma) suppressMessages(evaluate()) else evaluate()

  # Now more boilerplate stuff...
  fitted_values = attr(evaluated, "fitted")
  dpars_values = attr(evaluated, "dpars")
  response_data_values = attr(evaluated, "response_data")
  attr(evaluated, "fitted") = NULL
  attr(evaluated, "dpars") = NULL
  attr(evaluated, "response_data") = NULL
  restore_order = order(evaluation_order)
  evaluated = evaluated[restore_order]
  if (!is.null(fitted_values)) fitted_values = fitted_values[restore_order]
  if (!is.null(dpars_values)) dpars_values = lapply(dpars_values, function(v) v[restore_order])
  if (!is.null(response_data_values)) response_data_values = lapply(response_data_values, function(v) v[restore_order])
  if (type == "predict" && any(!is.na(imputed_response))) {
    response_data = get_family_response_data(fit$family, model_tables$segments, data = as.list(draws))
    imputed_return = fit$family$response$observed(imputed_response, response_data, rate)
    evaluated[!is.na(imputed_response)] = imputed_return[!is.na(imputed_response)]
  }
  draws[[type]] = evaluated

  # Plotting can request fitted and predicted values from the same evaluated
  # parameter rows and model evaluation.
  if (.include_fitted)
    draws$fitted = fitted_values

  if (!is.null(response_return))
    draws[[data_columns$response]] = response_return[[data_columns$response]][draws$.mcp_data_row]

  draws = draws %>% dplyr::select(-dplyr::starts_with(".pred_"))

  # Missing outcomes are latent in the fitted JAGS model, but they are not
  # observed-data likelihood contributions. Retain them while evaluating
  # GARMA histories above, then remove them from returned log likelihoods.
  if (type == "loglik") {
    observed_rows = which(!is.na(newdata[, data_columns$response]))
    if (length(observed_rows) == 0)
      stop("Log-likelihood evaluation requires at least one observed response.")
    draws = dplyr::filter(draws, .data$.mcp_data_row %in% observed_rows)
    newdata_return = dplyr::filter(
      newdata_return,
      .data$.mcp_data_row %in% observed_rows
    )
  }


  # Optionally compute residuals
  if (type == "residuals")
    draws = dplyr::mutate(draws, !!type := .data[[data_columns$response]] - .data[[type]])

  # Fail early if group-level joins or another evaluation step duplicated
  # or dropped any joint draw/evaluation-row combinations.
  validate_eval_draws(draws, type)

  # Optionally summarise
  if (summary == TRUE) {
    df_return = draws %>%
      # Summarise for each row in newdata
      dplyr::group_by(.data$.mcp_data_row) %>%
      dplyr::summarise(.groups = "drop",
                       sd = stats::sd(.data[[type]]),
                       !!type := mean(.data[[type]])
      ) %>%

      # Apply original order and put newdata as the first columns
      dplyr::arrange(.data$.mcp_data_row) %>%
      dplyr::left_join(newdata_return, by = ".mcp_data_row", relationship = "one-to-one") %>%
      dplyr::select(dplyr::one_of(colnames(newdata_return)), dplyr::all_of(type), "sd")


    # Quantiles
    if (!isFALSE(probs)) {
      val_col = if (type == "predict" && !is.null(fit$family$r$cdf)) ".predicted" else type
      quantiles = if (type == "predict" && !is.null(fit$family$r$cdf)) {
        get_mixture_quantiles(draws, probs, fit$family, keep = NULL, rate = rate, dpars = dpars_values, response_data = response_data_values)
      } else {
        get_quantiles(draws, probs, type, na.rm = type == "residuals")
      }
      quantiles = quantiles %>%
        dplyr::mutate(quantile = 100 * .data$quantile) %>%
        tidyr::pivot_wider(names_from = "quantile", names_prefix = "Q", values_from = dplyr::all_of(val_col))

      df_return = dplyr::left_join(df_return, quantiles, by = ".mcp_data_row", relationship = "one-to-one")
    }
    return(data.frame(dplyr::select(df_return, -".mcp_data_row")))
  } else if (draws_format == "tidy") {
    value_col = switch(type,
      fitted = ".epred",
      predict = ".prediction",
      residuals = ".residual",
      loglik = ".loglik",
      type
    )
    if (.include_fitted && "fitted" %in% colnames(draws)) {
      draws = dplyr::rename(draws, .epred = "fitted")
    }
    draws = dplyr::rename(draws, !!value_col := dplyr::all_of(type))
    if (.include_dpars) {
      if (!is.null(dpars_values)) attr(draws, "dpars") = dpars_values
      if (!is.null(response_data_values)) attr(draws, "response_data") = response_data_values
    }
    draws$data_row = draws$.mcp_data_row
    draws$.mcp_data_row = NULL
    return(draws)
  } else if (draws_format == "matrix") {
    df_return = tidy_to_matrix(draws, type)
    return(df_return)
  }
}



#' Fitted and predicted values of `mcp` models fits
#'
#' Evaluate the model on data, either summarised (per data-row) or per draw. You
#' can use draws from the prior (`prior = TRUE`), select a distributional
#' parameter with `dpar`, and choose the response or linear-predictor scale with
#' `scale` where applicable.
#'
#' @details
#' `residuals(fit)` is equivalent to `fit$data[[mcp_columns(fit)$response]] - fitted(fit, ...)` (or `newdata[[mcp_columns(fit)$response]] - fitted(fit, ...)`),
#' but with fixed arguments for `fitted`: `rate = FALSE, dpar = 'epred', draws_format = 'tidy'`.
#'
#' `log_lik()` defaults to an unsummarised draws-by-observation matrix, as used
#' by `loo` and other posterior workflows. Non-default `varying` and `arma`
#' settings evaluate conditional or counterfactual log-likelihoods (e.g.,
#' omitting random effects or serial correlation); they cannot be used in
#' `loo()` or `waic()` because estimating information criteria for reduced
#' models requires refitting.
#'
#' Missing responses in the original data remain missing in the response column.
#' `fitted()` returns their expected responses, while `predict()` uses retained
#' JAGS imputations for their posterior response draws. In GARMA models these
#' imputations also supply the history used for later fitted and predicted rows.
#'
#' @inheritParams pp_eval
#' @param ... Must be empty. Reserved for future use.
#' @inherit pp_eval return
#' @seealso \code{\link{fitted.mcpfit}} \code{\link{predict.mcpfit}} \code{\link{residuals.mcpfit}} \code{\link{log_lik.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' head(fitted(demo_fit))  # Expected response for each row of demo_fit$data
#' head(residuals(demo_fit))  # Residuals for each row of demo_fit$data
#' log_lik(demo_fit)[1:3, 1:3]  # Log-likelihood at each demo_fit$data
#'
#' # All of the above take a range of arguments. E.g.,:
#' \donttest{
#' head(predict(demo_fit))  # Pointwise posterior predictive
#' head(predict(demo_fit, probs = c(0.1, 0.5, 0.9)))  # Median and 80% posterior predictive interval.
#' head(predict(demo_fit, prior = TRUE))  # Prior predictive
#' head(fitted(demo_fit, summary = FALSE))  # Draws. Useful for plotting distributions.
#' head(fitted(demo_fit, dpar = "sigma"))  # Another model parameter
#'
#' # Evaluate at novel data
#' novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
#' head(predict(demo_fit, newdata = novel_data, probs = c(0.025, 0.5, 0.975)))
#'
#' # Work with missing responses
#' missing_fit = mcp_example("missing", plot = FALSE)
#' fitted(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Expected responses for missing y
#' fitted(missing_fit, summary = FALSE) |> dplyr::filter(is.na(y)) |> head()  # Same, but draws
#' predict(missing_fit) |> dplyr::filter(is.na(y)) |> head()  # Posterior predictive for missing y
#'}
#' @name execute-mcp-model
NULL


#' @aliases predict predict.mcpfit
#' @describeIn execute-mcp-model Predictive Distribution
#' @export
predict.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "predict.mcpfit")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "predict.mcpfit")
  dots = list(...)
  warn_which_y(dots, "predict")
  dots$which_y = NULL
  if (length(dots) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(names(dots)), call. = FALSE)

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "predict",
    probs = probs,
    rate = rate,
    prior = prior,
    dpar = NULL,
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    draws_format = draws_format
  )
}


#' @aliases fitted fitted.mcpfit
#' @describeIn execute-mcp-model Expected response
#' @export
fitted.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  dpar = "epred",
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "tidy",
  scale = "response",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "fitted.mcpfit")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "fitted.mcpfit")
  dots = list(...)
  warn_which_y(dots, "fitted")
  if ("which_y" %in% names(dots) && missing(dpar))
    dpar = dots$which_y
  dots$which_y = NULL
  if (length(dots) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(names(dots)), call. = FALSE)

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "fitted",
    probs = probs,
    rate = rate,
    prior = prior,
    dpar = dpar,
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    draws_format = draws_format,
    scale = scale
  )
}


#' Posterior prediction draws for `mcpfit` objects
#'
#' Methods for the `{rstantools}` posterior-prediction generics. They return a
#' draws-by-observation matrix and enable `{tidybayes}` workflows such as
#' `add_epred_draws()`, `add_predicted_draws()`, and `add_linpred_draws()`.
#' These methods and workflows require the suggested package `{rstantools}`.
#'
#' @param object An `mcpfit` object.
#' @param newdata Optional data frame at which to evaluate the model. For GARMA
#'   `posterior_predict()`, only predictors and required response auxiliaries are
#'   needed: each response series is generated recursively without conditioning
#'   on an observed response column.
#' @param draws,ndraws Number of posterior draws to return. `draws` follows the
#'   `{rstantools}` convention; `ndraws` is the mcp spelling. Supply at most one.
#' @param re.form,re_formula Group-level effects to include. `NULL` includes all
#'   effects and `NA` excludes them.
#' @param dpar Distributional parameter for `posterior_epred()` and
#'   `posterior_linpred()`; `NULL` uses the expected response.
#' @param transform For `posterior_linpred()`, return the inverse-link
#'   transformed expected response instead of the linear predictor.
#' @param seed Optional integer seed for draw selection and posterior prediction.
#' @param ... Must be empty. Reserved for future use.
#' @return A numeric `N_draws` by `nrow(newdata)` matrix.
#' @details For GARMA models, `posterior_predict()` generates each replicated
#'   response series recursively. It does not condition later predictions on
#'   the observed response history, unlike `fitted()` and `predict()`. These
#'   methods require posterior draws. For prior prediction, use `fitted()` or
#'   `predict()` with `prior = TRUE`.
#'
#'   For binomial models, `posterior_epred()` and `posterior_predict()` (and
#'   corresponding `{tidybayes}` workflows such as `add_epred_draws()`) follow
#'   `{brms}` and `{rstantools}` conventions by returning values on the outcome
#'   count scale (`rate = FALSE`), i.e., expected counts \eqn{E[Y] = n\mu} and
#'   simulated counts in \eqn{\{0, \dots, n\}}. In contrast, `fitted()` and `predict()`
#'   default to proportions (`rate = TRUE`). To obtain the success probability
#'   parameter \eqn{\mu} on the \eqn{[0, 1]} scale regardless of trial counts, pass
#'   `dpar = "mu"`.
#' @seealso [fitted.mcpfit()], [predict.mcpfit()]
posterior_epred.mcpfit = function(
  object,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  dpar = NULL,
  seed = NULL,
  ...
) {
  posterior_prediction_matrix(
    object = object,
    newdata = newdata,
    type = "fitted",
    draws = draws,
    ndraws = ndraws,
    re.form = re.form,
    re_formula = re_formula,
    dpar = dpar,
    scale = "response",
    seed = seed,
    ...
  )
}


#' @rdname posterior_epred.mcpfit
posterior_predict.mcpfit = function(
  object,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  seed = NULL,
  ...
) {
  posterior_prediction_matrix(
    object = object,
    newdata = newdata,
    type = "predict",
    draws = draws,
    ndraws = ndraws,
    re.form = re.form,
    re_formula = re_formula,
    seed = seed,
    ...
  )
}


#' @rdname posterior_epred.mcpfit
posterior_linpred.mcpfit = function(
  object,
  transform = FALSE,
  newdata = NULL,
  draws = NULL,
  ndraws = NULL,
  re.form = NULL,
  re_formula = NULL,
  dpar = NULL,
  seed = NULL,
  ...
) {
  checkmate::assert_flag(transform)
  posterior_prediction_matrix(
    object = object,
    newdata = newdata,
    type = "fitted",
    draws = draws,
    ndraws = ndraws,
    re.form = re.form,
    re_formula = re_formula,
    dpar = dpar,
    scale = if (transform) "response" else "linear",
    rate = TRUE,
    seed = seed,
    ...
  )
}


# Evaluate posterior prediction draws for rstantools-compatible methods
posterior_prediction_matrix = function(
  object,
  newdata,
  type,
  draws,
  ndraws,
  re.form,
  re_formula,
  dpar = NULL,
  scale = "response",
  rate = FALSE,
  seed = NULL,
  ...
) {
  checkmate::assert_class(object, "mcpfit")
  mcmclist_draws(object, message = FALSE, fallback_to_prior = FALSE)
  dots = list(...)
  if (length(dots) > 0)
    stop("Unrecognized argument(s): ", and_collapse(names(dots)), call. = FALSE)
  if (!is.null(seed))
    checkmate::assert_int(seed, lower = 1)
  if (!is.null(draws) && !is.null(ndraws))
    stop("Use only one of `draws` and `ndraws`.", call. = FALSE)
  if (!is.null(draws))
    ndraws = draws
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)

  varying = resolve_re_formula(re.form, re_formula)
  if (is.null(dpar))
    dpar = "epred"
  if (!is.null(seed))
    set.seed(seed)
  pp_eval(
    object,
    newdata = newdata,
    summary = FALSE,
    type = type,
    probs = FALSE,
    rate = rate,
    prior = FALSE,
    dpar = if (type == "fitted") dpar else NULL,
    varying = varying,
    arma = TRUE,
    ndraws = ndraws,
    draws_format = "matrix",
    scale = scale,
    .garma_replicate = type == "predict"
  )
}


# Convert rstantools group-effect syntax to mcp's group-effect selector
resolve_re_formula = function(re.form, re_formula) {
  if (!is.null(re_formula) && !is.null(re.form))
    stop("Use only one of `re.form` and `re_formula`.", call. = FALSE)
  formula = if (!is.null(re_formula)) re_formula else re.form
  if (is.null(formula))
    return(TRUE)
  if (length(formula) == 1 && is.na(formula))
    return(FALSE)
  stop("`re.form`/`re_formula` must be NULL or NA for mcpfit objects.", call. = FALSE)
}


#' @export
log_lik = function(object, ...) UseMethod("log_lik")

#' @aliases log_lik log_lik.mcpfit
#' @describeIn execute-mcp-model Pointwise log-likelihood
#' @export
log_lik.mcpfit = function(
  object,
  newdata = NULL,
  summary = FALSE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  draws_format = "matrix",
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "log_lik.mcpfit")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "log_lik.mcpfit")
  dots = list(...)
  warn_which_y(dots, "log_lik")
  dots$which_y = NULL
  if (length(dots) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(names(dots)), call. = FALSE)

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "loglik",
    probs = probs,
    rate = rate,
    prior = prior,
    dpar = NULL,
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    draws_format = draws_format,
    scale = "response"
  )
}


#' @aliases residuals residuals.mcpfit
#' @describeIn execute-mcp-model Residual distribution
#' @export
residuals.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "residuals.mcpfit")
  dots = list(...)
  warn_which_y(dots, "residuals")
  dots$which_y = NULL
  if (length(dots) > 0)
    stop("Unrecognized argument(s) passed in `...`: ", and_collapse(names(dots)), call. = FALSE)

  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "residuals",
    probs = probs,
    rate = FALSE,
    prior = prior,
    dpar = NULL,
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    draws_format = "tidy"
  )
}


# Add loo if not already present
#
# - fit: An mcpfit object
# - save_psis: Logical. See documentation of loo::loo
# - info: Optional message if adding loo
# - varying,arma: Evaluation settings passed to `loo.mcpfit()`.
# Returns: An mcpfit object with loo.
with_loo = function(fit, save_psis = FALSE, info = NULL,
                    varying = TRUE, arma = TRUE) {
  checkmate::assert_class(fit, "mcpfit")
  settings = get_loglik_settings(fit, varying, arma, ndraws = NULL)
  settings_match = identical(attr(fit$loo, "mcp_settings"), settings)
  needs_psis = save_psis == TRUE &&
    loo::is.loo(fit$loo) && is.null(fit$loo$psis_object)

  # Add loo if absent or needs psis
  if (is.null(fit$loo) || !settings_match || needs_psis) {
    if (is.character(info))
      message(info)
    fit$loo = loo(
      fit,
      save_psis = save_psis,
      varying = varying,
      arma = arma
    )
  }

  fit
}
