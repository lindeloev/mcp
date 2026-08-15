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
#' User-provided information (see \code{\link{mcp}} for more details):
#' @slot call The matched call to `mcp()`.
#' @slot model A list of formulas, making up the model.
#'   Provided by user. See \code{\link{mcp}} for more details.
#' @slot data A data frame.
#'   Provided by user. See \code{\link{mcp}} for more details.
#' @slot family An `mcpfamily` object.
#'   Provided by user. See \code{\link{mcp}} for more details.
#' @slot prior A named list.
#'   Provided by user. See \code{\link{mcp}} for more details.
#' @slot mcmc_post An \code{\link[coda]{mcmc.list}} object with posterior draws.
#' @slot mcmc_prior An \code{\link[coda]{mcmc.list}} object with prior draws.
#' @slot loglik An (Nchains * Ndraws) by N-observed-responses matrix of
#'   log-likelihoods.
#' @slot jags_code A string with jags code. Use `cat(fit$jags_code)` to show it.
#' @slot simulate A method to simulate and predict data.
#' @slot .internal Information that is used internally by mcp.
NULL


#' Internal function for summary.mcpfit, fixef.mcpfit, and ranef.mcpfit
#'
#' @aliases get_summary get_summary.mcpfit
#' @keywords internal
#' @noRd
#' @inheritParams summary.mcpfit
#' @param fit An \code{\link{mcpfit}}` object.
#' @param scope Which parameter scope to summarise: population-level parameters
#'   or group-level deviations.
#' @param role Optional parameter role to select within `scope`.
#' @param verbose Logical. Include the `segment` and `dpar` columns.
#' @return A data.frame with summaries for each model parameter. With
#'   `verbose = TRUE`, rows are labeled with `segment` and `dpar` columns (see
#'   `summary.mcpfit`).
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_summary = function(fit, width, scope = c("population", "group"), role = NULL,
                       prior = FALSE, verbose = FALSE) {
  # Check arguments
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_number(width, lower = 0, upper = 1)
  scope = rlang::arg_match0(scope, c("population", "group"))
  checkmate::assert_character(role, null.ok = TRUE)
  checkmate::assert_flag(prior)
  checkmate::assert_flag(verbose)

  if (scope == "group" & nrow(mcp_pars(fit, scope = "group")) == 0)
    return(NULL)

  draws = posterior_draws(fit, prior = prior)

  # Select by the independent scope and role dimensions of the parameter table.
  all_cols = posterior::variables(draws)
  pars = mcp_pars(fit)
  selected_names = pars$name[pars$scope == scope]
  if (!is.null(role))
    selected_names = selected_names[pars$role[pars$scope == scope] %in% role]

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
      name = character(), mean = numeric(), sd = numeric(), lower = numeric(),
      upper = numeric(), rhat = numeric(), ess_bulk = numeric(), ess_tail = numeric()
    )
    if (verbose) {
      estimates$segment = integer()
      estimates$dpar = character()
    }
    if (!is.null(attr(fit$data[, mcp_columns(fit)$response], "simulated"))) {
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
    dplyr::rename(name = "variable") %>%
    dplyr::mutate(
      ess_bulk = round(.data$ess_bulk),
      ess_tail = round(.data$ess_tail)
  )

  # Order rows and add `segment`/`dpar` using the canonical parameter table
  # built in mcp(). Group-level columns (e.g. "cp_1_id[A]") are matched by
  # their base name; ties (i.e., levels of the same group-level effect) are
  # broken alphabetically by the full column name.
  base_name = sub("\\[.*\\]$", "", estimates$name)
  match_idx = match(base_name, pars$name)
  estimates$segment = pars$segment[match_idx]
  estimates$dpar = pars$dpar[match_idx]
  estimates = estimates[order(match_idx, estimates$name), ]

  # Add simulation parameters if the data is simulated
  sim_list = attr(fit$data[, mcp_columns(fit)$response], "simulated")
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
      name = names(simulated),
      sim = as.numeric(simulated),  # without row names
      stringsAsFactors = FALSE
    )

    # Add simulation values for comparison with the fitted parameters.
    estimates = estimates %>%
      dplyr::left_join(simulated, by = "name", relationship = "one-to-one") %>%
      dplyr::mutate(
        cp_width = ifelse(stringr::str_detect(.data$name, "^cp_[0-9]+"), .data$upper - .data$lower, 0),
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
      "name", "mean", "sd", "lower", "upper", "rhat", "ess_bulk", "ess_tail",
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
#' @param prior TRUE/FALSE. Summarise prior instead of posterior?
#' @param verbose Logical. Include the `segment` and `dpar` columns. Defaults
#'   to `FALSE` for a compact, v0.3.4-compatible summary.
#' @inheritParams mcp
#' @param ... Currently ignored
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
#'     `"ar1"`, `"ma1"`, etc.) the parameter belongs to.
#'   * `mean` is the posterior mean
#'   * `sd` is the posterior standard deviation, i.e., the width of the posterior.
#'   * `lower` and `upper` are the bounds of the central posterior interval
#'     given in `width`.
#'   * `rhat` is the rank-normalized split-Rhat convergence diagnostic.
#'   * `ess_bulk` and `ess_tail` are the bulk and tail effective sample sizes.
#'     Low effective sample sizes are also obvious as poor mixing in trace plots
#'     (see `plot_pars(fit)`). Read how to deal with such problems [here](https://lindeloev.github.io/mcp/articles/tips.html)
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
  cat("Family: ", fit$family$family, "(link = '", fit$family$link, "')\n", sep = "")
  if (!is.null(draws))
    cat("Iterations: ", coda::niter(draws), " from ", coda::nchain(draws), " chains.\n", sep="")
  cat("Segments:\n")
  for (i in 1:length(fit$model)) {
    cat("  ", i, ": ", formula_to_char(fit$model[[i]]), "\n", sep = "")
  }

  # Data
  if (!is.null(draws)) {
    # Print and return population-level summaries invisibly.
    result = get_summary(fit, width, scope = "population", prior = prior, verbose = verbose)
    pars = mcp_pars(fit)
    cp_names = pars$name[pars$part == "cp" & pars$scope == "population"]
    is_cp = result$name %in% cp_names

    # Format before splitting, so both printed tables share column widths.
    display = data.frame(lapply(result, format, digits = digits), check.names = FALSE)
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



#' @export
fixef = function(object, ...) UseMethod("fixef")

#' @export
ranef = function(object, ...) UseMethod("ranef")


#' @aliases fixef fixef.mcpfit
#' @describeIn summary.mcpfit Population-level fixed effects (regression coefficients) of `mcpfit`.
#' @export
fixef.mcpfit = function(object, width = 0.95, prior = FALSE, verbose = FALSE, ...) {
  rlang::check_dots_empty()
  get_summary(object, width, scope = "population", role = "fixed_effect", prior = prior, verbose = verbose)
}

#' @aliases ranef ranef.mcpfit
#' @describeIn summary.mcpfit Group-level deviations (random effects) of `mcpfit`.
#' @export
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
#' @param ... Currently unused.
#' @return `formula()` returns the complete list of segment formulas, or one
#'   formula when `segment` is supplied. `family()` returns an `mcpfamily`.
#'   `model.frame()` returns the data retained in the fit. `nobs()` returns the
#'   number of fitting-data rows.
#' @name model-accessors-mcpfit
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
  nrow(model.frame(object))
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


#' Posterior Covariance and Central Intervals for `mcpfit` Objects
#'
#' Summarise the joint and marginal posterior uncertainty of population-level
#' model parameters.
#'
#' @param object An `mcpfit` object.
#' @param correlation Return the posterior correlation matrix instead of the
#'   covariance matrix?
#' @param pars Optional names of population-level parameters to extract.
#' @param parm Optional names or positions of population-level parameters to
#'   include in the intervals.
#' @param level Width of the central posterior interval.
#' @param ... Currently unused.
#' @return `vcov()` returns a posterior covariance or correlation matrix.
#'   `confint()` returns a two-column matrix of central posterior intervals.
#' @name posterior-uncertainty-mcpfit
NULL


#' @rdname posterior-uncertainty-mcpfit
#' @export
vcov.mcpfit = function(object, correlation = FALSE, pars = NULL, ...) {
  rlang::check_dots_empty()
  checkmate::assert_flag(correlation)
  population = mcp_pars(object, scope = "population")$name
  pars = if (is.null(pars)) population else as.character(pars)
  pars = intersect(population, pars)
  if (length(pars) == 0)
    return(NULL)

  draws = posterior::as_draws_matrix(posterior_draws(object))
  if (correlation)
    return(stats::cor(draws[, pars, drop = FALSE]))

  stats::cov(draws[, pars, drop = FALSE])
}


#' @rdname posterior-uncertainty-mcpfit
#' @export
confint.mcpfit = function(object, parm, level = 0.95, ...) {
  rlang::check_dots_empty()
  checkmate::assert_number(level, lower = 0, upper = 1)
  checkmate::assert_true(level > 0 && level < 1, .var.name = "level")

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
  draws = posterior::as_draws_matrix(posterior_draws(object))
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
#' @export
is.mcpfit = function(x) {
  inherits(x, "mcpfit")
}


#' Internal function to get draws.
#'
#' Returns posterior draws, if available. If not, then prior draws. If not,
#' then throw an informative error. This is useful for summary and plotting, that
#' works on both.
#'
#' @aliases mcmclist_draws mcmclist_draws.mcpfit
#' @keywords internal
#' @inheritParams summary.mcpfit
#' @param fit An \code{\link{mcpfit}} object
#' @param message TRUE: gives a message if returning prior draws. FALSE = no message
#' @param error TRUE: err if there are no draws. FALSE: return NULL
mcmclist_draws = function(fit, prior = FALSE, message = TRUE, error = TRUE) {
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
  } else if (coda::is.mcmc.list(mcmc_prior)) {
    if (message)
      message("Posterior was not drawn. Using prior draws. Set `prior = TRUE` to mute this message.")
    return(mcmc_prior)
  } else if (error == TRUE) {
    stop("This mcpfit contains no posterior or prior draws.")
  }

  NULL
}


#' Get draws as a posterior draws array
#'
#' This is the single internal conversion from the stored
#' \code{\link[coda]{mcmc.list}} representation to a posterior draws object.
#'
#' @keywords internal
#' @noRd
#' @inheritParams mcmclist_draws
posterior_draws = function(fit, prior = FALSE, message = TRUE, error = TRUE) {
  draws = mcmclist_draws(
    fit,
    prior = prior,
    message = message,
    error = error
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
#' @param prior Logical. Extract prior draws (`TRUE`) instead of posterior draws (`FALSE`)?
#' @param ... Passed to \pkg{posterior} or \pkg{tidybayes} format conversion functions.
#' @return A \pkg{posterior} `draws` object or a \pkg{coda} `mcmc.list` object.
#' @exportS3Method posterior::as_draws
as_draws.mcpfit = function(x, prior = FALSE, ...) {
  posterior_draws(x, prior = prior)
}

#' @exportS3Method posterior::as_draws_df
as_draws_df.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_df(posterior_draws(x, prior = prior), ...)
}

#' @exportS3Method posterior::as_draws_array
as_draws_array.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_array(posterior_draws(x, prior = prior), ...)
}

#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_matrix(posterior_draws(x, prior = prior), ...)
}

#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.mcpfit = function(x, prior = FALSE, ...) {
  posterior::as_draws_rvars(posterior_draws(x, prior = prior), ...)
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

#' @exportS3Method tidybayes::tidy_draws
tidy_draws.mcpfit = function(x, ...) {
  posterior::as_draws_df(x, ...)
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
  }
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

#' Get information about group-level parameters
#'
#' Returns parameters, data columns, and effect metadata given parameter
#' name(s), model part(s), or column(s).
#'
#' @aliases unpack_group_effects unpack_group_effects.mcpfit
#' @keywords internal
#' @noRd
#' @param pars `NULL`/`FALSE` for nothing. `TRUE` for all. A character vector
#'   containing `"cp"`, `"predictor"`, or exact group-level parameter names.
#' @param cols `NULL`/`FALSE` for nothing. `TRUE` for all. A vector of grouping
#'   column names for specifics. Usually provided via `facet_by` elsewhere.
#' @return A list. See details.
#'
#' @details
#' Returns a list with
#' @slot pars Character vector of parameter names. `NULL` if empty.
#' @slot cols Character vector of data column names. `NULL` if empty.
#' @slot indices Logical vector indexing the group-effects table.
#' @slot effects The selected rows of the group-effects table.
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



#' Resolve the deprecated `nsamples` argument
#'
#' @keywords internal
#' @noRd
resolve_ndraws = function(ndraws, nsamples, ndraws_missing, what,
                         env = rlang::caller_env(),
                         user_env = rlang::caller_env(2)) {
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


#' Resolve the deprecated `samples_format` argument
#'
#' @keywords internal
#' @noRd
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
#' @inheritParams mcmclist_draws
#' @inheritParams pp_eval
#' @param population
#'   * `TRUE` All population-level model parameters.
#'   * `FALSE` No population-level effects. Same as `c()`.
#'   * Character vector: Only include specified population-level parameters.
#' @param varying One of:
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
    # draws[, absolute] = draws[, absolute] + draws[, absolute_cps]
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


#' Deprecated internal helper for MCMC draw extraction
#' @keywords internal
#' @noRd
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
#' @inheritParams mcp_draws
#' @param object An `mcpfit` object.
#' @param newdata A `tibble` or a `data.frame` containing predictors in the model. Weighted
#'   Gaussian predictions and log-likelihoods also require the weights column. If `NULL`
#'   (default), the original data is used.
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
#' @param rate Boolean. For binomial models, return counts (`rate = FALSE`) or
#'   the observed or expected success proportion (`rate = TRUE`). If `FALSE`, linear
#'   interpolation on trial number is used to infer trials at a particular x.
#' @param prior TRUE/FALSE. Plot using prior draws? Useful for `mcp(..., sample = "both")`
#' @param dpar What distributional parameter to evaluate. This is only relevant when `type == "fitted"`. E.g.,
#'
#'   * `"epred"` (default): Expected response from the full model (or `NULL` for compatibility with brms etc.).
#'   * `"mu"`: The central tendency which is often the mean after applying the
#'     link function.
#'   * `"sigma"`: The standard deviation of the residuals.
#'   * `"ar1"`, `"ar2"`, `"ma1"`, `"ma2"`, etc. depending on which AR or MA
#'     coefficient you want to evaluate.
#' @param arma Whether to include AR and MA effects.
#'   * `TRUE` Compute the GARMA residual recurrence. Requires the response variable in `newdata`.
#'   * `FALSE` Disregard AR and MA effects. For `family = gaussian()`, `predict()` uses only `sigma` for residuals.
#' @param ndraws Integer or `NULL`. Number of posterior draws to return/summarise.
#'   If there are group-level effects, this is the number of draws from each group.
#'   `NULL` means "all". Ignored if both are `FALSE`. More draws trade speed for accuracy.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @param draws_format One of "tidy" or "matrix". Controls the output format when `summary == FALSE`.
#' @param samples_format Deprecated. Use `draws_format` instead.
#'   See more under "value"
#' @param scale One of
#'   * `"response"`: return on the observed scale, i.e., after applying the inverse link function.
#'   * `"linear"`: return on the linear-predictor (link) scale, where the linear
#'     trends are modeled.
#'     A linear scale is only applicable when `type == "fitted"` and `dpar` is not `NULL`.
#' @param .include_fitted Internal. Include fitted values with unsummarised predictions.
#' @return
#'   * If `summary = TRUE`: A data frame with the draw mean and SD (`error`) for
#'     each row in `newdata`. With posterior draws (the default), `error` is the
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
#'      - Predictors from `newdata`.
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
  nsamples = lifecycle::deprecated(),
  samples_format = lifecycle::deprecated()
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "pp_eval")
  draws_format = resolve_draws_format(draws_format, samples_format, missing(draws_format), "pp_eval")

  # Recode
  fit = object
  checkmate::assert_class(fit, "mcpfit")
  if (!is.mcpfamily(fit$family))
    fit$family = mcpfamily(fit$family)
  dpar = assert_dpar(dpar, fit = fit, type = type)

  if (is.null(newdata))
    newdata = fit$data

  data_columns = mcp_columns(fit)
  assert_arma_series(newdata, data_columns$series)


  ###############
  # FIX NEWDATA #
  ###############
  group_info = unpack_group_effects(fit, pars = varying)
  model_tables = get_fit_model_tables(fit)
  group_cols = unique(stats::na.omit(model_tables$group_effects$group_col))
  exclude_group_cols = setdiff(group_cols, c(group_info$cols, data_columns$series))
  required_cols = colnames(fit$data)  # Only predictive columns were saved in fit$data
  operation = switch(type, predict = "rng", loglik = "log_lik", fitted = "epred", residuals = "epred")
  aux_operations = c(operation, if (arma && is_arma(fit)) "garma")
  aux_columns = get_family_aux_columns(fit$family, model_tables$segments)
  aux_used = names(get_family_aux_columns(fit$family, model_tables$segments, aux_operations))
  unused_aux_columns = unname(aux_columns[names(aux_columns) %notin% aux_used])
  required_cols = required_cols[required_cols %notin% unused_aux_columns]
  required_cols = required_cols[required_cols %notin% exclude_group_cols]
  if ((arma == FALSE || is_arma(fit) == FALSE) & type %in% c("fitted", "predict")) {
    required_cols = required_cols[required_cols != data_columns$response]
  } else if (data_columns$response %notin% colnames(newdata)) {
    stop("`newdata` must contain a response column named '", data_columns$response, "' for when `arma == TRUE` and/or `type == 'residuals'`")
  }
  assert_data_cols(newdata, required_cols)  # Helpful error if something is missing
  newdata = data.frame(newdata[, required_cols, drop = FALSE])
  #colnames(newdata) = required_cols  # Special case for when there's only one predictor
  newdata$data_row = seq_len(nrow(newdata))  # Evaluation key throughout summaries, matrices, plots, and metrics
  point_size_col = aux_columns[fit$family$response$point_size]
  point_size_col = unname(point_size_col[!is.na(point_size_col)])
  newdata_return = dplyr::select(newdata, -dplyr::any_of(point_size_col))

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
  if (is.logical(probs) & all(probs == TRUE))
    probs = c(0.025, 0.975)
  checkmate::assert_flag(rate)
  checkmate::assert_flag(prior)
  checkmate::assert_flag(arma)
  checkmate::assert_flag(.include_fitted)
  if (.include_fitted && (type != "predict" || summary))
    stop_github("`.include_fitted` requires `type = 'predict'` and `summary = FALSE`.")
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)


  ########################
  # GET FITS/PREDICTIONS #
  ########################
  simulate_type = ifelse(type == "residuals", yes = "fitted", no = type)
  if (length(group_info$cols) > 0) {
    # Match group-level draws to each row of data.
    draws_predictors = dplyr::left_join(
      add_rhs_predictors(newdata, fit),
      mcp_draws(fit, population = TRUE, varying = varying, prior = prior, ndraws = ndraws),
      by = unique(group_info$cols),
      relationship = "many-to-many"
    )
  } else {
    # Without group-level effects, use all draws for each row of data.
    draws = tibble::as_tibble(mcp_draws(fit, population = TRUE, varying = varying, prior = prior, ndraws = ndraws))
    predictors = tibble::as_tibble(add_rhs_predictors(newdata, fit))
    draws_predictors = dplyr::bind_cols(
      draws[rep(seq_len(nrow(draws)), each = nrow(predictors)), , drop = FALSE],
      predictors[rep(seq_len(nrow(predictors)), times = nrow(draws)), , drop = FALSE]
    )
  }

  draws = draws_predictors
  evaluated = rlang::exec(simulate_vectorized, fit, !!!draws_predictors, .type = simulate_type, .rate = rate, .dpar = dpar, .arma = arma, .scale = scale, .include_fitted = .include_fitted)
  fitted_values = attr(evaluated, "fitted")
  attr(evaluated, "fitted") = NULL
  draws[[type]] = evaluated

  # Plotting can request fitted and predicted values from the same evaluated
  # parameter rows and model evaluation.
  if (.include_fitted)
    draws$fitted = fitted_values

  draws = draws %>% dplyr::select(-dplyr::starts_with(".pred_"), -dplyr::any_of(point_size_col))

  # Missing outcomes are latent in the fitted JAGS model, but they are not
  # observed-data likelihood contributions. Retain them while evaluating
  # GARMA histories above, then remove them from returned log likelihoods.
  if (type == "loglik") {
    observed_rows = which(!is.na(newdata[, data_columns$response]))
    if (length(observed_rows) == 0)
      stop("Log-likelihood evaluation requires at least one observed response.")
    draws = dplyr::filter(draws, .data$data_row %in% observed_rows)
    newdata_return = dplyr::filter(
      newdata_return,
      .data$data_row %in% observed_rows
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
      dplyr::group_by(.data$data_row) %>%
      dplyr::summarise(.groups = "drop",
                       error = stats::sd(.data[[type]]),
                       !!type := mean(.data[[type]])
      ) %>%

      # Apply original order and put newdata as the first columns
      dplyr::arrange(.data$data_row) %>%
      dplyr::left_join(newdata_return, by = "data_row", relationship = "one-to-one") %>%
      dplyr::select(dplyr::one_of(colnames(newdata_return)), dplyr::all_of(type), "error")


    # Quantiles
    if (all(probs != FALSE)) {
      quantiles = get_quantiles(draws, probs, type) %>%
        dplyr::mutate(quantile = 100 * .data$quantile) %>%
        tidyr::pivot_wider(names_from = "quantile", names_prefix = "Q", values_from = dplyr::all_of(type))

      df_return = dplyr::left_join(df_return, quantiles, by = "data_row", relationship = "one-to-one")
    }
    return(data.frame(dplyr::select(df_return, -"data_row")))
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
#' `residuals(fit)` is equivalent to `fit$data[, fit$data$yvar] - fitted(fit, ...)` (or `newdata[, fit$data$yvar] - fitted(fit, ...)`),
#' but with fixed arguments for `fitted`: `rate = FALSE, dpar = 'epred', draws_format = 'tidy'`.
#'
#' `log_lik()` defaults to an unsummarised draws-by-observation matrix, as used
#' by `loo` and other posterior workflows.
#'
#' @inheritParams pp_eval
#' @param ... Currently ignored.
#' @inherit pp_eval return
#' @seealso \code{\link{fitted.mcpfit}} \code{\link{predict.mcpfit}} \code{\link{residuals.mcpfit}} \code{\link{log_lik.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' fitted(demo_fit)  # Expected response for each row of demo_fit$data
#' residuals(demo_fit)  # Residuals for each row of demo_fit$data
#' log_lik(demo_fit)  # Log-likelihood at each demo_fit$data
#'
#' # All of the above take a range of arguments. E.g.,:
#' \donttest{
#' predict(demo_fit)  # Pointwise posterior predictive
#' predict(demo_fit, probs = c(0.1, 0.5, 0.9))  # Median and 80% posterior predictive interval.
#' predict(demo_fit, prior = TRUE)  # Prior predictive
#' fitted(demo_fit, summary = FALSE)  # Draws instead of summary. Useful for plotting distributions.
#' fitted(demo_fit, dpar = "sigma")  # Another model parameter
#'
#' # Evaluate at novel data
#' novel_data = data.frame(time = c(-5, 20, 300))  # Only predictors are needed
#' predict(demo_fit, newdata = novel_data, probs = c(0.025, 0.5, 0.975))
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
#'
#' @param object An `mcpfit` object.
#' @param newdata Optional data frame at which to evaluate the model.
#' @param draws,ndraws Number of posterior draws to return. `draws` follows the
#'   `{rstantools}` convention; `ndraws` is the mcp spelling. Supply at most one.
#' @param re.form,re_formula Group-level effects to include. `NULL` includes all
#'   effects and `NA` excludes them.
#' @param dpar Distributional parameter for `posterior_epred()` and
#'   `posterior_linpred()`; `NULL` uses the expected response.
#' @param transform For `posterior_linpred()`, return the inverse-link
#'   transformed expected response instead of the linear predictor.
#' @param seed Accepted for `{tidybayes}` compatibility. Randomness is controlled
#'   by the calling context.
#' @param ... Currently ignored.
#' @return A numeric `N_draws` by `nrow(newdata)` matrix.
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
    seed = seed,
    ...
  )
}


#' Evaluate posterior prediction draws for rstantools-compatible methods
#' @keywords internal
#' @noRd
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
  seed = NULL,
  ...
) {
  checkmate::assert_class(object, "mcpfit")
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
  pp_eval(
    object,
    newdata = newdata,
    summary = FALSE,
    type = type,
    probs = FALSE,
    rate = TRUE,
    prior = FALSE,
    dpar = if (type == "fitted") dpar else NULL,
    varying = varying,
    arma = TRUE,
    ndraws = ndraws,
    draws_format = "matrix",
    scale = scale
  )
}


#' Convert rstantools group-effect syntax to mcp's group-effect selector
#' @keywords internal
#' @noRd
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


#' Add loo if not already present
#'
#' @aliases with_loo
#' @keywords internal
#' @noRd
#' @param fit An mcpfit object
#' @param save_psis Logical. See documentation of loo::loo
#' @param info Optional message if adding loo
#' @param varying,arma Evaluation settings passed to `loo.mcpfit()`.
#' @return An mcpfit object with loo.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
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
