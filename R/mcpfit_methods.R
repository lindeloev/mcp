# ABOUT: These are non-plotting functions that take an mcpfit as the first argument
# -----------------

#' Class `mcpfit` of Models Fitted with the \pkg{mcp} Package
#'
#' Models fitted with the \code{\link[mcp:mcp]{mcp}} function are represented as
#' an `mcpfit` object which contains the user input (model, data, family),
#' derived model characteristics (prior, parameter names, and jags code), and
#' the fit (prior and/or posterior mcmc samples).
#'
#' @name mcpfit-class
#' @aliases mcpfit
#' @docType class
#'
#' @details
#' See `methods(class = "mcpfit")` for an overview of available methods.
#'
#' User-provided information (see \code{\link{mcp}} for more details):
#' @slot model A list of formulas, making up the model.
#'   Provided by user. See \code{\link{mcp}} for more details.
#' @slot data A data frame.
#'   Provided by user. See \code{\link{mcp}} for more details.
#' @slot family An `mcpfamily` object.
#'   Provided by user. See \code{\link{mcp}} for more details.
#' @slot prior A named list.
#'   Provided by user. See \code{\link{mcp}} for more details.
#' @slot mcmc_post An \code{\link[coda]{mcmc.list}} object with posterior samples.
#' @slot mcmc_prior An \code{\link[coda]{mcmc.list}} object with prior samples.
#' @slot loglik An (Nchains * Ndraws) by N-observed-responses matrix of
#'   log-likelihoods.
#' @slot pars A list of character vectors of model parameter names.
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
#' @param varying Boolean. Get results for varying (TRUE) or population (FALSE)?
#' @return A data.frame with summaries for each model parameter, ordered and
#'   labeled with `segment` and `dpar` columns (see `summary.mcpfit`).
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_summary = function(fit, width, varying = FALSE, prior = FALSE) {
  # Check arguments
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_number(width, lower = 0, upper = 1)
  checkmate::assert_flag(varying)
  checkmate::assert_flag(prior)

  if (varying == TRUE & is.null(fit$pars$varying))
    return(NULL)

  samples = posterior_draws(fit, prior = prior)

  # Select only group-indexed or only population-scope columns.
  all_cols = posterior::variables(samples)
  if (varying == FALSE) {
    get_cols = all_cols[all_cols %in% fit$pars$population]
  } else {
    get_cols = all_cols[vapply(
      all_cols,
      function(column) any(startsWith(column, paste0(fit$pars$varying, "["))),
      logical(1)
    )]
    if (length(get_cols) == 0)
      stop("There were no matching parameters in the model.")
  }

  samples = posterior::subset_draws(samples, variable = get_cols)

  # Get parameter estimates and diagnostics
  tail_prob = (1 - width) / 2
  estimates = posterior::summarise_draws(
    samples,
    mean = base::mean,
    lower = function(x) stats::quantile(x, tail_prob, names = FALSE),
    upper = function(x) stats::quantile(x, 1 - tail_prob, names = FALSE),
    Rhat = posterior::rhat,
    ess_bulk = function(x) suppressWarnings(posterior::ess_bulk(x)),
    ess_tail = function(x) suppressWarnings(posterior::ess_tail(x))
  ) %>%
    dplyr::rename(name = "variable") %>%
    dplyr::mutate(
      ess_bulk = round(.data$ess_bulk),
      ess_tail = round(.data$ess_tail)
  )

  # Order rows and add `segment`/`dpar` using the canonical parameter table
  # built in mcp(). Varying-effect columns (e.g. "cp_1_id[A]") are matched by
  # their base name; ties (i.e., levels of the same varying effect) are
  # broken alphabetically by the full column name.
  pars = get_fit_model_tables(fit)$pars
  base_name = sub("\\[.*\\]$", "", estimates$name)
  match_idx = match(base_name, pars$name)
  estimates$segment = pars$segment[match_idx]
  estimates$dpar = pars$dpar[match_idx]
  estimates = estimates[order(match_idx, estimates$name), ]

  # Add simulation parameters if the data is simulated
  sim_list = attr(fit$data[, fit$pars$y], "simulated")
  if(!is.null(sim_list)) {
    simulated = as.list(sim_list)  # Get as oroper list
    simulated = simulated[sapply(simulated, is.numeric)]  # Remove non-numeric

    # Handle varying effects. Finds the matching labels
    for (this_varying in fit$pars$varying) {
      if (!is.null(simulated[[this_varying]])) {
        # Find the needed values and labels
        value = simulated[[this_varying]]  # Extract simulation values
        group_effects = get_fit_model_tables(fit)$group_effects
        label_col = group_effects$group_col[group_effects$name == this_varying]
        labs = fit$data[[label_col]]  # Find the labels. Same length as `value`
        if (length(value) != length(labs)) {
          warning("This is simulated data, but the labels for varying effect '", label_col, "' in data does not have the same length as the numeric params used for simulation.")
          next
        }

        # Name like the MCMC columns and use only unique combinations (assuming identical value for each level)
        value = value[!duplicated(value)]
        names(value) = unique(paste0(this_varying, "[", labs, "]"))

        # Delete the simulation vector and add the new label-value pairs to list
        simulated[[this_varying]] = NULL
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

    # Add simulation to the beginning of the list
    estimates = estimates %>%
      dplyr::left_join(simulated, by = "name", relationship = "one-to-one") %>%
      dplyr::mutate(match = ifelse(.data$sim > .data$lower & .data$sim < .data$upper, yes = "OK", no = "")) %>%
      dplyr::relocate("name", "segment", "dpar", "match", "sim")
  } else {
    estimates = dplyr::relocate(estimates, "name", "segment", "dpar")
  }

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
#' @param digits a non-null value for digits specifies the minimum number of
#'   significant digits to be printed in values. The default, NULL, uses
#'   getOption("digits"). (For the interpretation for complex numbers see signif.)
#'   Non-integer values will be rounded down, and only values greater than or
#'   equal to 1 and no greater than 22 are accepted.
#' @param prior TRUE/FALSE. Summarise prior instead of posterior?
#' @param ... Currently ignored
#'
#' @return A data frame with parameter estimates and MCMC diagnostics. Rows
#'   are ordered by change point first, then `mu`, then the other
#'   distributional parameters, then `ar`/`ma` components - each ascending by
#'   segment. OBS: The change point distributions are often not unimodal and
#'   symmetric so the intervals can be deceiving. Plot them using
#'   `plot_pars(fit)`.
#'
#'   * `segment` is the segment the parameter belongs to.
#'   * `dpar` is the distributional parameter (`"cp"`, `"mu"`, `"sigma"`,
#'     `"ar1"`, `"ma1"`, etc.) the parameter belongs to.
#'   * `mean` is the posterior mean
#'   * `lower` and `upper` are the bounds of the central posterior interval
#'     given in `width`.
#'   * `Rhat` is the rank-normalized split-Rhat convergence diagnostic.
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
#' # Varying (random) effects
#' # ranef(my_fit)
#'
#' # Summarise prior
#' summary(demo_fit, prior = TRUE)
summary.mcpfit = function(object, width = 0.95, digits = 2, prior = FALSE, ...) {
  fit = object  # Standard name in mcp
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_number(width, lower = 0, upper = 1)
  checkmate::assert_int(digits, lower = 0)
  checkmate::assert_flag(prior)
  rlang::check_dots_empty()

  samples = mcmclist_samples(fit, prior = prior, error = FALSE)

  # Model info
  cat("Family: ", fit$family$family, "(link = '", fit$family$link, "')\n", sep = "")
  if (!is.null(samples))
    cat("Iterations: ", niterations(fit, prior = prior), " from ", nchains(fit, prior = prior), " chains.\n", sep="")
  cat("Segments:\n")
  for (i in 1:length(fit$model)) {
    cat("  ", i, ": ", formula_to_char(fit$model[[i]]), "\n", sep = "")
  }

  # Data
  if (!is.null(samples)) {
    # Print and return invisibly
    cat("\nPopulation-level parameters:\n")
    result = get_summary(fit, width, varying = FALSE, prior = prior)
    print(data.frame(result), digits = digits, row.names = FALSE)

    if (!is.null(fit$pars$varying))
      cat("\nUse `ranef(fit)` to summarise the varying effect(s):", paste0(fit$pars$varying, collapse = ", "))

    return(invisible(result))
  }
  else {
    cat("\nNo samples. Nothing to summarise.")
    return(invisible(NULL))
  }
}



#' @export
fixef = function(object, ...) UseMethod("fixef")

#' @export
ranef = function(object, ...) UseMethod("ranef")


#' @aliases fixef fixef.mcpfit
#' @describeIn summary.mcpfit Fixed (population-level) effects of `mcpfit`.
#' @export
fixef.mcpfit = function(object, width = 0.95, prior = FALSE, ...) {
  rlang::check_dots_empty()
  get_summary(object, width, varying = FALSE, prior = prior)
}

#' @aliases ranef ranef.mcpfit
#' @describeIn summary.mcpfit Random (varying) effects of `mcpfit`.
#' @export
ranef.mcpfit = function(object, width = 0.95, prior = FALSE, ...) {
  rlang::check_dots_empty()
  get_summary(object, width, varying = TRUE, prior = prior)
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


#' Check if this is an AR/MA model
#'
#' @aliases is_arma
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return TRUE or FALSE
is_arma = function(fit) {
  length(fit$pars$arma) > 0
}


#' Internal function to get samples.
#'
#' Returns posterior samples, if available. If not, then prior samples. If not,
#' then throw an informative error. This is useful for summary and plotting, that
#' works on both.
#'
#' @aliases mcmclist_samples mcmclist_samples.mcpfit
#' @keywords internal
#' @inheritParams summary.mcpfit
#' @param fit An \code{\link{mcpfit}} object
#' @param message TRUE: gives a message if returning prior samples. FALSE = no message
#' @param error TRUE: err if there are no samples. FALSE: return NULL
mcmclist_samples = function(fit, prior = FALSE, message = TRUE, error = TRUE) {
  if (prior == TRUE) {
    if (coda::is.mcmc.list(fit$mcmc_prior)) {
      return(fit$mcmc_prior)
    } else {
      stop("Prior requested but the prior was not sampled.")
    }
  }

  if (coda::is.mcmc.list(fit$mcmc_post)) {
    return(fit$mcmc_post)
  } else if (coda::is.mcmc.list(fit$mcmc_prior)) {
    if (message)
      message("Posterior was not sampled. Using prior samples. Set `prior = TRUE` to mute this message.")
    return(fit$mcmc_prior)
  } else if (error == TRUE) {
    stop("This mcpfit contains no posterior or prior samples.")
  }

  NULL
}


#' Get samples as a posterior draws array
#'
#' This is the single internal conversion from the stored
#' \code{\link[coda]{mcmc.list}} representation to a posterior draws object.
#'
#' @keywords internal
#' @noRd
#' @inheritParams mcmclist_samples
posterior_draws = function(fit, prior = FALSE, message = TRUE, error = TRUE) {
  samples = mcmclist_samples(
    fit,
    prior = prior,
    message = message,
    error = error
  )
  if (is.null(samples))
    return(NULL)

  posterior::as_draws_array(samples)
}


#' Index \code{mcpfit} objects
#'
#' Index variables, iterations, chains, and draws.
#'
#' @inheritParams fitted.mcpfit
#' @name draws-index-mcp
#' @examples
#' niterations(demo_fit)
#' nchains(demo_fit, prior = TRUE)
NULL


#' @aliases niterations niterations.mcpfit
#' @describeIn draws-index-mcp Total number of iterations of an `mcpfit` object.
#' @export
niterations.mcpfit = function(object, prior = FALSE, ...) {
  samples = mcmclist_samples(object, prior = prior, error = FALSE)
  sum(sapply(samples, nrow))
}

#' @aliases nchains nchains.mcpfit
#' @describeIn draws-index-mcp Number of chains of an `mcpfit` object.
#' @export
nchains.mcpfit = function(object, prior = FALSE, ...) {
  samples = mcmclist_samples(object, prior = prior, error = FALSE)
  length(samples)
}

#' @export
niterations = function(object, ...) UseMethod("niterations")

#' @export
nchains = function(object, ...) UseMethod("nchains")


#' Get relevant info about varying parameters
#'
#' Returns parameters, data columns, and effect metadata given parameter
#' name(s), model part(s), or column(s).
#'
#' @aliases unpack_varying unpack_varying.mcpfit
#' @keywords internal
#' @noRd
#' @param pars `NULL`/`FALSE` for nothing. `TRUE` for all. A character vector
#'   containing `"cp"`, `"predictor"`, or exact varying parameter names.
#' @param cols `NULL`/`FALSE` for nothing. `TRUE` for all. A vector of varying column names for specifics. Usually provided via "facet_by" argument in other functions.
#' @return A list. See details.
#'
#' @details
#' Returns a list with
#' @slot pars Character vector of parameter names. `NULL` if empty.
#' @slot cols Character vector of data column names. `NULL` if empty.
#' @slot indices Logical vector indexing the group-effects table.
#' @slot effects The selected rows of the group-effects table.
unpack_varying = function(fit, pars = NULL, cols = NULL) {
  checkmate::assert_multi_class(pars, c("logical", "character"), null.ok = TRUE)
  checkmate::assert_multi_class(cols, c("logical", "character"), null.ok = TRUE)
  if (is.logical(pars))
    checkmate::assert_flag(pars)
  if (is.logical(cols))
    checkmate::assert_flag(cols)
  group_effects = get_fit_model_tables(fit)$group_effects
  use_varying = rep(FALSE, nrow(group_effects))

  # If everything is NULL, just return NULLs
  if ((is.null(pars) && is.null(cols))) {
    return(list(
      pars = NULL,
      cols = NULL,
      indices = use_varying,
      effects = group_effects[use_varying, , drop = FALSE]
    ))
  } else if (!is.null(pars) && !is.null(cols)) {
    stop("One of `pars` and `cols` must be NULL.")
  }


  if (!is.null(pars)) {
    if (all(pars == FALSE)) {
      # Select no varying effects
      use_varying[] = FALSE
    } else if (all(pars == TRUE)) {
      # Select all varying effects
      use_varying[] = TRUE
    } else if (is.character(pars)) {
      allowed = c("cp", "predictor", group_effects$name)
      unknown = pars[pars %notin% allowed]
      if (length(unknown) > 0)
        stop(
          "Unknown `varying` selection: ", and_collapse(unknown), ". ",
          "Use TRUE, FALSE, \"cp\", \"predictor\", or names from fit$pars$varying."
        )
      use_varying = group_effects$part %in% pars | group_effects$name %in% pars
    }
  } else if (!is.null(cols)) {
    if (all(cols == TRUE)) {
      use_varying[] = TRUE
    } else if (!all(cols == FALSE)) {
      use_varying = group_effects$group_col %in% cols
    }
  }

  # Return
  list(
    pars = logical0_to_null(group_effects$name[use_varying]),
    cols = logical0_to_null(group_effects$group_col[use_varying]),
    indices = use_varying,
    effects = group_effects[use_varying, , drop = FALSE]
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


#' Get tidy samples with or without varying effects
#'
#' Returns in a format useful for `fit$simulate()` with population parameters in wide format
#' and varying effects in long format (the number of rows will be `ndraws * n_levels_in_varying`).
#'
#' @aliases tidy_samples tidy_samples.mcpfit
#' @keywords internal
#' @inheritParams mcmclist_samples
#' @inheritParams pp_eval
#' @param population
#'   * `TRUE` All population effects. Same as `fit$pars$population`.
#'   * `FALSE` No population effects. Same as `c()`.
#'   * Character vector: Only include specified population parameters - see `fit$pars$population`.
#' @param varying One of:
#'   * `TRUE` All varying effects (`fit$pars$varying`).
#'   * `FALSE` No varying effects (`c()`).
#'   * `"cp"` or `"predictor"`: All varying effects belonging to that part of
#'     the model.
#'   * Character vector: Only include specified varying parameters - see
#'     `fit$pars$varying`.
#' @param absolute
#'   * `TRUE` Returns the absolute location of all varying change points.
#'   * `FALSE` Just returns the varying effects.
#'   * Character vector: Only do absolute transform for these varying parameters - see `fit$pars$varying`.
#'
#' @return `tibble` of posterior draws in `tidybayes` format.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
tidy_samples = function(
  fit,
  population = TRUE,
  varying = TRUE,
  absolute = FALSE,
  prior = FALSE,
  ndraws = NULL,
  nsamples = lifecycle::deprecated()
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "tidy_samples")

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
  # Varying parameters. Result is `terms_varying`
  varying_info = unpack_varying(fit, pars = varying)
  terms_varying = paste0(varying_info$pars, "[", varying_info$cols, "]")  # for tidybayes
  if (all(terms_varying == "[]")) terms_varying = ""  # quick fix

  # Population parameters. Result is `pars_population`.
  if (all(population == FALSE)) {
    pars_population = c()  # Empty if no absolute varying
  } else if (all(population == TRUE)) {
    pars_population = fit$pars$population
  } else if (is.character(population)) {
    if (!all(population %in% fit$pars$population))
      stop("Not all `population` are population parameters (see fit$pars$population).")

    pars_population = population
  }

  # Absolute effects. Results is `absolute_cps` and `absolute` (recoded to varying cp names)
  if (all(absolute == TRUE)) {
    cp_effects = varying_info$effects[varying_info$effects$part == "cp", , drop = FALSE]
    absolute = cp_effects$name
    absolute_cps = cp_effects$population_name
  } else if (all(absolute == FALSE)) {
    absolute_cps = NULL
  } else if (is.character(absolute)) {
    # Check
    is_in_varying = absolute %in% varying_info$pars
    if (any(!is_in_varying))
      stop("The following parameter names in `absolute` are not in `varying`: ", and_collapse(absolute[!is_in_varying]))
    absolute_effects = varying_info$effects[
      varying_info$effects$name %in% absolute, , drop = FALSE
    ]
    if (any(absolute_effects$part != "cp"))
      stop("`absolute` can select change-point group-level effects only.")

    absolute_cps = absolute_effects$population_name
  }

  # ----- GET THESE PARAMETERS AS TIDY DRAWS -----
  # Select posterior/prior samples
  samples = mcmclist_samples(fit, prior = prior)

  # Build code for tidybayes::spread_draws() and execute it
  all_terms = unique(c(pars_population, terms_varying, absolute_cps))
  code = paste0("tidybayes::spread_draws(samples, ", paste0(all_terms, collapse = ", "), ", ndraws = ndraws)")
  samples = eval(str2lang(code))

  # Make varying columns factor if they are factors in fit$data
  if (length(varying_info$cols) > 0) {
    is_factor = lapply(fit$data, is.factor)[varying_info$cols]
    cols_to_factorize = varying_info$cols[as.logical(is_factor)]
    samples = dplyr::mutate_at(samples, cols_to_factorize, as.factor)
  }

  # Add population cp to varying and delete population cols only included for this purpose.
  if (length(absolute_cps) > 0) {
    #samples[, absolute] = samples[, absolute] + samples[, absolute_cps]
    samples[, absolute_cps] = samples[, absolute_cps] + samples[, absolute]
    samples = dplyr::select(samples, -dplyr::all_of(absolute))
  }

  # Unassigned varying effects are just simulated as zero (the population mean)
  remaining_varying_cols = dplyr::setdiff(fit$pars$varying, colnames(samples))
  samples[, remaining_varying_cols] = 0

  # Return with chain etc. first
  samples %>%
    dplyr::relocate(".chain", ".iteration", ".draw")
}




#' Fits and Predictions given Draws and data
#'
#' @aliases pp_eval pp_eval.mcpfit
#' @keywords internal
#' @inheritParams tidy_samples
#' @param object An `mcpfit` object.
#' @param newdata A `tibble` or a `data.frame` containing predictors in the model. Weighted
#'   Gaussian predictions and log-likelihoods also require the weights column. If `NULL`
#'   (default), the original data is used.
#' @param summary Summarise at each x-value
#' @param type One of:
#'   - `"fitted"`: return expected values. When `dpar` is the name of a dpar
#'     (e.g., `"mu"` or `"sigma"`), the expected value for just this dpar is returned.
#'     See also `fitted()`.
#'   - `"predict"`: return predicted values (e.g., `y_predict = rnorm(N, y_fitted, sigma_fitted)` for `family = gaussian()`).
#'     See also `predict()`.
#'   - `"residuals"`: observed y-values minus the fitted values. See also `residuals()`.
#'   - `"loglik"`: return the log-likelihood for each sample for each data point. See also `log_lik()`.
#'     Requires `scale = "response"`.
#' @param probs Vector of quantiles. Only in effect when `summary == TRUE`.
#' @param rate Boolean. For binomial models, plot on raw data (`rate = FALSE`) or
#'   response divided by number of trials (`rate = TRUE`). If FALSE, linear
#'   interpolation on trial number is used to infer trials at a particular x.
#' @param prior TRUE/FALSE. Plot using prior samples? Useful for `mcp(..., sample = "both")`
#' @param dpar What distributional parameter to evaluate. This is only relevant when `type == "fitted"`. E.g.,
#'
#'   * `"epred"` (default): Expected value of the full model (or `NULL` for compatibility with brms etc.).
#'   * `"mu"`: The central tendency which is often the mean after applying the
#'     link function.
#'   * `"sigma"`: The standard deviation of the residuals.
#'   * `"ar1"`, `"ar2"`, `"ma1"`, `"ma2"`, etc. depending on which AR or MA
#'     coefficient you want to evaluate.
#' @param arma Whether to include AR and MA effects.
#'   * `TRUE` Compute the GARMA residual recurrence. Requires the response variable in `newdata`.
#'   * `FALSE` Disregard AR and MA effects. For `family = gaussian()`, `predict()` uses only `sigma` for residuals.
#' @param ndraws Integer or `NULL`. Number of posterior draws to return/summarise.
#'   If there are varying effects, this is the number of draws from each varying group.
#'   `NULL` means "all". Ignored if both are `FALSE`. More samples trade speed for accuracy.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @param samples_format One of "tidy" or "matrix". Controls the output format when `summary == FALSE`.
#'   See more under "value"
#' @param scale One of
#'   * `"response"`: return on the observed scale, i.e., after applying the inverse link function.
#'   * `"linear"`: return on the parameter scale (where the linear trends are modelled).
#'     A linear scale is only applicable when `type == "fitted"` and `dpar` is not `NULL`.
#' @param .include_fitted Internal. Include fitted values with unsummarised predictions.
#' @return
#'   * If `summary = TRUE`: A `tibble` with the posterior mean for each row in `newdata`,
#'     If `newdata` is `NULL`, the data in `fit$data` is used.
#'
#'   * If `summary = FALSE` and `samples_format = "tidy"`: A `tidybayes` `tibble` with all the posterior
#'     draws (`Nd`) evaluated at each row in `newdata` (`Nn`), i.e., with `Nd x Nn` rows. If there are
#'     varying effects, the returned data is expanded with the relevant levels for each row.
#'
#'     The return columns are:
#'
#'      - Predictors from `newdata`.
#'      - Draw descriptors: ".chain", ".iteration", ".draw" (see the `posterior` and `tidybayes` packages), and `data_row`, the row number in the evaluated `newdata`.
#'      - Draw values: one column for each parameter in the model.
#'      - The estimate. Either "predict" or "fitted", i.e., the name of the `type` argument.
#'
#'   * If `summary = FALSE` and `samples_format = "matrix"`: An `N_draws` X `nrows(newdata)` matrix with fitted/predicted
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
  samples_format = "tidy",
  scale = 'response',
  .include_fitted = FALSE,
  nsamples = lifecycle::deprecated()
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "pp_eval")

  # Recode
  fit = object
  checkmate::assert_class(fit, "mcpfit")
  if (!is.mcpfamily(fit$family))
    fit$family = mcpfamily(fit$family)
  dpar = assert_dpar(dpar, fit = fit, type = type)

  if (is.null(newdata))
    newdata = fit$data


  ###############
  # FIX NEWDATA #
  ###############
  varying_info = unpack_varying(fit, pars = varying)
  model_tables = get_fit_model_tables(fit)
  varying_cols = unique(stats::na.omit(model_tables$group_effects$group_col))
  exclude_varying = setdiff(varying_cols, varying_info$cols)
  required_cols = colnames(fit$data)  # Only predictive columns were saved in fit$data
  operation = switch(type, predict = "rng", loglik = "log_lik", fitted = "epred", residuals = "epred")
  aux_operations = c(operation, if (arma && is_arma(fit)) "garma")
  aux_columns = get_family_aux_columns(fit$family, model_tables$segments)
  aux_used = names(get_family_aux_columns(fit$family, model_tables$segments, aux_operations))
  unused_aux_columns = unname(aux_columns[names(aux_columns) %notin% aux_used])
  required_cols = required_cols[required_cols %notin% unused_aux_columns]
  required_cols = required_cols[required_cols %notin% exclude_varying]
  if ((arma == FALSE || is_arma(fit) == FALSE) & type %in% c("fitted", "predict")) {
    required_cols = required_cols[required_cols != fit$pars$y]
  } else if (fit$pars$y %notin% colnames(newdata)) {
    stop("`newdata` must contain a response column named '", fit$pars$y, "' for when `arma == TRUE` and/or `type == 'residuals'`")
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
  samples_format = rlang::arg_match0(samples_format, c("tidy", "matrix"))


  ########################
  # GET FITS/PREDICTIONS #
  ########################
  simulate_type = ifelse(type == "residuals", yes = "fitted", no = type)
  if (length(varying_info$cols) > 0) {
    # If there are varying effects: use varying-matching samples for each row of data
    samples_predictors = dplyr::left_join(
      add_rhs_predictors(newdata, fit),
      tidy_samples(fit, population = TRUE, varying = varying, prior = prior, ndraws = ndraws),
      by = unique(varying_info$cols),
      relationship = "many-to-many"
    )
  } else {
    # No varying effects: use all samples for each row of data
    samples_predictors = tidy_samples(fit, population = TRUE, varying = varying, prior = prior, ndraws = ndraws) %>%
      tidyr::expand_grid(add_rhs_predictors(newdata, fit))
  }

  samples = samples_predictors %>%
    dplyr::mutate("{type}" := rlang::exec(simulate_vectorized, fit, !!!samples_predictors, .type = simulate_type, .rate = rate, .dpar = dpar, .arma = arma, .scale = scale))

  # Plotting can request fitted and predicted values from the same evaluated
  # parameter rows, guaranteeing identical joint draw IDs without rebuilding
  # model predictors in the plotting layer.
  if (.include_fitted) {
    samples$fitted = rlang::exec(simulate_vectorized, fit, !!!samples_predictors, .type = "fitted", .rate = rate, .dpar = dpar, .arma = arma, .scale = scale)
  }

  samples = samples %>% dplyr::select(-dplyr::starts_with(".pred_"), -dplyr::any_of(point_size_col))

  # Missing outcomes are latent in the fitted JAGS model, but they are not
  # observed-data likelihood contributions. Retain them while evaluating
  # GARMA histories above, then remove them from returned log likelihoods.
  if (type == "loglik") {
    observed_rows = which(!is.na(newdata[, fit$pars$y]))
    if (length(observed_rows) == 0)
      stop("Log-likelihood evaluation requires at least one observed response.")
    samples = dplyr::filter(samples, .data$data_row %in% observed_rows)
    newdata_return = dplyr::filter(
      newdata_return,
      .data$data_row %in% observed_rows
    )
  }


  # Optionally compute residuals
  if (type == "residuals")
    samples = dplyr::mutate(samples, !!type := .data[[fit$pars$y]] - .data[[type]])

  # Fail early if varying-effect joins or another evaluation step duplicated
  # or dropped any joint draw/evaluation-row combinations.
  validate_eval_draws(samples, type)

  # Optionally summarise
  if (summary == TRUE) {
    df_return = samples %>%
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
      quantiles = get_quantiles(samples, probs, type) %>%
        dplyr::mutate(quantile = 100 * .data$quantile) %>%
        tidyr::pivot_wider(names_from = "quantile", names_prefix = "Q", values_from = dplyr::all_of(type))

      df_return = dplyr::left_join(df_return, quantiles, by = "data_row", relationship = "one-to-one")
    }
    return(data.frame(dplyr::select(df_return, -"data_row")))
  } else if (samples_format == "tidy") {
    return(samples)
  } else if (samples_format == "matrix") {
    df_return = tidy_to_matrix(samples, type)
    return(df_return)
  }
}



#' Fitted and predicted values of `mcp` models fits
#'
#' Evaluate the model on data, either summarised (per data-row) or per draw. You
#' can use draws from the prior (`prior = TRUE`), select the parameter to predict
#' from (``)
#'
#' @details
#' `residuals(fit)` is equivalent to `fit$data[, fit$data$yvar] - fitted(fit, ...)` (or `newdata[, fit$data$yvar] - fitted(fit, ...)`),
#' but with fixed arguments for `fitted`: `rate = FALSE, dpar = 'epred', samples_format = 'tidy'`.
#'
#' @inheritParams pp_eval
#' @param ... Currently ignored.
#' @inherit pp_eval return
#' @seealso \code{\link{fitted.mcpfit}} \code{\link{predict.mcpfit}} \code{\link{residuals.mcpfit}} \code{\link{log_lik.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' fitted(demo_fit)  # Expected value at each demo_fit$data at response-level
#' residuals(demo_fit)  # Residuals at each demo_fit$data at response-level
#' log_lik(demo_fit)  # Log-likelihood at each demo_fit$data
#'
#' # All of the above take a range of arguments. E.g.,:
#' \donttest{
#' predict(demo_fit)  # Pointwise posterior predictive
#' predict(demo_fit, probs = c(0.1, 0.5, 0.9))  # With median and 80% credible interval.
#' predict(demo_fit, prior = TRUE)  # Prior predictive
#' fitted(demo_fit, summary = FALSE)  # Samples instead of summary. Useful for plotting distributions.
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
  samples_format = "tidy",
  nsamples = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "predict.mcpfit")
  rlang::check_dots_empty()
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
    samples_format = samples_format
  )
}


#' @aliases fitted fitted.mcpfit
#' @describeIn execute-mcp-model Expected distribution
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
  samples_format = "tidy",
  scale = "response",
  nsamples = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "fitted.mcpfit")
  rlang::check_dots_empty()
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
    samples_format = samples_format,
    scale = scale
  )
}


#' @export
log_lik = function(object, ...) UseMethod("log_lik")

#' @aliases log_lik log_lik.mcpfit
#' @describeIn execute-mcp-model Pointwise log-likelihood
#' @export
log_lik.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  ndraws = NULL,
  samples_format = "tidy",
  nsamples = lifecycle::deprecated(),
  ...
) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "log_lik.mcpfit")
  rlang::check_dots_empty()
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
    samples_format = samples_format,
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
  rlang::check_dots_empty()
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
    samples_format = "tidy"
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
