# ABOUT: These are non-plotting functions that take an mcpfit as the first argument
# -----------------

#' Class `mcpfit` of models fitted with the \pkg{mcp} package
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
#' @slot mcmc_loglik An \code{\link[coda]{mcmc.list}} object with samples of log-likelihood.
#' @slot pars A list of character vectors of model parameter names.
#' @slot jags_code A string with jags code. Use `cat(fit$jags_code)` to show it.
#' @slot simulate A method to simulate and predict data.
#' @slot .other Information that is used internally by mcp.
#'
NULL


#' Internal function for summary.mcpfit, fixef.mcpfit, and ranef.mcpfit
#'
#' @aliases get_summary get_summary.mcpfit
#' @keywords internal
#' @inheritParams summary.mcpfit
#' @param fit An \code{\link{mcpfit}}` object.
#' @param varying Boolean. Get results for varying (TRUE) or population (FALSE)?
#' @return A data.frame with summaries for each model parameter.
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_summary = function(fit, width, varying = FALSE, prior = FALSE) {
  # Check arguments
  assert_mcpfit(fit)
  assert_numeric(width, lower = 0, upper = 1)
  assert_logical(varying)
  assert_logical(prior)

  # Get posterior/prior samples
  # TO DO: transition this to tidy_samples at some point if rhat and n.eff
  #    can be computed on this data. Would avoid some hacky solutions in the code below.
  samples = mcmclist_samples(fit, prior = prior)

  # Select only varying or only population-level columns in data
  if (varying == FALSE) {
    regex_pars = paste0(paste0("^", fit$pars$population, "$"), collapse = "|")
  } else {
    regex_pars = paste0("^", paste(fit$pars$varying, collapse="|^"))
    if (regex_pars == "^")
      stop("There were no matching parameters in the model.")
  }

  # HACK: If there is just one parameter, it's not a matrix. Add two in to make the code run.
  all_cols = colnames(samples[[1]])
  get_cols = all_cols[stringr::str_detect(all_cols, regex_pars)]
  if (!stringr::str_detect(regex_pars, "\\|") && varying == FALSE) {
    samples = lapply(samples, function(x) x[, c(get_cols, "cp_0", "cp_1")])
  } else {
    samples = lapply(samples, function(x) x[, get_cols])
  }
  class(samples) = "mcmc.list"  # Return to original class

  # Get parameter estimates
  estimates = samples %>%
    # Get ready to compute stuff for each parameter
    tidybayes::tidy_draws() %>%
    tidyr::pivot_longer(-tidyselect::starts_with(".")) %>%
    dplyr::group_by(.data$name) %>%

    # Compute mean and HDI intervals and name appropriately
    tidybayes::mean_hdci(.data$value, .width = width) %>%
    dplyr::rename(mean = .data$value,
                  lower = .data$.lower,
                  upper = .data$.upper) %>%

    # Remove unneeded stuf
    dplyr::select(-tidyselect::starts_with("."))

  # Revert HACK and continue
  if (!stringr::str_detect(regex_pars, "\\|") && varying == FALSE) {
    estimates = dplyr::filter(estimates, !.data$name %in% c("cp_0", "cp_1"))
    samples = lapply(samples, function(x) x[, get_cols])
  }

  # Diagnostics: Gelman-Rubin and effective sample size
  Rhat = try(coda::gelman.diag(samples, multivariate = FALSE)$psrf[, 1], TRUE)
  if (!is.numeric(Rhat)) {
    warning("Rhat computation failed: ", Rhat)
    Rhat = rep(NA, nrow(estimates))
  }
  diagnostics = data.frame(
    Rhat = Rhat,
    n.eff = round(coda::effectiveSize(samples))
  )

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
        which_label_col = which(fit$.other$ST$cp_group == this_varying)
        label_col = fit$.other$ST$cp_group_col[which_label_col]  # What column is the labels in data.
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
      dplyr::select("name", "match", "sim", "mean", "lower", "upper")
  }

  # Merge them and return
  result = dplyr::bind_cols(estimates, diagnostics)
  data.frame(result, row.names = NULL)
}


#' Summarise mcpfit objects
#'
#' Summarise parameter estimates and model diagnostics.
#'
#' @aliases summary summary.mcpfit
#' @param object An \code{\link{mcpfit}} object.
#' @param width Float. The width of the highest posterior density interval
#'   (between 0 and 1).
#' @param digits a non-null value for digits specifies the minimum number of
#'   significant digits to be printed in values. The default, NULL, uses
#'   getOption("digits"). (For the interpretation for complex numbers see signif.)
#'   Non-integer values will be rounded down, and only values greater than or
#'   equal to 1 and no greater than 22 are accepted.
#' @param prior TRUE/FALSE. Summarise prior instead of posterior?
#' @param ... Currently ignored
#'
#' @return A data frame with parameter estimates and MCMC diagnostics.
#'   OBS: The change point distributions are often not unimodal and symmetric so
#'   the intervals can be deceiving Plot them using `plot_pars(fit)`.
#'
#'   * `mean` is the posterior mean
#'   * `lower` is the lower quantile of the highest-density interval (HDI) given in `width`.
#'   * `upper` is the upper quantile.
#'   * `Rhat` is the Gelman-Rubin convergence diagnostic which is often taken to
#'     be acceptable if < 1.1. It is computed using \code{\link[coda]{gelman.diag}}.
#'   * `n.eff` is the effective sample size computed using \code{\link[coda]{effectiveSize}}.
#'     Low effective sample sizes are also obvious as poor mixing in trace plots
#'     (see `plot_pars(fit)`). Read how to deal with such problems [here](https://lindeloev.github.io/mcp/articles/tips.html)
#'   * `ts_err` is the time-series error, taking autoregressive correlation
#'     into account. It is computed using \code{\link[coda]{spectrum0.ar}}.
#'
#'  For simulated data, the summary contains two additional columns so that it
#'  is easy to inspect whether the model can recover the parameters. Run
#'  simulation and summary multiple times to get a sense of the robustness.
#'
#'   * `sim` is the value used to generate the data.
#'   * `match` is `"OK"` if `sim` is contained in the HDI interval (`lower` to
#'     `upper`).
#'
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' # Typical usage
#' summary(demo_fit)
#' summary(demo_fit, width = 0.8, digits = 4)  # Set HDI width
#'
#' # Get the results as a data frame
#' results = summary(demo_fit)
#'
#' # Varying (random) effects
#' # ranef(my_fit)
#'
#' # Summarise prior
#' summary(demo_fit, prior = TRUE)
#'
summary.mcpfit = function(object, width = 0.95, digits = 2, prior = FALSE, ...) {
  fit = object  # Standard name in mcp
  assert_mcpfit(fit)
  assert_numeric(width, lower = 0, upper = 1)
  assert_integer(digits, lower = 0)
  assert_logical(prior)
  assert_ellipsis(...)

  samples = mcmclist_samples(fit, prior = prior, error = FALSE)

  # Model info
  cat("Family: ", fit$family$family, "(link = '", fit$family$link, "')\n", sep = "")
  if (!is.null(samples))
    cat("Iterations: ", nrow(samples[[1]]) * length(samples), " from ", length(samples), " chains.\n", sep="")
  cat("Segments:\n")
  for (i in 1:length(fit$model)) {
    cat("  ", i, ": ", format(fit$model[i]), "\n", sep = "")
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


#' @aliases fixef fixef.mcpfit fixed.effects
#' @describeIn summary.mcpfit Get population-level ("fixed") effects of an \code{\link{mcpfit}} object.
#' @export
fixef = function(object, width = 0.95, prior = FALSE, ...) {
  assert_ellipsis(...)
  get_summary(object, width, varying = FALSE, prior = prior)
}


#' @aliases ranef ranef.mcpfit random.effects
#' @describeIn summary.mcpfit Get varying ("random") effects of an \code{\link{mcpfit}} object.
#' @export
ranef = function(object, width = 0.95, prior = FALSE, ...) {
  assert_ellipsis(...)
  get_summary(object, width, varying = TRUE, prior = prior)
}


#' @aliases print print.mcpfit
#' @describeIn summary.mcpfit Print the posterior summary of an \code{\link{mcpfit}} object.
#' @param x An \code{\link{mcpfit}} object.
#' @export
print.mcpfit = function(x, ...) {
  summary(x, ...)
}


#' Checks if argument is an `mcpfit` object
#'
#' @aliases is.mcpfit
#' @param x An `R` object.
#' @export
#'
is.mcpfit = function(x) {
  inherits(x, "mcpfit")
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
    message("Posterior was not sampled. Using prior samples.")
    return(fit$mcmc_prior)
  } else if (error == TRUE) {
    stop("This mcpfit contains no posterior or prior samples.")
  }

  return(NULL)
}



#' Get relevant info about varying parameters
#'
#' Returns parameters, data columns, and implicated segments given parameter name(s) or column(s).
#'
#' @aliases unpack_varying unpack_varying.mcpfit
#' @keywords internal
#' @param pars `NULL`/`FALSE` for nothing. `TRUE` for all. A vector of varying parameter names for specifics.
#' @param cols `NULL`/`FALSE` for nothing. `TRUE` for all. A vector of varying column names for specifics. Usually provided via "facet_by" argument in other functions.
#' @return A list. See details.
#'
#' @details
#' Returns a list with
#' @slot pars Character vector of parameter names. `NULL` if empty.
#' @slot cols Character vector of data column names. `NULL` if empty.
#' @slot indices Logical vector of segments in the segment table that contains the varying effect
#'
unpack_varying = function(fit, pars = NULL, cols = NULL) {
  assert_types(pars, "null", "logical", "character")
  assert_types(cols, "null", "logical", "character")

  # If everything is NULL, just return NULLs
  if ((is.null(pars) && is.null(cols))) {
    return(list(
      pars = NULL,
      cols = NULL,
      indices = rep(FALSE, nrow(fit$.other$ST))
    ))
  } else if (!is.null(pars) && !is.null(cols)) {
    stop("One of `pars` and `cols` must be NULL.")
  }


  if (!is.null(pars)) {
    if (all(pars == FALSE)) {
      # Select no varying effects
      use_varying = NULL
    } else if (all(pars == TRUE)) {
      # Select all varying effects
      use_varying = !is.na(fit$.other$ST$cp_group)
    } else if (is.character(pars)) {
      if (!all(pars %in% fit$pars$varying))
        stop("Not all `pars` are varying parameters (see fit$pars$varying).")
      # Select only specified varying effects
      use_varying = fit$.other$ST$cp_group %in% pars
    }
  } else if (!is.null(cols)) {
    use_varying = tidyr::replace_na(fit$.other$ST$cp_group_col == cols, FALSE)
  }

  return(list(
    pars = fit$.other$ST$cp_group[use_varying],
    cols = fit$.other$ST$cp_group_col[use_varying],
    indices = use_varying
  ))
}



#' Get tidy samples with or without varying effects
#'
#' Returns in a format useful for `fit$simulate()` with population parameters in wide format
#' and varying effects in long format (the number of rows will be `nsamples * n_levels_in_varying`).
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
#'   * Character vector: Only include specified varying parameters - see `fit$pars$varying`.
#' @param absolute
#'   * `TRUE` Returns the absolute location of all varying change points.
#'   * `FALSE` Just returns the varying effects.
#'   * Character vector: Only do absolute transform for these varying parameters - see `fit$pars$varying`.
#'
#'   OBS: This currently only applies to varying change points. It is not implemented for `rel()` regressors yet.
#' @return `tibble` of posterior draws in `tidybayes` format.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
tidy_samples = function(
  fit,
  population = TRUE,
  varying = TRUE,
  absolute = FALSE,
  prior = FALSE,
  nsamples = NULL
) {
  # General argument checks
  assert_mcpfit(fit)
  assert_types(population, "logical", "character")
  assert_types(varying, "null", "logical", "character")
  assert_types(absolute, "null", "logical", "character")
  assert_logical(prior)
  if (!is.null(nsamples))
    assert_integer(nsamples, lower = 1)

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
    absolute = fit$.other$ST$cp_group[varying_info$indices]
    absolute_cps = fit$.other$ST$cp_name[varying_info$indices]
  } else if (all(absolute == FALSE)) {
    absolute_cps = NULL
  } else if (is.character(absolute)) {
    # Check
    is_in_varying = absolute %in% varying_info$pars
    if (any(!is_in_varying))
      stop("The following parameter names in `absolute` are not in `varying`:", absolute[is_in_varying])

    use_absolute = fit$.other$ST$cp_group %in% absolute
    absolute_cps = fit$.other$ST$cp_name[use_absolute]
  }

  # ----- GET THESE PARAMETERS AS TIDY DRAWS -----
  # Select posterior/prior samples
  samples = mcmclist_samples(fit, prior = prior)

  # Build code for tidybayes::spread_draws() and execute it
  all_terms = unique(c(pars_population, terms_varying, absolute_cps))
  code = paste0("tidybayes::spread_draws(samples, ", paste0(all_terms, collapse = ", "), ", ndraws = nsamples)")
  samples = eval(parse(text = code))

  # Make varying columns factor if they are factors in fit$data
  if (length(varying_info$cols) > 0) {
    is_factor = lapply(fit$data, is.factor)[varying_info$cols]
    cols_to_factorize = varying_info$cols[as.logical(is_factor)]
    samples = dplyr::mutate_at(samples, cols_to_factorize, as.factor)
  }

  # Add population cp to varying and delete population cols only included for this purpose.
  if (length(absolute_cps) > 0) {
    samples[, absolute] = samples[, absolute] + samples[, absolute_cps]
    samples = dplyr::select(samples, -!!dplyr::setdiff(absolute_cps, pars_population))
  }

  # Return with chain etc. first
  samples = dplyr::select(samples, ".chain", ".iteration", ".draw", dplyr::everything())
  return(samples)
}




#' Fits and predictions from samples and newdata
#'
#' @aliases pp_eval pp_eval.mcpfit
#' @keywords internal
#' @inheritParams tidy_samples
#' @param object An `mcpfit` object.
#' @param newdata A `tibble` or a `data.frame` containing predictors in the model. If `NULL` (default),
#'   the original data is used.
#' @param summary Summarise at each x-value
#' @param type One of:
#'   - "fitted": return fitted values. See also `fitted()`
#'   - "predict": return predicted values, using random dispersion around the central tendency
#'     (e.g., `y_predict = rnorm(N, y_fitted, sigma_fitted)` for `family = gaussian()`).
#'     See also `predict()`.
#'   - "residuals": same as "predict" but the observed y-values are subtracted. See also `residuals()`
#' @param probs Vector of quantiles. Only in effect when `summary == TRUE`.
#' @param rate Boolean. For binomial models, plot on raw data (`rate = FALSE`) or
#'   response divided by number of trials (`rate = TRUE`). If FALSE, linear
#'   interpolation on trial number is used to infer trials at a particular x.
#' @param prior TRUE/FALSE. Plot using prior samples? Useful for `mcp(..., sample = "both")`
#' @param which_y What to plot on the y-axis. One of
#'
#'   * `"ct"`: The central tendency which is often the mean after applying the
#'     link function.
#'   * `"sigma"`: The variance
#'   * `"ar1"`, `"ar2"`, etc. depending on which order of the autoregressive
#'     effects you want to plot.
#' @param arma Whether to include autoregressive effects.
#'   * `TRUE` Compute autoregressive residuals. Requires the response variable in `newdata`.
#'   * `FALSE` Disregard the autoregressive effects. For `family = gaussian()`, `predict()` just use `sigma` for residuals.
#' @param nsamples Integer or `NULL`. Number of samples to return/summarise.
#'   If there are varying effects, this is the number of samples from each varying group.
#'   `NULL` means "all". Ignored if both are `FALSE`. More samples trade speed for accuracy.
#' @param samples_format One of "tidy" or "matrix". Controls the output format when `summary == FALSE`.
#'   See more under "value"
#' @param scale One of
#'   * "response": return on the observed scale, i.e., after applying the inverse link function.
#'   * "linear": return on the parameter scale (where the linear trends are modelled).
#' @param ... Currently unused
#' @return
#'   * If `summary = TRUE`: A `tibble` with the posterior mean for each row in `newdata`,
#'     If `newdata` is `NULL`, the data in `fit$data` is used.
#'
#'   * If `summary = FALSE` and `samples_format = "tidy"`: A `tidybayes` `tibble` with all the posterior
#'     samples (`Ns`) evaluated at each row in `newdata` (`Nn`), i.e., with `Ns x Nn` rows. If there are
#'     varying effects, the returned data is expanded with the relevant levels for each row.
#'
#'     The return columns are:
#'
#'      - Predictors from `newdata`.
#'      - Sample descriptors: ".chain", ".iter", ".draw" (see the `tidybayes` package for more), and "data_row" (`newdata` rownumber)
#'      - Sample values: one column for each parameter in the model.
#'      - The estimate. Either "predict" or "fitted", i.e., the name of the `type` argument.
#'
#'   * If `summary = FALSE` and `samples_format = "matrix"`: An `N_draws` X `nrows(newdata)` matrix with fitted/predicted
#'       values (depending on `type`). This format is used by `brms` and it's useful as `yrep` in
#'      `bayesplot::ppc_*` functions.
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @seealso \code{\link{fitted.mcpfit}} \code{\link{predict.mcpfit}} \code{\link{residuals.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
pp_eval = function(
  object,
  newdata = NULL,
  summary = TRUE,
  type = "fitted",
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  which_y = "ct",
  varying = TRUE,
  arma = TRUE,
  nsamples = NULL,
  samples_format = "tidy",
  scale = 'response',
  ...
) {
  # Recodings
  fit = object
  assert_mcpfit(fit)
  xvar = rlang::sym(fit$pars$x)
  yvar = rlang::sym(fit$pars$y)
  returnvar = rlang::sym(type)
  is_arma = length(fit$pars$arma) > 0

  # What varying cols to use (varying = TRUE keeps them all as is)
  varying_info = unpack_varying(fit, pars = varying)
  required_cols = unique(c(fit$pars$x, fit$pars$trials, varying_info$cols))  # May be amended later

  # R CMD Check wants a global definition of ".". The formal way of doing it is
  # if(getRversion() >= "2.15.1") utils::globalVariables(".")
  # but that makes the tests fail.
  . = "ugly fix to please R CMD check"


  ########################
  # ASSERTS AND RECODING #
  ########################
  assert_logical(summary)
  assert_value(type, allowed = c("fitted", "predict", "residuals"))
  assert_types(probs, "logical", "numeric")
  if (is.numeric(probs))
    assert_numeric(probs, lower = 0, upper = 1)
  if (all(probs == TRUE))
    probs = c(0.025, 0.975)
  assert_logical(rate)
  assert_logical(prior)
  assert_logical(arma)
  assert_types(nsamples, "null", "numeric")
  if (!is.null(nsamples))
    assert_integer(nsamples, lower = 1)
  assert_value(samples_format, allowed = c("tidy", "matrix"))

  if (type == "residuals") {
    if (!is.null(newdata) && !(fit$pars$y %in% colnames(newdata)))
      stop("`newdata` must contain a response column named '", fit$pars$y, "' when `type = 'residuals'`")

    if (!(fit$pars$y %in% required_cols))
      required_cols = c(required_cols, fit$pars$y)
  }
  if (scale == "linear" && (type != "fitted"))
    stop("`scale = 'linear'` is only meaningful when `type = 'fitted'`.")


  # Check and build stuff related to newdata
  if (!is.null(newdata)) {
    assert_types(newdata, "data.frame", "tibble")

    # Add response column to required_cols for ARMA models
    if (is_arma == TRUE && arma == TRUE) {
      if (fit$pars$y %in% colnames(newdata) && !(fit$pars$y %in% required_cols)) {
        required_cols = c(required_cols, fit$pars$y)
      } else if (".ydata" %in% colnames(newdata)) {
        required_cols = c(required_cols, ".ydata")
      }
    }

    cols_in_newdata = required_cols %in% colnames(newdata)
    if (any(cols_in_newdata == FALSE))
      stop("Missing columns '", paste0(required_cols[!cols_in_newdata], collapse = "', '"), "' in `newdata`.")
  } else {
    # Use original data. Remember the response variable too for ARMA
    newdata = fit$data
    if (is_arma == TRUE && arma == TRUE && !(fit$pars$y %in% required_cols))
      required_cols = c(required_cols, fit$pars$y)
  }
  newdata = data.frame(newdata[, required_cols])
  colnames(newdata) = required_cols  # Special case for when there's only one predictor
  newdata$data_row = 1:nrow(newdata)  # to maintain order in the output when summary == TRUE


  ########################
  # GET FITS/PREDICTIONS #
  ########################
  type_for_simulate = ifelse(type == "residuals", yes = "fitted", no = type)
  if (length(varying_info$cols) > 0) {
    # If there are varying effects: use varying-matching samples for each row of data
    samples = dplyr::left_join(
      newdata,
      tidy_samples(fit, population = TRUE, varying = varying, prior = prior, nsamples = nsamples),
      by = unique(varying_info$cols),
      relationship = "many-to-many"
    ) %>%
      dplyr::mutate(!!returnvar := rlang::exec(fit$simulate, !!!., type = type_for_simulate, rate = rate, which_y = which_y, arma = arma, add_attr = FALSE, scale = scale))
  } else {
    # No varying effects: use all samples for each row of data
    samples = tidy_samples(fit, population = TRUE, varying = varying, prior = prior, nsamples = nsamples) %>%
      tidyr::expand_grid(newdata) %>%
      dplyr::mutate(!!returnvar := rlang::exec(fit$simulate, !!!., type = type_for_simulate, rate = rate, which_y = which_y, arma = arma, add_attr = FALSE, scale = scale))
  }

  # Optionally compute residuals
  if (type == "residuals")
    samples = dplyr::mutate(samples, !!returnvar := !!returnvar - !!yvar)  # returnvar should be "residuals" in this case

  # Optionally summarise
  if (summary == TRUE) {
    df_return = samples %>%
      # Summarise for each row in newdata
      dplyr::group_by(.data$data_row) %>%
      dplyr::summarise(.groups = "drop",
                       error = stats::sd(!!returnvar),
                       !!returnvar := mean(!!returnvar)
      ) %>%

      # Apply original order and put newdata as the first columns
      dplyr::arrange(.data$data_row) %>%
      dplyr::left_join(newdata, by = "data_row", relationship = "one-to-one") %>%
      dplyr::select(dplyr::one_of(colnames(newdata)), !!returnvar, "error", -"data_row")


    # Quantiles
    if (all(probs != FALSE)) {
      quantiles_fit = samples %>%
        get_quantiles(probs, xvar, returnvar) %>%
        dplyr::mutate(quantile = 100 * .data$quantile) %>%
        tidyr::pivot_wider(names_from = "quantile", names_prefix = "Q", values_from = "y")

      df_return = dplyr::left_join(df_return, quantiles_fit, by = as.character(xvar), relationship = "many-to-one")
    }
    return(data.frame(df_return))
  } else if (samples_format == "tidy") {
    return(samples)
  } else if (samples_format == "matrix") {
    df_return = tidy_to_matrix(samples, type)
    return(df_return)
  }
}


#' Samples from the Posterior Predictive Distribution
#'
#' @aliases predict predict.mcpfit
#' @inheritParams pp_eval
#' @inherit pp_eval return
#' @seealso \code{\link{pp_eval}} \code{\link{fitted.mcpfit}} \code{\link{residuals.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \donttest{
#' predict(demo_fit)  # Evaluate at each demo_fit$data
#' predict(demo_fit, probs = c(0.1, 0.5, 0.9))  # With median and 80% credible interval.
#' predict(demo_fit, summary = FALSE)  # Samples instead of summary.
#' predict(
#'   demo_fit,
#'   newdata = data.frame(time = c(-5, 20, 300)),  # Evaluate
#'   probs = c(0.025, 0.5, 0.975)
#' )
#'}
#'
predict.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  which_y = "ct",
  varying = TRUE,
  arma = TRUE,
  nsamples = NULL,
  samples_format = "tidy",
  ...
) {
  assert_ellipsis(...)
  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "predict",
    probs = probs,
    rate = rate,
    prior = prior,
    which_y = which_y,
    varying = varying,
    arma = arma,
    nsamples = nsamples,
    samples_format = samples_format
  )
}


#' Expected Values from the Posterior Predictive Distribution
#'
#' @aliases fitted fitted.mcpfit
#' @inheritParams pp_eval
#' @inherit pp_eval return
#' @seealso \code{\link{pp_eval}} \code{\link{predict.mcpfit}} \code{\link{residuals.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \donttest{
#' fitted(demo_fit)
#' fitted(demo_fit, probs = c(0.1, 0.5, 0.9))  # With median and 80% credible interval.
#' fitted(demo_fit, summary = FALSE)  # Samples instead of summary.
#' fitted(demo_fit,
#'        newdata = data.frame(time = c(-5, 20, 300)),  # New data
#'        probs = c(0.025, 0.5, 0.975))
#'}
#'
fitted.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  rate = TRUE,
  prior = FALSE,
  which_y = "ct",
  varying = TRUE,
  arma = TRUE,
  nsamples = NULL,
  samples_format = "tidy",
  scale = "response",
  ...
) {
  assert_ellipsis(...)
  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "fitted",
    probs = probs,
    rate = rate,
    prior = prior,
    which_y = which_y,
    varying = varying,
    arma = arma,
    nsamples = nsamples,
    samples_format = samples_format,
    scale = scale
  )
}


#' Compute Residuals From Mcpfit Objects
#'
#' Equivalent to  `fitted(fit, ...) - fit$data[, fit$data$yvar]` (or `fitted(fit, ...) - newdata[, fit$data$yvar]`),
#' but with fixed arguments for `fitted`: `rate = FALSE, which_y = 'ct', samples_format = 'tidy'`.
#'
#' @aliases residuals residuals.mcpfit resid resid.mcpfit
#' @inheritParams pp_eval
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @seealso \code{\link{pp_eval}} \code{\link{fitted.mcpfit}} \code{\link{predict.mcpfit}}
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \donttest{
#' residuals(demo_fit)
#' residuals(demo_fit, probs = c(0.1, 0.5, 0.9))  # With median and 80% credible interval.
#' residuals(demo_fit, summary = FALSE)  # Samples instead of summary.
#'}
#'
residuals.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  probs = TRUE,
  prior = FALSE,
  varying = TRUE,
  arma = TRUE,
  nsamples = NULL,
  ...
) {
  assert_ellipsis(...)
  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    type = "residuals",
    probs = probs,
    rate = FALSE,
    prior = prior,
    which_y = "ct",
    varying = varying,
    arma = arma,
    nsamples = nsamples,
    samples_format = "tidy"
  )
}


#' Convert from tidy to matrix
#'
#' Converts from the output of `tidy_samples()` or `pp_eval(fit, samples_format = "tidy")`
#' to an `N_draws` X `nrows(newdata)` matrix with fitted/predicted values. This format is
#' used y `brms` and it's useful as `yrep` in `bayesplot::ppc_*` functions.
#'
#' @aliases tidy_to_matrix
#' @keywords internal
#' @param samples Samples in tidy format
#' @param returnvar An `rlang::sym()` object.
#' @return An  `N_draws` X `nrows(newdata)` matrix.
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
tidy_to_matrix = function(samples, returnvar) {
  returnvar = rlang::sym(returnvar)
  samples %>%
    dplyr::select(".draw", "data_row", !!returnvar) %>%
    tidyr::pivot_wider(names_from = "data_row", values_from = !!returnvar) %>%
    dplyr::select(-".draw") %>%
    as.matrix()
}


#' Add loo if not already present
#'
#' @aliases with_loo
#' @keywords internal
#' @param fit An mcpfit object
#' @param save_psis Logical. See documentation of loo::loo
#' @param info Optional message if adding loo
#' @return An mcpfit object with loo.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
with_loo = function(fit, save_psis = FALSE, info = NULL) {
  assert_mcpfit(fit)

  # Add loo if absent or needs psis
  if (is.null(fit$loo) || (save_psis == TRUE && loo::is.loo(fit$loo) && is.null(fit$loo$psis_object))) {
    if (is.character(info))
      message(info)
    fit$loo = loo(fit, save_psis = save_psis)
  }

  return(fit)
}
