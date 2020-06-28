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
#' @aliases get_summary
#' @keywords internal
#' @inheritParams summary.mcpfit
#' @param fit An \code{\link{mcpfit}}` object.
#' @param varying Boolean. Get results for varying (TRUE) or population (FALSE)?
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_summary = function(fit, width, varying = FALSE, prior = FALSE) {
  # Check arguments
  if (class(fit) != "mcpfit")
    stop("`object`` must be an mcpfit object.")

  if (!is.numeric(width) | width < 0 | width > 1 | length(width) != 1)
    stop("`width`` must be a float between 0 and 1. Got: ", width)

  if (!is.logical(varying))
    stop("`varying` must be TRUE or FALSE. Got: ", varying)

  if (!is.logical(prior))
    stop("`prior` must be TRUE or FALSE. Got: ", prior)


  # Get posterior/prior samples
  # TO DO: transition this to tidy_samples at some point if rhat and n.eff
  #    can be computed on this data. Would avoid some hacky solutions in the code below.
  samples = get_samples(fit, prior = prior)

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
  if (!stringr::str_detect(regex_pars, "\\|") & varying == FALSE) {
    samples = lapply(samples, function(x) x[, c(get_cols, "cp_0", "cp_1")])
  } else {
    samples = lapply(samples, function(x) x[, get_cols])
  }

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
  if (!stringr::str_detect(regex_pars, "\\|") & varying == FALSE) {
    estimates = dplyr::filter(estimates, !.data$name %in% c("cp_0", "cp_1"))
    samples = lapply(samples, function(x) x[, get_cols])
  }

  # Diagnostics
  Rhat = try(coda::gelman.diag(samples, multivariate = FALSE)$psrf[, 1], TRUE)
  if (!is.numeric(Rhat)) {
    warning("Rhat computation failed: ", Rhat)
    Rhat = rep(NA, nrow(estimates))
  }
  diagnostics = data.frame(
    Rhat = Rhat,  # Gelman-Rubin
    n.eff = round(coda::effectiveSize(samples))  # Effective sample size
  )

  # Add simulation parameters if the data is simulated
  sim_list = attr(fit$data[, fit$pars$y], "simulated")
  if(!is.null(sim_list)) {
    simulated = as.list(sim_list)  # Get as oroper list
    simulated = simulated[sapply(simulated, is.numeric)]  # Remove non-numeric

    # Handle varying effects. Finds the matching labels
    for (varying in fit$pars$varying) {
      if (!is.null(simulated[[varying]])) {
        # Find the needed values and labels
        value = simulated[[varying]]  # Extract simulation values
        label_col = stats::na.omit(fit$.other$ST$cp_group_col)[1]  # What column is the labels in data. TO DO: only works if there is just one grouping factor!
        labs = fit$data[[label_col]]  # Find the labels. Same length as `value`
        if (length(value) != length(labs)) {
          warning("This is simulated data, but the labels for varying effect '", label_col, "' in data does not have the same length as the numeric params used for simulation.")
          next
        }

        # Name like the MCMC columns and use only unique combinations (assuming identical value for each level)
        value = value[!duplicated(value)]
        names(value) = unique(paste0(varying, "[", labs, "]"))

        # Delete the simulation vector and add the new label-value pairs to list
        simulated[[varying]] = NULL
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
      dplyr::left_join(simulated, by = "name") %>%
      dplyr::mutate(match = ifelse(.data$sim > .data$lower & .data$sim < .data$upper, yes = "OK", no = "")) %>%
      dplyr::select(.data$name, .data$match, .data$sim, .data$mean, .data$lower, .data$upper)
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
#'     (see `plot_pars(fit)`). Read how to deal with such problems [here](https://lindeloev.github.io/mcp/articles/debug.html)
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
#' summary(ex_fit)
#' summary(ex_fit, width = 0.8, digits = 4)  # Set HDI width
#'
#' # Get the results as a data frame
#' results = summary(ex_fit)
#'
#' # Varying (random) effects
#' # ranef(my_fit)
#'
#' # Summarise prior
#' summary(ex_fit, prior = TRUE)
#'
summary.mcpfit = function(object, width = 0.95, digits = 2, prior = FALSE, ...) {
  # Standard name in mcp
  fit = object
  samples = get_samples(fit, prior = prior, error = FALSE)

  if (class(object) != "mcpfit")
    stop("`object`` must be an mcpfit object.")

  if (digits != floor(digits) | digits < 0)
    stop("`digits`` must be a positive integer.")


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
  get_summary(object, width, varying = FALSE, prior = prior)
}


#' @aliases ranef ranef.mcpfit random.effects
#' @describeIn summary.mcpfit Get varying ("random") effects of an \code{\link{mcpfit}} object.
#' @export
ranef = function(object, width = 0.95, prior = FALSE, ...) {
  get_summary(object, width, varying = TRUE, prior = prior)
}


#' @aliases print print.mcpfit
#' @describeIn summary.mcpfit Print the posterior summary of an \code{\link{mcpfit}} object.
#' @param x An \code{\link{mcpfit}} object.
#' @export
print.mcpfit = function(x, ...) {
  summary(x, ...)
}

#' Print mcpprior
#'
#' The mcpprior is just a list, but it can be displayed in a more condensed
#' way using cbind.
#' @aliases print.mcpprior
#' @inheritParams print.mcpfit
print.mcpprior = function(x, ...) {
  to_print = cbind(x)
  colnames(to_print) = "prior"
  print(to_print)
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
#' @aliases get_samples
#' @keywords internal
#' @inheritParams summary.mcpfit
#' @param fit An \code{\link{mcpfit}} object
#' @param message TRUE: gives a message if returning prior samples. FALSE = no message
#' @param error TRUE: err if there are no samples. FALSE: return NULL
get_samples = function(fit, prior = FALSE, message = TRUE, error = TRUE) {
  if (!is.logical(prior))
    stop("`prior` must be TRUE or FALSE.")

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


#' Get tidy samples with or without varying effects
#'
#' Returns in a format useful for `fit$simulate()` with population parameters in wide format
#' and varying effects in long format (the number of rows will be `nsamples * n_levels_in_varying`).
#'
#' @aliases get_tidy_samples
#' @keywords internal
#' @inheritParams get_samples
#' @param population
#'   * \strong{TRUE:} All population effects. Same as `fit$pars$population`.
#'   * \strong{FALSE:} No population effects. Same as `c()`.
#'   * \strong{Character vector:} Only include specified population parameters - see `fit$pars$population`.
#' @param varying
#'   * \strong{TRUE:} All varying effects. Same as `fit$pars$varying`.
#'   * \strong{FALSE:} No varying efects. Same as `c()`.
#'   * \strong{Character vector:} Only include specified varying parameters - see `fit$pars$varying`.
#' @param absolute
#'   * \strong{TRUE:} Returns the absolute location of all varying change points.
#'   * \strong{FALSE:} Just returns the varying effects.
#'   * \strong{Character vector:} Only do absolute transform for these varying parameters - see `fit$pars$varying`.
#'
#'   OBS: This currently only applies to varying change points. It is not implemented for `rel()` regressors yet.
#' @return `tibble` of posterior draws in `tidybayes` format.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @example
tidy_samples = function(fit, population = TRUE, varying = TRUE, absolute = FALSE, prior = FALSE) {
  # General argument checks
  if (all(population == FALSE) & all(varying == FALSE))
    stop("At least one TRUE or one parameter must be provided through either the `varying` or the `population` arguments.")

  # ----- IDENTIFY PARAMETERS -----
  # Varying parameters. Result is `terms_varying`
  if (all(varying == FALSE)) {
    # Select no varying effects
    use_varying = c()
  } else if (all(varying == TRUE)) {
    # Select all varying effects
    use_varying = !is.na(fit$.other$ST$cp_group)
  } else if (is.character(varying)) {
    if (!all(varying %in% fit$pars$varying))
      stop("Not all `varying` are varying parameters (see fit$pars$varying).")
    # Select only specified varying effects
    use_varying = fit$.other$ST$cp_group %in% varying
  } else {
    stop("`varying` must be TRUE, FALSE, or a character vector.")
  }
  pars_varying = fit$.other$ST$cp_group[use_varying]
  cols_varying = fit$.other$ST$cp_group_col[use_varying]
  terms_varying = paste0(pars_varying, "[", cols_varying, "]")  # for tidybayes
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
  } else {
    stop("`population` must be TRUE, FALSE, or a character vector.")
  }

  # Absolute effects. Results is `absolute_cps` and `absolute` (recoded to varying cp names)
  if (all(absolute == TRUE)) {
    absolute = fit$.other$ST$cp_group[use_varying]
    absolute_cps = fit$.other$ST$cp_name[use_varying]
  } else if (all(absolute == FALSE)) {
    absolute_cps = c()
  } else if (is.character(absolute)) {
    # Check
    is_in_varying = absolute %in% pars_varying
    if (any(!is_in_varying))
        stop("The following parameter names in `absolute` are not in `varying`:", absolute[is_in_varying])

    use_absolute = fit$.other$ST$cp_group %in% absolute
    absolute_cps = fit$.other$ST$cp_name[use_absolute]
  } else {
    stop("`absolute` must be TRUE, FALSE, or a character vector.")
  }

  # ----- GET THESE PARAMETERS AS TIDY DRAWS -----
  # Select posterior/prior samples
  samples = get_samples(fit, prior = prior)

  # Build code for tidybayes::spread_draws() and execute it
  all_terms = unique(c(pars_population, terms_varying, absolute_cps))
  code = paste0("tidybayes::spread_draws(samples, ", paste0(all_terms, collapse = ", "), ")")
  tidy_samples = eval(parse(text = code))

  # Add population cp to varying and delete population cols only included for this purpose.
  if (length(absolute_cps) > 0) {
    tidy_samples[, absolute] = tidy_samples[, absolute] + tidy_samples[, absolute_cps]
    tidy_samples = dplyr::select(tidy_samples, -!!dplyr::setdiff(absolute_cps, pars_population))
  }

  # Return with chain etc. first
  tidy_samples = dplyr::select(tidy_samples, .chain, .iteration, .draw, everything())
  return(tidy_samples)
}




#' Predictions from samples
#'
#' @aliases pp_eval pp_eval.mcpfit
#' @keywords internal
#' @inheritParams plot.mcpfit
#' @param newdata A `tibble` or a `data.frame` containing predictors in the model. If `NULL` (default), the original data is used.
#' @param summary Summarise at each x-value
#' @param q_fit Vector of quantiles. Only in effect when `summary == TRUE`.
#' @param q_predict Vector of quantiles. Only in effect when `summary == TRUE`.
#' @return
#'   * If `summary = TRUE`: A `tibble` with the posterior mean for each row in `newdata`, If `newdata` is `NULL`, the data in `fit$data` is used.
#'   * If `summary = FALSE`: A `tidybayes` `tibble` with all the posterior samples (`Ns`) evaluated at each row in `newdata` (`Nn`), i.e., with `Ns x Nn` rows.
#' @seealso fitted.mcpfit predict.mcpfit
#'
pp_eval = function(
  object,
  newdata = NULL,
  summary = TRUE,
  q_fit = FALSE,
  q_predict = FALSE,
  rate = TRUE,
  prior = FALSE,
  which_y = "ct",
  varying = TRUE
) {
  # Naming convention in mcp
  fit = object

  # Check inputs
  if (length(fit$pars$arma) != 0)
    stop("`predict.mcpfit()` is not implemented for ARMA models yet. Raise an issue on GitHub if you need it.")

  if (!is.null(newdata)) {
    if (!is.data.frame(newdata) & !tibble::is.tibble(newdata))
      stop("`newdata` must be a data.frame or a tibble. Got '", paste0(class(newdata), collapse = "', '"), "'")

    required_cols = c(fit$pars$x, fit$pars$trials, fit$pars$varying)
    cols_in_newdata = required_cols %in% colnames(newdata)
    if (any(cols_in_newdata == FALSE)) {
      stop("Missing columns '", paste0(required_cols[!cols_in_newdata], collapse = "', '"), "' in `newdata`.")
    }
  } else {
    newdata = fit$data
  }

  if (all(q_fit == TRUE))
    q_fit = c(0.025, 0.975)

  if (all(q_predict == TRUE))
    q_predict = c(0.025, 0.975)

  if (!is.logical(q_fit) & !is.numeric(q_fit))
    stop("`q_fit` has to be TRUE, FALSE, or a vector of numbers.")

  if (is.numeric(q_fit) & (any(q_fit > 1) | any(q_fit < 0)))
    stop ("All `q_fit` have to be between 0 (0%) and 1 (100%).")

  if (!is.logical(q_predict) & !is.numeric(q_predict))
    stop("`q_predict` has to be TRUE, FALSE, or a vector of numbers.")

  if (is.numeric(q_predict) & (any(q_predict > 1) | any(q_predict < 0)))
    stop ("All `q_predict` have to be between 0 (0%) and 1 (100%).")

  if (!is.logical(rate))
    stop("`rate` has to be TRUE or FALSE.")


  samples = get_tidy_samples(fit, prior, varying)

  # Summarise
  xvar = rlang::sym(fit$pars$x)
  yvar = rlang::sym(fit$pars$y)
  if (summary == TRUE) {
    df_return = suppressMessages(samples %>%  # TO DO: For dplyr 1.0.0 summarise... until it gets more widely adopted and .groups is not experimental anymore.
                                   dplyr::group_by(!!xvar) %>%
                                   dplyr::summarise(
                                     mean = mean(!!yvar),
                                     Std
                                   ))

    # Quantiles of estimated value
    if (all(q_fit != FALSE)) {
      quantiles_fit = samples %>%
        get_quantiles(q_fit, xvar, yvar, NULL) %>%
        dplyr::mutate(quantile = 100*quantile) %>%
        tidyr::pivot_wider(names_from = quantile, names_prefix = "fit", values_from = y)

      df_return = dplyr::left_join(df_return, quantiles_fit, by = as.character(xvar))
    }

    # Quantiles of prediction
    if(all(q_predict != FALSE)) {
      quantiles_predict = samples %>%
        dplyr::mutate(predicted_ := rlang::exec(simulate, !!!., type = "predict", rate = rate, which_y = which_y, add_attr = FALSE)) %>%
        get_quantiles(q_predict, xvar, rlang::sym("predicted_"), NULL) %>%
        dplyr::mutate(quantile = 100*quantile) %>%
        tidyr::pivot_wider(names_from = quantile, names_prefix = "predict", values_from = y)

      df_return = dplyr::left_join(df_return, quantiles_predict, by = as.character(xvar))
    }

    return(df_return)
  } else {
    return(samples)
  }
}


#' Samples from the Posterior Predictive Distribution
#'
#' @aliases predict predict.mcpfit
#' @inheritParams pp_eval
#' @inherit pp_eval return
#' @seealso pp_eval fitted.mcpfit
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' predict(ex_fit)  # Evaluate at each ex_fit$data
#' predict(ex_fit, q_predict = TRUE)  # Same, but add prediction intervals
#' fitted(ex_fit, summary = FALSE)  # Samples instead of summary.
#' predict(
#'   ex_fit,
#'   newdata = data.frame(x = c(-5, 20, 300)),  # Evaluate
#'   q_predict = c(0.025, 0.5, 0.975)
#' )
#'
predict.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  q_predict = TRUE,
  rate = TRUE,
  prior = FALSE,
  which_y = "ct"
) {
  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    q_predict = q_predict,
    rate = rate,
    prior = prior,
    which_y = which_y
  )
}


#' Expected Values from the Posterior Predictive Distribution
#'
#' @aliases fitted fitted.mcpfit
#' @inheritParams pp_eval
#' @inherit pp_eval return
#' @seealso pp_eval predict.mcpfit
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' fitted(ex_fit)
#' fitted(ex_fit, q_fit = TRUE)  # With 95% credible interval.
#' fitted(ex_fit, summary = FALSE)  # Samples instead of summary.
#' fitted(ex_fit,
#'        newdata = data.frame(x = c(-5, 20, 300)),  # New data
#'        q_predict = c(0.025, 0.5, 0.975))
#'
fitted.mcpfit = function(
  object,
  newdata = NULL,
  summary = TRUE,
  q_fit = TRUE,
  rate = TRUE,
  prior = FALSE,
  which_y = "ct"
) {
  pp_eval(
    object,
    newdata = newdata,
    summary = summary,
    q_fit = q_fit,
    rate = rate,
    prior = prior,
    which_y = which_y
  )
}
