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
  # Get posterior/prior samples
  samples = get_samples(fit, prior = prior)

  # Select only varying or only population-level columns in data
  if (varying == FALSE) {
    regex_pars = paste0(paste0("^", fit$pars$population, "$"), collapse = "|")
  } else {
    regex_pars = paste0("^", paste(fit$pars$varying, collapse="|^"))
    if (regex_pars == "^")
      stop("There were no matching parameters in the model.")
  }

  # HACK: If there is just one parameter, add two in to make the code run.
  all_cols = colnames(samples[[1]])
  get_cols = all_cols[stringr::str_detect(all_cols, regex_pars)]
  if (!stringr::str_detect(regex_pars, "\\|") & varying == FALSE) {
    samples = lapply(samples, function(x) x[, c(get_cols, "cp_0", "cp_1")])
  } else {
    samples = lapply(samples, function(x) x[, get_cols])
  }

  # Get parameter estimates
  # TO DO: Temporarily suppress warnings about tidyr1.0::unnest() until tidybayes 1.2 is out.
  suppressWarnings(
    estimates <- samples %>%  # whoa, a use case for <- over = !!!1!1one
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
  )

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
    #simulated = sapply(simulated, language_to_numeric)
    simulated = simulated[sapply(simulated, is.numeric)]  # Remove non-numeric

    # Handle varying effects. Finds the matching labels
    for (varying in fit$pars$varying) {
      if (!is.null(simulated[[varying]])) {
        # Find the needed values and labels
        value = simulated[[varying]]  # Extract simulation values
        label_col = stats::na.omit(fit$.other$ST$cp_group_col)[1]  # What column is the labels in data. TO DO: only works if there is just one grouping factor!
        labs = fit$data[[label_col]]  # Find the labels. Same length as `value`
        if (length(value) != length(labs)) {
          warning("This is simulated data, but the labels for varying effects in data does not have the same length as the numeric params used for simulation.")
          next
        }

        # Name like the MCMC columns and use only unique combinations (assuming identical value for each level)
        names(value) = paste0(varying, "[", labs, "]")
        value = value[!duplicated(value)]

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
  data.frame(result)
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

  if (width < 0 | width > 1)
    stop("`width`` must be between 0 and 1")

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
      cat("\nUse `ranef(object)` to summarise the varying effect(s):", paste0(fit$pars$varying, collapse = ", "))

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
