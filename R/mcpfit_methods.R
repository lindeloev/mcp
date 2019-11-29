
#' Internal function for summary.mcpfit, fixef.mcpfit, and ranef.mcpfit
#'
#' @aliases get_summary
#' @inheritParams summary.mcpfit
#' @param fit An mcpfit object.
#' @param varying Boolean. Get results for varying (TRUE) or population (FALSE)?
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
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
    n.eff = round(coda::effectiveSize(samples)),  # Effective sample size
    ts_se = lapply(samples, coda::spectrum0.ar) %>%  # Time-series SE
      data.frame() %>%
      dplyr::select(tidyselect::starts_with("spec")) %>%
      rowMeans()
  )

  # Merge them and return
  result = dplyr::bind_cols(estimates, diagnostics)
  data.frame(result)
}




#' Summarise mcpfit objects
#'
#' Computes posterior means and highest density intervals, and model
#' diagnostics. Get them in a data frame by doing `result = summary(fit)`.
#'
#' @aliases summary summary.mcpfit
#' @param object An `mcpfit` object returned by \code{\link{mcp}}.
#' @param width Float. The width of the highest posterior density interval
#'   (between 0 and 1).
#' @param digits Positive integer. Number of digits to print
#' @param prior TRUE/FALSE. Summarise prior instead of posterior?
#' @param ... Currently ignored
#'
#' @return Posterior means and HDI intervals. `Rhat` is the Gelman-Rubin
#'   convergence diagnostic which is often taken to be acceptable if < 1.1. It
#'   is computed using \code{\link[coda]{gelman.diag}}.
#'   `n.eff` is the effective sample size computed using
#'   \code{\link[coda]{effectiveSize}}. Low effective sample sizes are
#'   also obvious as poor mixing in trace plots (see `plot(fit, "combo")`).
#'   `ts_err` is the time-series error, taking autoregressive correlation
#'   into account. It is computed using \code{\link[coda]{spectrum0.ar}}.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \dontrun{
#' result = summary(fit, width = 0.8, digits = 4)
#' ranef(fit)  # varying (random) effects
#' }
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
  for (i in 1:length(fit$segments)) {
    cat("  ", i, ": ", format(fit$segments[i]), "\n", sep = "")
  }

  # Data
  if (!is.null(samples)) {
    # Print and return invisibly
    cat("\nPopulation-level parameters:\n")
    result = get_summary(fit, width, varying = FALSE, prior = prior)
    print(data.frame(result), digits = digits + 1, row.names = FALSE)

    if (!is.null(fit$pars$varying))
      cat("\nUse `ranef(object)` to summarise the varying effect(s):", paste0(fit$pars$varying, collapse = ", "))

    return(invisible(result))
  }
  else {
    cat("\nNo samples. Nothing to summarise.")
    return(invisible(NULL))
  }
}



#' Get population-level ("fixed") effects of mcpfit.
#'
#' This is identical to what `summary(object)` does.
#'
#' @aliases fixef fixef.mcpfit fixed.effects
#' @inheritParams summary.mcpfit
#' @export
fixef = function(object, width = 0.95, prior = FALSE, ...) {
  get_summary(object, width, varying = FALSE, prior = prior)
}

#' Get varying ("random") effects of mcpfit.
#'
#' @aliases ranef ranef.mcpfit random.effects
#' @inheritParams summary.mcpfit
#' @export
ranef = function(object, width = 0.95, prior = FALSE, ...) {
  get_summary(object, width, varying = TRUE, prior = prior)
}


#' Print mcpfit
#'
#' Use \code{\link{summary.mcpfit}} for greater control.
#'
#' @aliases print print.mcpfit
#' @param x `mcpfit` object.
#' @param ... Currently ignored.
#' @export
print.mcpfit = function(x, ...) {
  summary(x)
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
#' @inheritParams summary.mcpfit
#' @param fit An mcpfit object
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
