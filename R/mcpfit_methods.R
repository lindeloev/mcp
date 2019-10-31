
#' Internal function for summary.mcpfit, fixef.mcpfit, and ranef.mcpfit
#'
#' @aliases get_summary
#' @inheritParams summary.mcpfit
#' @param fit An mcpfit object.
#' @param varying Boolean. Get results for varying (TRUE) or population (FALSE)?
#' @importFrom magrittr %>%
#' @importFrom dplyr .data
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_summary = function(fit, width, varying = FALSE) {
  # Select only varying or only population-level columns in data
  if (varying == FALSE) {
    regex_pars = paste0(paste0("^", fit$pars$population, "$"), collapse = "|")
  } else {
    regex_pars = paste0("^", paste(fit$pars$varying, collapse="|^"))
    if (regex_pars == "^")
      stop("There were no matching parameters in the model.")
  }


  all_cols = colnames(fit$samples[[1]])
  get_cols = all_cols[stringr::str_detect(all_cols, regex_pars)]
  samples = lapply(fit$samples, function(x) x[, get_cols])
  # HACK: If there is just one parameter, add two in to make the code run.
  if (!stringr::str_detect(regex_pars, "\\|") & varying == FALSE) {
    samples_org = samples
    samples = lapply(fit$samples, function(x) x[, c(get_cols, "cp_0", "cp_1")])
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
      dplyr::rename(mean = .data$value) %>%
      dplyr::rename(!!as.character(100*(1 - width) / 2) := .data$.lower,
                    !!as.character(100*(width + (1 - width)/2)) := .data$.upper) %>%

      # Remove unneeded stuf
      dplyr::select(-tidyselect::starts_with("."))
  )

  # Revert HACK and continue
  if (!stringr::str_detect(regex_pars, "\\|") & varying == FALSE) {
    estimates = dplyr::filter(estimates, !.data$name %in% c("cp_0", "cp_1"))
    samples = samples_org
  }

  # Diagnostics
  rhat = try(coda::gelman.diag(samples)$psrf[, 1], TRUE)
  if (!is.numeric(rhat)) {
    warning("rhat computation failed: ", rhat)
    rhat = rep(NA, nrow(estimates))
  }
  diagnostics = data.frame(
    rhat = rhat,  # Gelman-Rubin
    eff = round(coda::effectiveSize(samples)),  # Effective sample size
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
#' diagnostics. Get them in a data frame by doing result = summary(fit).
#'
#' @aliases summary summary.mcpfit
#' @param object An \code{mcpfit} object returned by \code{\link{mcp}}.
#' @param width Float. The width of the highest posterior density interval
#'   (between 0 and 1).
#' @param digits Positive integer. Number of digits to print
#' @param ... Currently ignored
#'
#' @return Posterior means and HDI intervals. \code{rhat} is the Gelman-Rubin
#'   convergence diagnostic which is often taken to be acceptable if < 1.1. It
#'   is computed using \code{\link[coda]{gelman.diag}}.
#'   \code{eff} is the effective sample size computed using
#'   \code{\link[coda]{effectiveSize}}. Low effective sample sizes are
#'   also obvious as poor mixing in trace plots (see \code{plot(fit, "combo")}).
#'   \code{ts_err} is the time-series error, taking autoregressive correlation
#'   into account. It is computed using \code{\link[coda]{spectrum0.ar}}.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \dontrun{
#' result = summary(fit, width = 0.8, digits = 4)
#' ranef(fit)  # varying (random) effects
#' }
summary.mcpfit = function(object, width = 0.95, digits = 2, ...) {
  # Standard name in mcp
  fit = object

  if (class(object) != "mcpfit")
    stop("`object`` must be an mcpfit object.")

  if (width < 0 | width > 1)
    stop("`width`` must be between 0 and 1")

  if (digits != floor(digits) | digits < 0)
    stop("`digits`` must be a positive integer.")


  # Model info
  cat("Family: ", fit$family$family, "(link = '", fit$family$link, "')\n", sep = "")
  if (!is.null(fit$samples))
    cat("Iterations: ", nrow(fit$samples[[1]]) * length(fit$samples), " from ", length(fit$samples), " chains.\n", sep="")
  cat("Model:\n")
  for (segment in fit$segments) {
    cat("  ", format(segment), "\n")
  }

  # Data
  if (!is.null(fit$samples)) {
    # Print and return invisibly
    cat("\nPopulation-level parameters:\n")
    result = get_summary(fit, width, varying = FALSE)
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
#' This is identical to what \code{summary(object)} does.
#'
#' @aliases fixef fixef.mcpfit fixed.effects
#' @inheritParams summary.mcpfit
#' @export
fixef = function(object, width = 0.95, ...) {
  get_summary(object, width, varying = FALSE)
}

#' Get varying ("random") effects of mcpfit.
#'
#' @aliases ranef ranef.mcpfit random.effects
#' @inheritParams summary.mcpfit
#' @export
ranef = function(object, width = 0.95, ...) {
  get_summary(object, width, varying = TRUE)
}


#' Print mcpfit
#'
#' Use \code{\link{summary.mcpfit}} for greater control.
#'
#' @aliases print print.mcpfit
#' @param x \code{mcpfit} object.
#' @param ... Currently ignored.
#' @export
print.mcpfit = function(x, ...) {
  summary(x)
}
