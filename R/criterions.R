#' Compute information criteria for model comparison
#'
#' Takes a \code{mcpfit} as input and computes information criteria using loo or
#' WAIC. Compare models using \code{\link[loo]{loo_compare}} and \code{\link[loo]{loo_model_weights}}.
#' more in \code{\link[loo]{loo}}.
#'
#' @aliases criterion
#' @param fit An mcpfit object.
#' @param criterion One of "loo" (calls loo::loo) or "waic" (calls loo::waic).
#' @return a \code{loo} or \code{psis_loo} object. See more in \code{\link[loo]{loo}}.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \dontrun{
#' # Compute loos
#' fit1$loo = loo(fit1)
#' fit2$loo = loo(fit2)
#' fit2$waic = waic(fit2)  # Just for fun
#' fit1$loo  # view it
#'
#' # Compare loos. Top is best. Should be several SDs better than others.
#' loo::loo_compare(fit1$loo, fit2$loo)
#'
#' # Compute model weights. Higher weight is better. See help for details..
#' loo::loo_model_weights(list(fit1$loo, fit2$loo))
#'}

criterion = function(fit, criterion = "loo") {
  if (!class(fit) == "mcpfit") {
    stop("class(fit) must be 'mcpfit'")
  }
  if (!criterion %in% c("loo", "waic")) {
    stop("criterion must be one of 'loo' or 'waic'")
  }

  # Log-likelihood MCMC samples as matrix
  loglik = as.matrix(do.call(rbind.data.frame, fit$loglik))

  # Add LOO
  if (any(criterion == "loo")) {
    # Compute relative effective sample size (for each loglik col)
    chain_id = rep(seq_len(length(fit$samples)), each = nrow(fit$samples[[1]]))
    r_eff = loo::relative_eff(exp(loglik), chain_id)  # Likelihood = exp(log-likelihood)

    # Add LOO
    return(loo::loo(loglik, r_eff = r_eff))
  }

  # Add WAIC
  if (any(criterion == "waic")) {
    return(loo::waic(loglik))
  }
}



#' @aliases loo LOO loo.mcpfit
#' @describeIn criterion Computes loo on mcpfit objects
#' @param x \code{mcpfit} object.
#' @param ... Currently ignored
#' @importFrom loo loo
#' @export loo
#' @export
loo.mcpfit = function(x, ...) {
  criterion(x, "loo")
}

#' @aliases waic WAIC waic.mcpfit
#' @describeIn criterion Computes WAIC on mcpfit objects
#' @param x \code{mcpfit} object.
#' @param ... Currently ignored
#' @importFrom loo waic
#' @export waic
#' @export
waic.mcpfit = function(x, ...) {
  criterion(x, "waic")
}
