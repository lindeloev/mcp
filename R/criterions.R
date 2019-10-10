#' Compute information criteria for model comparison
#'
#' Takes a `mcpfit`` as input and computes information criteria using loo or
#' WAIC. Interpret directly or use `loo::loo_compare` to compare models.
#'
#' @param fit An mcmc.list with log-density column(s)
#' @param criterion One of "loo" (calls loo::loo) or "waic" (calls loo:waic).
#' @keywords information, loo, waic, mcmc
#' @import loo
#' @export
#' @examples
#' # Compute loos
#' fit1$loo = loo(fit1)
#' fit2$loo = loo(fit2)
#' fit2$waic = waic(fit2)  # Just for fun
#' fit1$loo  # view it
#'
#' # Compare loos. Top is best. Should be several SDs better than others.
#' loo::loo_compare(fit1$loo, fit2$loo)
#'
#' # Compute model weights. Higher weight is best. See help for details..
#' loo::loo_model_weights(list(fit1$loo, fit2$loo))
#'

criterion = function(fit, criterion = "loo") {
  if(!class(fit) == "mcpfit") {
    stop("class(fit) must be 'mcpfit'")
  }
  if(!criterion %in% c("loo", "waic")) {
    stop("criterion must be one of 'loo' or 'waic'")
  }

  # Log-likelihood MCMC samples as matrix
  loglik = as.matrix(do.call(rbind.data.frame, fit$loglik))

  # Add LOO
  if(any(criterion == "loo")) {
    # Compute relative effective sample size (for each loglik col)
    chain_id = rep(1:length(fit$samples), each = nrow(fit$samples[[1]]))
    r_eff = loo::relative_eff(exp(loglik), chain_id)  # Likelihood = exp(log-likelihood)

    # Add LOO
    return(loo::loo(loglik, r_eff = r_eff))
  }

  # Add WAIC
  if(any(criterion == 'waic')) {
    return(loo::waic(loglik))
  }
}



#' Computes loo on mcpfit objects
#' @export
#' @examples
#' fit$loo = loo(fit)
loo.mcpfit = function(fit) {
  criterion(fit, "loo")
}

#' Computes WAIC loo on mcpfit objects
#' @export
#' @examples
#' fit$waic = waic(fit)
waic.mcpfit = function(fit) {
  criterion(fit, "waic")
}
