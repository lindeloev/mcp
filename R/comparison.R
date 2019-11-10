#' Compute information criteria for model comparison
#'
#' Takes a \code{mcpfit} as input and computes information criteria using loo or
#' WAIC. Compare models using \code{\link[loo]{loo_compare}} and \code{\link[loo]{loo_model_weights}}.
#' more in \code{\link[loo]{loo}}.
#'
#' @aliases criterion
#' @param fit An mcpfit object.
#' @param criterion One of "loo" (calls \code{\link[loo]{loo}}) or "waic" (calls \code{\link[loo]{waic}}).
#' @return a \code{loo} or \code{psis_loo} object.
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
  loglik = as.matrix(do.call(rbind.data.frame, fit$mcmc_loglik))

  # Add LOO
  if (criterion == "loo") {
    # Compute relative effective sample size (for each loglik col)
    chain_id = rep(seq_len(length(fit$mcmc_post)), each = nrow(fit$mcmc_post[[1]]))
    r_eff = loo::relative_eff(exp(loglik), chain_id)  # Likelihood = exp(log-likelihood)

    # Add LOO
    return(loo::loo(loglik, r_eff = r_eff))
  }

  # Add WAIC
  if (criterion == "waic") {
    return(loo::waic(loglik))
  }
}



#' @aliases loo LOO loo.mcpfit
#' @describeIn criterion Computes loo on mcpfit objects
#' @param x \code{mcpfit} object.
#' @param ... Currently ignored
#' @seealso \link{criterion}
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
#' @seealso \link{criterion}
#' @export waic
#' @export
waic.mcpfit = function(x, ...) {
  criterion(x, "waic")
}


#' Test hypotheses on mcp objects.
#'
#' This function is highly inspired by \code{\link[brms]{hypothesis}}. For
#' directional hypotheses, \code{hypothesis} executes the hypothesis string in
#' a \code{tidybayes} environment and summerises the proportion of samples where
#' the expression evaluates to TRUE. For equals-hypothesis, a Savage-Dickey
#' ratio is computed.
#'
#' @aliases hypothesis hypothesis.mcpfit
#' @param fit An mcpfit object
#' @param hypothesis String representation of a logical test involving model parameters.
#'   Takes R code that evaluates to TRUE or FALSE in a vectorized way.
#'
#'   Directional hypotheses are specified using <, >, <=, or >=. \code{hypothesis}
#'   returns the posterior probability and odds in favor of the stated hypothesis.
#'   The odds can be interpreted as a Bayes Factor. For example:
#'
#'   \itemize{
#'     \item \code{"cp_1 > 30"}:  the first change point is above 30.
#'     \item \code{"int_1 > int_2"}: the intercept is greater in segment 1 than 2.
#'     \item \code{"x_2 - x_1 <= 3"}: the difference between slope 1 and 2 is less than 3.
#'     \item \code{"int_1 > 20 & int_1 < 30"}: int_1 is between 20 and 30
#'     \item \code{"cp_1^2 < 30 | (log(x_1) + log(x_2)) > 5"}: many options.
#'     \item \code{"`cp_1_id[1]` > `cp_1_id[2]``}: id1 is greater than id2, as estimated
#'       through the varying-by-"id" change point in segment 1. Note that \code{``}
#'       required for varying effects.
#'   }
#'   Hypotheses can also test equality using the equal sign (=). This runs a
#'   Savage-Dickey test, i.e., the proportion by which the probability density
#'   has increased from the prior to the posterior at a given value. Therefore,
#'   it requires \code{mcp(sample = "both")}. There are two requirements:
#'   First, there can only be one equal sign, so don't use and (&) or or (|).
#'   Second, the point to test has to be on the right, and the variables on the left.
#'   .
#'   \itemize{
#'     \item \code{"cp_1 = 30"}: is the first change point at 30? Or to be more precise:
#'       by what factor has the credence in cp_1 = 30 risen/fallen when
#'       conditioning on the data, relative to the prior credence?
#'     \item \code{"int_1 + int_2 = 0"}: Is the sum of two intercepts zero?
#'     \item \code{"`cp_1_id[John]`/`cp_1_id[Erin]` = 2"}: is the varying change
#'       point for John (which is relative to \code{cp_1}) double that of Erin?
#'   }
#' @return A data.frame with a row per hypothesis
#' @export
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
hypothesis = function(fit, hypothesis) {
  ###########################
  # CHECK HYPOTHESIS FORMAT #
  ###########################
  expression = hypothesis
  expression = gsub(" ", "", expression)
  n_equals = stringi::stri_count(expression, regex="(?<!(<|>))=")
  n_directional = stringi::stri_count(expression, regex="<|<=|>|>=")

  if (n_equals > 1)
    stop("Only one equals-test (Savage-Dickey ratio) allowed in each hypothesis.")

  if (n_equals == 1 & n_directional > 0)
    stop("Equals cannot be combined with directional tests.")

  if (n_equals + n_directional == 0) {
    stop("At least one operator must be present: <, >, =, <=, or >=")
  }

  # If & or |, test that each side is legit
  sub_tests = unlist(strsplit(expression, "\\||&"))
  for (sub_test in sub_tests) {
    n_directional = stringi::stri_count(sub_test, regex="<|<=|>|>=|(?<!(<|>))=")
    if (n_directional != 1)
      stop("The hypothesis should contain exactly one directional operator (<, >, <=, or >=) at each side of & or |")
  }

  #################
  # SAVAGE-DICKEY #
  #################
  # Savage-Dickey
  if (n_equals == 1) {
    if (!coda::is.mcmc.list(fit$mcmc_prior) | !coda::is.mcmc.list(fit$mcmc_post))
      stop("Both prior and posterior samples are needed to compute Savage-Dickey density ratios. Run mcp(..., sample = 'both'")

    # Split into expression on the left and value on the right.
    expression_str = strsplit(expression, "=")[[1]][1]
    value_str = strsplit(expression, "=")[[1]][2]

    # Check them
    value = suppressWarnings(na.omit(as.numeric(value_str)))
    if (length(value) == 0)
      stop("Right-hand side has to be numeric. Was: ", value_str)

    # Finally, let's compute those densities
    dens_prior = get_density(fit$mcmc_prior, expression_str, value)
    dens_post = get_density(fit$mcmc_post, expression_str, value)
    odds = dens_post / dens_prior

    # If there is almost no density. somehow we get negative values.
    if (dens_post < 0 & dens_prior > 0)
      odds = 0
    if (dens_post > 0 & dens_prior < 0)
      odds = Inf

    return(data.frame(
      hypothesis = hypothesis,
      p = odds / (odds + 1),
      odds = odds
    ))


  ###############
  # DIRECTIONAL #
  ###############
  } else {
    # Get variable names
    # hyp_pars = strsplit(hyp, "<|>|=|\\(|\\)|\\+|\\-|(?<!_)[0-9.]", perl=T)
    # hyp_pars = fixit(hyp_pars)
    # hyp_pars = unique(hyp_pars)

    result = get_samples(fit) %>%
      #tidybayes::spread_draws(!!!dplyr::syms(hyp_pars)) %>%
      tidybayes::tidy_draws() %>%
      dplyr::mutate(result = eval(parse(text = expression))) %>%  # this is where the magic happens
      dplyr::summarise(
        prob = sum(result == TRUE) / dplyr::n()
      )

    return(data.frame(
      hypothesis = hypothesis,
      p = result$prob,
      odds = result$prob / (1 - result$prob)
    ))
  }
}


# fixit = function(stuff) {
#   x = unlist(stuff)
#   x = x[x != ""]  # remove empty
#   x = gsub(" ", "", x)
#   x
# }



#' Compute the density at a specific point.
#'
#' Used in \link{hypothesis}
#'
#' @aliases get_density
#' @param samples An mcmc.list
#' @param par_name Name of a column in the list
#' @param value What value to evaluate the density at
#' @return A float
#'
get_density = function(samples, eval_this, value) {
  samples = tidybayes::tidy_draws(samples) %>%
    dplyr::mutate(result = eval(parse(text = eval_this)))
  dens = density(dplyr::pull(samples, "result"))
  dens_point = spline(dens$x, dens$y, xout = value)$y
  return(dens_point)
}
