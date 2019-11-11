#' Compute information criteria for model comparison
#'
#' Takes a `mcpfit` as input and computes information criteria using loo or
#' WAIC. Compare models using \code{\link[loo]{loo_compare}} and \code{\link[loo]{loo_model_weights}}.
#' more in \code{\link[loo]{loo}}.
#'
#' @aliases criterion
#' @param fit An mcpfit object.
#' @param criterion One of `"loo"` (calls \code{\link[loo]{loo}}) or `"waic"` (calls \code{\link[loo]{waic}}).
#' @return a `loo` or `psis_loo` object.
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
    chain_id = rep(seq_along(fit$mcmc_post), each = nrow(fit$mcmc_post[[1]]))
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
#' @param x `mcpfit` object.
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
#' @param x `mcpfit` object.
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
#' This function is highly inspired by \code{\link[brms]{hypothesis}}. It returns
#' posterior probabilities and Bayes Factors for flexible hypotheses involving
#' model parameters. See the parameter `hypotheses`` for examples of how to
#' specify hypotheses.
#'
#' For directional hypotheses, `hypothesis`` executes the hypothesis string in
#' a `tidybayes`` environment and summerises the proportion of samples where
#' the expression evaluates to TRUE. For equals-hypothesis, a Savage-Dickey
#' ratio is computed. Savage-Dickey requires a prior too, so remember
#' `mcp(..., sample = "both")`.
#'
#' @aliases hypothesis hypothesis.mcpfit
#' @param fit An mcpfit object
#' @param hypotheses String representation of a logical test involving model parameters.
#'   Takes R code that evaluates to TRUE or FALSE in a vectorized way.
#'
#'   Directional hypotheses are specified using <, >, <=, or >=. `hypothesis`
#'   returns the posterior probability and odds in favor of the stated hypothesis.
#'   The odds can be interpreted as a Bayes Factor. For example:
#'
#'   * `"cp_1 > 30"`:  the first change point is above 30.
#'   * `"int_1 > int_2"`: the intercept is greater in segment 1 than 2.
#'   * `"x_2 - x_1 <= 3"`: the difference between slope 1 and 2 is less
#'       than or equal to 3.
#'   * `"int_1 > -2 & int_1 < 2"`: int_1 is between -2 and 2 (an interval hypothesis). This can be useful as a Region Of Practical Equivalence test (ROPE).
#'   * `"cp_1^2 < 30 | (log(x_1) + log(x_2)) > 5"`: be creative.
#'   * \code{"`cp_1_id[1]` > `cp_1_id[2]`"}: id1 is greater than id2, as estimated
#'       through the varying-by-"id" change point in segment 1. Note that \code{``}
#'       required for varying effects.
#'
#'   Hypotheses can also test equality using the equal sign (=). This runs a
#'   Savage-Dickey test, i.e., the proportion by which the probability density
#'   has increased from the prior to the posterior at a given value. Therefore,
#'   it requires `mcp(sample = "both")`. There are two requirements:
#'   First, there can only be one equal sign, so don't use and (&) or or (|).
#'   Second, the point to test has to be on the right, and the variables on the left.
#'
#'   * `"cp_1 = 30"`: is the first change point at 30? Or to be more precise:
#'       by what factor has the credence in cp_1 = 30 risen/fallen when
#'       conditioning on the data, relative to the prior credence?
#'   * `"int_1 + int_2 = 0"`: Is the sum of two intercepts zero?
#'   * ````"`cp_1_id[John]`/`cp_1_id[Erin]` = 2"````: is the varying change
#'       point for John (which is relative to `cp_1``) double that of Erin?
#' @param width Between 0 and 1. Width of the highest density interval.
#' @return A data.frame with a row per hypothesis and the following columns:
#'   * `hypothesis` is the hypothesis; often re-arranged to test against zero.
#'   * `mean` is the posterior mean of the left-hand side of the hypothesis.
#'   * `lower` is the lower bound of the (two-sided) highest-density interval of width `width`.
#'   * `upper` is the upper bound of ditto.
#'   * `p` Posterior probability.
#'       For "=" (Savage-Dickey), it is the BF converted to p.
#'       For directional hypotheses, it is the proportion of samples that returns TRUE.
#'   * `BF` Bayes Factor in favor  of the hypothesis.
#'       For "=" it is the Savage-Dickey density ratio.
#'       For directional hypotheses, it is p converted to odds.
#' @export
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
hypothesis = function(fit, hypotheses, width = 0.95) {
  # Loop through hypotheses and populate return_df
  return_df = data.frame()
  for (expression in hypotheses) {
    ####################
    # PREPARE FOR TEST #
    ####################
    # Check input
    n_equals = stringi::stri_count(expression, regex="(?<!(<|>))=")
    n_directional = stringi::stri_count(expression, regex="<|<=|>|>=")

    if (n_equals > 1)
      stop("Only one equals-test (Savage-Dickey ratio) allowed in each hypothesis: ", expression)

    if (n_equals == 1 & n_directional > 0)
      stop("Equals cannot be combined with directional tests: ", expression)

    if (n_equals + n_directional == 0) {
      stop("At least one operator must be present: <, >, =, <=, or >=: ", expression)
    }

    if (stringr::str_detect(expression, "\\[|\\]") & !stringr::str_detect(expression, "`"))
      stop("Needs `` around varying effects, e.g., `cp_1_id[2]`. Got this: ", expression)


    # If this is a single expression (does not contain & or |), we can estimate
    # the test value by putting everything on the LHS and zero on the RHS.
    if (!stringr::str_detect(expression, "\\||&")) {
      # Determine which comparator is used here
      comparators = c("=", "<", ">"," <=", ">=")
      which_comparator = stringr::str_detect(expression, comparators)
      this_comparator = comparators[which_comparator]

      # Re-arrange to LHS [comparator] 0.
      sides_split = strsplit(expression, "<|<=|>|>=|(?<!(<|>))=", perl = TRUE)[[1]]
      sides_split = stringr::str_trim(sides_split)
      if (stringr::str_detect(sides_split[2], "\\+|\\-"))
        sides_split[2] = paste0("(", sides_split[2], ")")
      LHS = paste0(sides_split[1], " - ", sides_split[2])
      expression = paste0(LHS, " ", this_comparator, " 0")

      # Get effect estimate
      samples = get_samples(fit) %>%
        tidybayes::tidy_draws() %>%
        dplyr::mutate(effect = eval(parse(text = LHS)))

      # TO DO: check need to suppress warnings when tidybayes 1.2 is out?
      estimate = suppressWarnings(tidybayes::mean_hdci(samples, .data$effect, .width = width))
    } else {
      samples = get_samples(fit) %>%
        tidybayes::tidy_draws()

      estimate = list(effect = NA, .lower = NA, .upper = NA)
    }

    # SAVAGE-DICKEY: compute p and BF
    if (n_equals == 1) {
      if (!coda::is.mcmc.list(fit$mcmc_prior) | !coda::is.mcmc.list(fit$mcmc_post))
        stop("Model contains '='. Both prior and posterior samples are needed to compute Savage-Dickey density ratios. Run mcp(..., sample = 'both'")

      # Finally, let's compute those densities
      dens_prior = get_density(fit$mcmc_prior, LHS, 0)
      dens_post = get_density(fit$mcmc_post, LHS, 0)
      BF = dens_post / dens_prior

      # If there is almost no density. somehow we get negative values.
      if (dens_post < 0 & dens_prior > 0)
        BF = 0
      if (dens_post > 0 & dens_prior < 0)
        BF = Inf

      p = BF / (BF + 1)
    }

    # DIRECTIONAL: compute p and BF
    if (n_directional != 0) {
      prob = samples %>%
        dplyr::mutate(result = eval(parse(text = expression))) %>%  # this is where the magic happens
        dplyr::summarise(
          prob = sum(result == TRUE) / dplyr::n()
        )

      p = prob$prob
      BF = prob$prob / (1 - prob$prob)
    }

    # Add to list
    new_row = data.frame(
      hypothesis = stringr::str_trim(expression),
      mean = estimate$effect,
      lower = estimate$.lower,
      upper = estimate$.upper,
      p = p,
      BF = BF,
      stringsAsFactors = FALSE
    )
    return_df = dplyr::bind_rows(return_df, new_row)
  }

  # Finally return
  return(return_df)
}


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
