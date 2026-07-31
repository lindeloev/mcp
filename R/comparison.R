# ABOUT: These functions compare models and/or hypotheses.
# -----------------

#' Information Criteria for Model Comparison
#'
#' Compare models using \code{\link[loo]{loo_compare}} and \code{\link[loo]{loo_model_weights}}.
#' more in \code{\link[loo]{loo}}.
#'
#' @aliases loo LOO loo.mcpfit
#' @inheritParams pp_eval
#' @param x An \code{\link{mcpfit}} object.
#' @param ... Further arguments passed to \code{\link[loo]{loo}}, e.g., `cores` or `save_psis`.
#' @param pointwise `TRUE` calls calls \code{\link[loo]{loo.function}} which is slower but more memory efficient.
#'   `FALSE` calls the default \code{\link[loo]{loo}}.
#' @param ndraws Integer or `NULL`. Number of posterior draws used for the
#'   log-likelihood or information criterion. `NULL` uses all draws.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @details Observationwise PSIS-LOO and WAIC are problematic for AR/MA models
#'   because both treat individual conditional likelihood terms as validation
#'   units. In PSIS-LOO, a held-out response also remains in the conditioning
#'   history of later terms. Prefer leave-future-out or blocked
#'   cross-validation, which are not currently implemented in mcp.
#' @return a `loo` or `psis_loo` object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export loo
#' @export
#' @examples
#' \donttest{
#' # Define two models and sample them
#' # future::plan(future::multisession, workers = 3)  # Uncomment for parallel sampling
#' data = mcp_example_data("intercepts")  # Get some simulated data.
#' model1 = list(y ~ 1 + x, ~ 1)
#' model2 = list(y ~ 1 + x)  # Without a change point
#' fit1 = mcp(model1, data)
#' fit2 = mcp(model2, data)
#'
#' # Compute LOO for each and compare (works for waic(fit) too)
#' fit1$loo = loo(fit1)
#' fit2$loo = loo(fit2)
#' loo::loo_compare(fit1$loo, fit2$loo)
#' }
loo.mcpfit = function(x, ..., pointwise = FALSE, varying = TRUE, arma = TRUE,
                      ndraws = NULL, nsamples = lifecycle::deprecated()) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "loo.mcpfit")
  fit = x
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_multi_class(varying, c("logical", "character"))
  checkmate::assert_flag(arma)
  warn_arma_check(fit, arma, "information_criterion")
  ndraws = validate_loglik_ndraws(fit, ndraws)
  mcmc_post = .subset2(fit, "mcmc_post")
  n_draws = sum(vapply(mcmc_post, nrow, integer(1)))
  settings = get_loglik_settings(fit, varying, arma, ndraws)
  if (length(settings$observed_rows) == 0)
    stop("LOO requires at least one observed response.")
  chain_id_all = rep(
    seq_along(mcmc_post),
    vapply(mcmc_post, nrow, integer(1))
  )


  # Matrix: Fast but memory-greedy matrix-based computation
  if (pointwise == FALSE) {
    if (!loglik_settings_match(fit$loglik, settings))
      fit = add_loglik(fit, varying = varying, arma = arma, ndraws = ndraws)
    chain_id = as.integer(rownames(fit$loglik))
    r_eff = if (is.null(ndraws)) {
      loo::relative_eff(exp(fit$loglik), chain_id)
    } else {
      1
    }
    result = loo::loo(fit$loglik, r_eff = r_eff, ...)

  # Pointwise: per-data-row computation
  } else {
    if (is.null(ndraws)) {
      n_eval = n_draws
      chain_id = chain_id_all
    } else {
      # Select draws once and reuse them for every observation. Treat the
      # arbitrary subset as one chain; relative_eff is set to 1 below.
      draw_indices = sort(sample.int(n_draws, ndraws))
      selected = as.matrix(.subset2(fit, "mcmc_post"))[draw_indices, , drop = FALSE]
      fit$mcmc_post = coda::mcmc.list(coda::mcmc(selected))
      n_eval = ndraws
      chain_id = rep(1L, ndraws)
    }
    predictors = get_fit_model_tables(fit)$predictors
    ar_order = get_arma_order(predictors, "ar")
    ma_order = get_arma_order(predictors, "ma")

    # For small models, the majority of the computation time will be pp_eval overhead
    llfun = function(data_i, draws = NULL, link_fun = identity) {
      original_row = data_i$.mcp_original_row
      if (is.na(ar_order) && is.na(ma_order)) {
        loglik_samples = pp_eval(
          fit, newdata = fit$data[original_row, , drop = FALSE],
          summary = FALSE, type = "loglik",
          varying = varying, arma = arma
        )
      } else {
        # AR needs only its direct response lags. MA innovations recursively
        # depend on all preceding observations, so retain the full history.
        first_row = if (is.na(ma_order)) max(1, original_row - ar_order) else 1
        data_rows = seq(first_row, original_row)
        lldata = fit$data[data_rows, ]
        loglik_samples = fit %>%
          pp_eval(newdata = lldata, summary = FALSE, type = "loglik", varying = varying, arma = arma) %>%
          dplyr::filter(.data$data_row == max(.data$data_row))  # last observed row
      }

      # Matrix conversion validates and explicitly orders the joint draws, so
      # every observation passed to loo uses the same draw identity.
      loglik = as.numeric(tidy_to_matrix(loglik_samples, type = ".loglik")[, 1])
      link_fun(loglik)
    }

    loo_data = data.frame(.mcp_original_row = settings$observed_rows)

    r_eff = if (is.null(ndraws)) {
      loo::relative_eff(
        llfun, data = loo_data, chain_id = chain_id,
        link_fun = exp, draws = seq_len(n_eval)
      )
    } else {
      1
    }
    result = loo::loo.function(
      llfun, data = loo_data, r_eff = r_eff,
      draws = seq_len(n_eval), ...
    )
  }

  attr(result, "mcp_settings") = settings
  result
}


#' @aliases waic WAIC waic.mcpfit
#' @describeIn loo.mcpfit Computes WAIC on mcpfit objects
#' @inheritParams loo.mcpfit
#' @param ... Currently ignored
#' @export waic
#' @export
waic.mcpfit = function(x, ..., varying = TRUE, arma = TRUE, ndraws = NULL,
                       nsamples = lifecycle::deprecated()) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "waic.mcpfit")
  rlang::check_dots_empty()
  fit = x
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_multi_class(varying, c("logical", "character"))
  checkmate::assert_flag(arma)
  warn_arma_check(fit, arma, "information_criterion")
  ndraws = validate_loglik_ndraws(fit, ndraws)
  settings = get_loglik_settings(fit, varying, arma, ndraws)
  if (!loglik_settings_match(fit$loglik, settings))
    fit = add_loglik(fit, varying = varying, arma = arma, ndraws = ndraws)

  loo::waic(fit$loglik)
}


#' Add Log-Likelihood to an mcpfit Object.
#'
#' @aliases add_loglik
#' @inheritParams loo.mcpfit
#' @seealso loo.mcpfit waic.mcpfit
#' @return An `mcpfit` object with `fit$loglik` filled as an (Nchains * Ndraws)
#'   by N-observed-responses data matrix with chain number as rownames. Rows
#'   with missing responses are excluded from the likelihood columns.
#' @export
#' @examples
#' \donttest{
#' demo_fit = add_loglik(demo_fit)
#' }
add_loglik = function(x, varying = TRUE, arma = TRUE, ndraws = NULL,
                      nsamples = lifecycle::deprecated()) {
  ndraws = resolve_ndraws(ndraws, nsamples, missing(ndraws), "add_loglik")
  fit = x
  ndraws = validate_loglik_ndraws(fit, ndraws)
  settings = get_loglik_settings(fit, varying, arma, ndraws)
  loglik_samples = pp_eval(
    fit, type = "loglik", summary = FALSE, probs = FALSE,
    varying = varying, arma = arma, ndraws = ndraws
  )
  draw_index = loglik_samples %>%
    dplyr::select(".draw", ".chain") %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$.draw)
  if (anyDuplicated(draw_index$.draw))
    stop_github("Chain metadata differs across evaluation rows for the same `.draw`.")
  fit$loglik = tidy_to_matrix(loglik_samples, type = ".loglik", data_rows = settings$observed_rows)

  # Chain info
  rownames(fit$loglik) = draw_index$.chain
  attr(fit$loglik, "mcp_settings") = settings

  fit
}


#' Validate the number of draws used for information criteria
#'
#' @keywords internal
#' @noRd
validate_loglik_ndraws = function(fit, ndraws) {
  checkmate::assert_int(ndraws, lower = 1, null.ok = TRUE)
  if (is.null(ndraws))
    return(NULL)

  mcmc_post = .subset2(fit, "mcmc_post")
  n_draws = sum(vapply(mcmc_post, nrow, integer(1)))
  if (ndraws > n_draws)
    stop("`ndraws` cannot exceed the ", n_draws, " available posterior draws.")
  as.integer(ndraws)
}


#' Record the settings represented by a log-likelihood or LOO object
#'
#' @keywords internal
#' @noRd
get_loglik_settings = function(fit, varying, arma, ndraws) {
  if (!is.null(ndraws))
    ndraws = as.integer(ndraws)
  list(
    varying = varying,
    arma = arma,
    ndraws = ndraws,
    observed_rows = which(!is.na(fit$data[, fit$pars$y]))
  )
}


#' Check whether cached log likelihoods match requested settings
#'
#' @keywords internal
#' @noRd
loglik_settings_match = function(loglik, settings) {
  !is.null(loglik) && identical(attr(loglik, "mcp_settings"), settings)
}


#' Test Hypotheses Concerning Individual Parameters
#'
#' Returns posterior probabilities and Bayes Factors for flexible hypotheses involving
#' model parameters. The documentation for the argument `hypotheses` below
#' shows examples of how to specify hypotheses, and [read worked examples on the mcp website](https://lindeloev.github.io/mcp/articles/comparison.html).
#' For directional hypotheses, `hypothesis` executes the hypothesis string in
#' a data-frame environment and summarises the proportion of posterior and
#' prior samples where the expression evaluates to TRUE. The Bayes factor is
#' the posterior odds divided by the prior odds. For equals-hypotheses, a
#' Savage-Dickey ratio is computed. Both kinds of Bayes factor require prior
#' samples, so remember `mcp(..., sample = "both")`. This function is heavily inspired by the
#' `hypothesis` function from the `brms` package.
#'
#' @aliases hypothesis hypothesis.mcpfit
#' @inheritParams summary.mcpfit
#' @param fit An \code{\link{mcpfit}} object.
#' @param hypotheses String representation of a logical test involving model parameters.
#'   Takes R code that evaluates to TRUE or FALSE in a vectorized way.
#'
#'   Directional hypotheses are specified using <, >, <=, or >=. `hypothesis`
#'   returns the posterior probability and the Bayes factor in favor of the
#'   stated hypothesis. The Bayes factor requires both prior and posterior
#'   samples from `mcp(sample = "both")`. For example:
#'
#'   * `"cp_1 > 30"`:  the first change point is above 30.
#'   * `"Intercept_1 > Intercept_2"`: the intercept is greater in segment 1 than 2.
#'   * `"x_2 - x_1 <= 3"`: the difference between slope 1 and 2 is less
#'       than or equal to 3.
#'   * `"Intercept_1 > -2 & Intercept_1 < 2"`: Intercept_1 is between -2 and 2 (an interval hypothesis). This can be useful as a Region Of Practical Equivalence test (ROPE).
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
#'   * `"Intercept_1 + Intercept_2 = 0"`: Is the sum of two intercepts zero?
#'   * ````"`cp_1_id[John]`/`cp_1_id[Erin]` = 2"````: is the varying change
#'       point for John (which is relative to `cp_1``) double that of Erin?
#' @return A data.frame with a row per hypothesis and the following columns:
#'
#'   * `hypothesis` is the hypothesis; often re-arranged to test against zero.
#'   * `mean` is the posterior mean of the left-hand side of the hypothesis.
#'   * `lower` is the lower bound of the central posterior interval of width `width`.
#'   * `upper` is the upper bound of ditto.
#'   * `p` Posterior probability.
#'       For "=" (Savage-Dickey), it is the BF converted to p.
#'       For directional hypotheses, it is the proportion of samples that returns TRUE.
#'   * `BF` Bayes Factor in favor  of the hypothesis.
#'       For "=" it is the Savage-Dickey density ratio.
#'       For directional hypotheses, it is the posterior odds divided by the
#'       prior odds.
#'
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
hypothesis = function(fit, hypotheses, width = 0.95, digits = 3, prior = FALSE) {
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_character(hypotheses)
  checkmate::assert_number(width, lower = 0, upper = 1)
  checkmate::assert_int(digits, lower = 0)
  checkmate::assert_flag(prior)

  # Loop through hypotheses and populate df_result
  df_result = data.frame()
  for (expression in hypotheses) {
    ####################
    # PREPARE FOR TEST #
    ####################
    # Check input
    n_equals = stringr::str_count(expression, "(?<!(<|>))=")
    n_directional = stringr::str_count(expression, "<|<=|>|>=")

    if (n_equals > 1)
      stop("Only one equals-test (Savage-Dickey ratio) allowed in each hypothesis: ", expression)

    if (n_equals == 1 && n_directional > 0)
      stop("Equals cannot be combined with directional tests: ", expression)

    if (n_equals + n_directional == 0)
      stop("At least one operator must be present: <, >, =, <=, or >=: ", expression)

    if (stringr::str_detect(expression, "\\[|\\]") && !stringr::str_detect(expression, "`"))
      stop("Needs `` around varying effects, e.g., `cp_1_id[2]`. Got this: ", expression)


    # If this is a single expression (does not contain & or |), we can estimate
    # the test value by putting everything on the LHS and zero on the RHS.
    if (!stringr::str_detect(expression, "\\||&")) {
      # Determine which comparator is used here
      #comparators = c("=", "<", ">"," <=", ">=")
      this_comparator = stringr::str_extract(expression, "<=|>=|<|>|=")

      # Re-arrange to LHS [comparator] 0.
      sides_split = strsplit(expression, "<=|>=|<|>|=", perl = TRUE)[[1]]
      sides_split = stringr::str_trim(sides_split)
      if (stringr::str_detect(sides_split[2], "\\+|\\-"))
        sides_split[2] = paste0("(", sides_split[2], ")")
      LHS = paste0(sides_split[1], " - ", sides_split[2])
      expression = paste0(LHS, " ", this_comparator, " 0")

      # Get effect estimate
      samples = posterior_draws(fit, prior = prior) %>%
        posterior::as_draws_df() %>%
        dplyr::mutate(effect = eval(str2lang(LHS)))

      tail_prob = (1 - width) / 2
      estimate = list(
        effect = mean(samples$effect),
        .lower = stats::quantile(samples$effect, tail_prob, names = FALSE),
        .upper = stats::quantile(samples$effect, 1 - tail_prob, names = FALSE)
      )
    } else {
      samples = posterior_draws(fit, prior = prior) %>%
        posterior::as_draws_df()

      estimate = list(effect = NA, .lower = NA, .upper = NA)
    }

    # SAVAGE-DICKEY: compute p and BF
    if (n_equals == 1) {
      if (!coda::is.mcmc.list(.subset2(fit, "mcmc_prior")) || !coda::is.mcmc.list(.subset2(fit, "mcmc_post")))
        stop("Model contains '='. Both prior and posterior samples are needed to compute Savage-Dickey density ratios. Run mcp(..., sample = 'both'")

      # Finally, let's compute those densities
      dens_prior = get_density(posterior_draws(fit, prior = TRUE), LHS, 0)
      dens_post = get_density(posterior_draws(fit), LHS, 0)
      BF = dens_post / dens_prior

      # If there is almost no density. somehow we get negative values.
      if (dens_post < 0 && dens_prior > 0)
        BF = 0
      if (dens_post > 0 && dens_prior < 0)
        BF = Inf

      p = BF / (BF + 1)
    }

    # DIRECTIONAL: compute p and BF
    if (n_directional != 0) {
      if (!coda::is.mcmc.list(.subset2(fit, "mcmc_prior")) || !coda::is.mcmc.list(.subset2(fit, "mcmc_post")))
        stop("Directional Bayes factors require both prior and posterior samples. Run mcp(..., sample = 'both').")

      # Evaluate the same hypothesis on the posterior and prior draws.
      prob_post = samples %>%
        dplyr::mutate(result = eval(str2lang(expression))) %>%  # this is where the magic happens
        dplyr::summarise(
          prob = sum(.data$result == TRUE) / dplyr::n()
        )

      prob_prior = posterior_draws(fit, prior = TRUE) %>%
        posterior::as_draws_df() %>%
        dplyr::mutate(result = eval(str2lang(expression))) %>%
        dplyr::summarise(
          prob = sum(.data$result == TRUE) / dplyr::n()
        )

      p = prob_post$prob

      # A Bayes factor is the update from prior odds to posterior odds.
      posterior_odds = prob_post$prob / (1 - prob_post$prob)
      prior_odds = prob_prior$prob / (1 - prob_prior$prob)
      BF = posterior_odds / prior_odds
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
    df_result = dplyr::bind_rows(df_result, new_row)
  }

  # Finally return
  df_result
}


#' Compute the density at a specific point.
#'
#' Used in \link{hypothesis}
#'
#' @aliases get_density
#' @keywords internal
#' @noRd
#' @param samples A posterior draws object.
#' @param LHS Expression to compute posterior
#' @param value What value to evaluate the density at
#' @return A float
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_density = function(samples, LHS, value) {
  samples = posterior::as_draws_df(samples) %>%
    dplyr::mutate(result = eval(str2lang(LHS)))
  dens = stats::density(dplyr::pull(samples, "result"), bw = "SJ")
  dens_point = stats::spline(dens$x, dens$y, xout = value)$y
  dens_point
}
