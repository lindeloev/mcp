# ABOUT: These functions compare models and/or hypotheses.
# -----------------

#' Information Criteria for Model Comparison
#'
#' Compare models using \code{\link[loo]{loo_compare}} and \code{\link[loo]{loo_model_weights}}.
#' more in \code{\link[loo]{loo}}.
#'
#' @aliases loo LOO loo.mcpfit
#' @param x An \code{\link{mcpfit}} object.
#' @param ... Further arguments passed to \code{\link[loo]{loo}}, e.g., `cores` or `save_psis`.
#' @param by_row `TRUE` calls \code{\link[loo]{loo.function}} to compute log-likelihood
#'   contributions row-by-row (observation-by-observation), which is slower but more
#'   memory efficient. `FALSE` (default) computes the full log-likelihood matrix at once.
#'   Note that both modes calculate pointwise (observation-level) PSIS-LOO cross-validation.
#' @param pointwise Deprecated alias for `by_row`.
#' @param ndraws Integer or `NULL`. Target number of posterior draws used for
#'   the log-likelihood or information criterion. Draws are balanced across
#'   chains, so the actual number may be rounded. `NULL` uses all draws.
#' @param nsamples Deprecated. Use `ndraws` instead.
#' @details Observationwise PSIS-LOO and WAIC are problematic for AR/MA models
#'   because both treat individual conditional likelihood terms as validation
#'   units. In PSIS-LOO, a held-out response also remains in the conditioning
#'   history of later terms. Prefer leave-future-out or blocked
#'   cross-validation, which are not currently implemented in mcp. When a
#'   missing response enters a later observed GARMA history, `loo()`, `waic()`,
#'   and `log_lik()` with `arma = TRUE` are unavailable: the observed-data
#'   likelihood requires integrating over that missing history, which mcp does
#'   not currently implement.
#'
#'   `loo()` and `waic()` evaluate the likelihood of the fitted model, including
#'   all group-level effects and AR/MA terms. Evaluating an information
#'   criterion with fitted components dropped post-hoc violates the PSIS
#'   identity because draws come from the full model's posterior; comparing a
#'   reduced model requires refitting it. `log_lik(fit, group_effects = FALSE)`
#'   and `log_lik(fit, arma = FALSE)` remain available as conditional or
#'   counterfactual diagnostics.
#'
#'   When `ndraws` is supplied to `loo()`, draws are balanced across chains and
#'   thinned at evenly spaced midpoint iterations. This preserves MCMC chain
#'   identities and chronological order, allowing `relative_eff()` to be
#'   computed directly.
#' @return a `loo` or `psis_loo` object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export loo
#' @export
#' @examples
#' \donttest{
#' # Define two models and sample them
#' # future::plan(future::multisession, workers = 3)  # Uncomment for parallel sampling
#' set.seed(42)
#' data = data.frame(x = seq(-1, 1, length.out = 100))
#' data$y = 1 + 2 * data$x + rnorm(100, sd = 0.3)
#' model1 = list(y ~ 1 + x)
#' model2 = list(y ~ 1)
#' fit1 = mcp(model1, data, seed = 42)
#' fit2 = mcp(model2, data, par_x = "x", seed = 42)
#'
#' # Compute LOO for each and compare (works for waic(fit) too)
#' loo1 = loo(fit1)
#' loo2 = loo(fit2)
#' loo::loo_compare(loo1, loo2)
#' }
loo.mcpfit = function(x, ..., by_row = FALSE, pointwise = lifecycle::deprecated(),
                      ndraws = NULL, nsamples = lifecycle::deprecated()) {
  if (lifecycle::is_present(pointwise)) {
    lifecycle::deprecate_soft("0.4.0", "loo(pointwise)", "loo(by_row)")
    by_row = pointwise
  }
  ndraws = legacy_arg(ndraws, nsamples, missing(ndraws), "loo.mcpfit", "ndraws", "nsamples")
  rlang::check_dots_used()  # `loo::loo()` silently ignores unknown arguments, e.g., `arma = FALSE`
  fit = x
  checkmate::assert_class(fit, "mcpfit")
  assert_loglik_garma_history(fit, fit$data, TRUE, "`loo()`")
  warn_arma_check(fit, TRUE, "information_criterion")

  # Extract posterior draws
  mcmc_post = mcmclist_draws(fit)
  ndraws = validate_loglik_ndraws(fit, ndraws)
  n_chains = length(mcmc_post)

  # Thin chains evenly across iterations because this is required by loo::relative_eff().
  if (!is.null(ndraws)) {
    n_iter_min = min(vapply(mcmc_post, nrow, integer(1)))
    iter_per_chain = max(1L, as.integer(round(ndraws / n_chains)))
    indices = pmax(1L, as.integer(round(
      (seq_len(iter_per_chain) - 0.5) * n_iter_min / iter_per_chain
    )))

    # Subset mcmc chains and any retained JAGS response imputations
    fit$mcmc_post = coda::mcmc.list(lapply(mcmc_post, function(ch) {
      coda::mcmc(ch[indices, , drop = FALSE])
    }))
    if (!is.null(fit$.internal$imputed_response)) {
      fit$.internal$imputed_response = coda::mcmc.list(lapply(
        fit$.internal$imputed_response,
        function(ch) coda::mcmc(ch[indices, , drop = FALSE])
      ))
    }

    n_draws = n_chains * iter_per_chain
    chain_id = rep(seq_len(n_chains), each = iter_per_chain)
  } else {
    n_draws = sum(vapply(mcmc_post, nrow, integer(1)))
    chain_id = rep(seq_along(mcmc_post), vapply(mcmc_post, nrow, integer(1)))
  }

  observed_rows = which(!is.na(fit$data[, mcp_columns(fit)$response]))
  if (length(observed_rows) == 0)
    stop("LOO requires at least one observed response.")

  # Matrix: Fast but memory-greedy matrix-based computation
  if (by_row == FALSE) {
    loglik = log_lik(fit, summary = FALSE)
    col_max = apply(loglik, 2, function(x) { m = max(x); if (is.finite(m)) m else 0 })
    r_eff = loo::relative_eff(exp(sweep(loglik, 2, col_max, "-")), chain_id)
    result = loo::loo(loglik, r_eff = r_eff, ...)

  # Pointwise: per-data-row computation
  } else {
    predictors = get_fit_model_tables(fit)$predictors
    ar_order = get_arma_order(predictors, "ar")
    ma_order = get_arma_order(predictors, "ma")

    # For small models, the majority of the computation time will be pp_eval overhead
    llfun = function(data_i, draws = NULL, link_fun = identity) {
      original_row = data_i$.mcp_original_row
      if (is.na(ar_order) && is.na(ma_order)) {
        loglik_samples = pp_eval(
          fit, newdata = fit$data[original_row, , drop = FALSE],
          summary = FALSE, type = "loglik"
        )
      } else {
        # AR needs only its direct response lags. MA innovations recursively
        # depend on all preceding observations, so retain the full history.
        first_row = if (is.na(ma_order)) max(1, original_row - ar_order) else 1
        data_rows = seq(first_row, original_row)
        lldata = fit$data[data_rows, ]
        loglik_samples = fit %>%
          pp_eval(newdata = lldata, summary = FALSE, type = "loglik") %>%
          dplyr::filter(.data$.data_row == max(.data$.data_row))  # last observed row
      }

      # Matrix conversion validates and explicitly orders the joint draws, so
      # every observation passed to loo uses the same draw identity.
      loglik = as.numeric(tidy_to_matrix(loglik_samples, type = ".loglik")[, 1])
      link_fun(loglik)
    }

    loo_data = data.frame(.mcp_original_row = observed_rows)

    r_eff = loo::relative_eff(
      llfun, data = loo_data, chain_id = chain_id,
      link_fun = function(x) { m = max(x); exp(if (is.finite(m)) x - m else x) },
      draws = seq_len(n_draws)
    )
    result = loo::loo.function(
      llfun, data = loo_data, r_eff = r_eff,
      draws = seq_len(n_draws), ...
    )
  }

  result
}


#' @aliases waic WAIC waic.mcpfit
#' @describeIn loo.mcpfit Computes WAIC on mcpfit objects
#' @param ... Must be empty. Reserved for future use.
#' @export waic
#' @export
waic.mcpfit = function(x, ..., ndraws = NULL, nsamples = lifecycle::deprecated()) {
  ndraws = legacy_arg(ndraws, nsamples, missing(ndraws), "waic.mcpfit", "ndraws", "nsamples")
  rlang::check_dots_empty()
  fit = x
  checkmate::assert_class(fit, "mcpfit")
  assert_loglik_garma_history(fit, fit$data, TRUE, "`waic()`")
  warn_arma_check(fit, TRUE, "information_criterion")
  mcmclist_draws(fit)
  ndraws = validate_loglik_ndraws(fit, ndraws)
  loglik = log_lik(fit, summary = FALSE, ndraws = ndraws)
  loo::waic(loglik)
}


# Validate ndraws argument for log-likelihood evaluation
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



#' Test Hypotheses Concerning Individual Parameters
#'
#' Returns posterior probabilities and Bayes Factors for flexible hypotheses involving
#' model parameters. The documentation for the argument `hypotheses` below
#' shows examples of how to specify hypotheses, and [read worked examples on the mcp website](https://lindeloev.github.io/mcp/articles/comparison.html).
#' For directional hypotheses, `hypothesis` executes the hypothesis string in
#' a data-frame environment and summarises the proportion of posterior and
#' prior draws where the expression evaluates to TRUE. The Bayes factor is
#' the posterior odds divided by the prior odds. For equality hypotheses, a
#' Savage-Dickey ratio is computed. Both kinds of Bayes factor require prior
#' draws, so remember `mcp(..., sample = "both")`. This function is heavily inspired by the
#' `hypothesis` function from the `brms` package.
#' When `prior = TRUE`, the summary is based on prior draws and `BF` is `NA`.
#'
#' @aliases hypothesis hypothesis.mcpfit
#' @inheritParams summary.mcpfit
#' @param fit An \code{\link{mcpfit}} object.
#' @param hypotheses String representation of a logical test involving model parameters.
#'   Takes R code that evaluates to TRUE or FALSE in a vectorized way.
#'
#'   **Directional hypotheses** are specified using <, >, <=, or >=. `hypothesis`
#'   returns the posterior probability \eqn{P(H \mid \text{data})}{P(H | data)}
#'   and the Bayes factor in favor of the stated hypothesis \eqn{H}:
#'   \deqn{\text{BF}_{10} = \frac{P(H \mid \text{data}) / (1 - P(H \mid \text{data}))}{P(H) / (1 - P(H))}}
#'   where \eqn{P(H)} is the prior probability of \eqn{H}. The Bayes factor requires both prior and posterior
#'   draws from `mcp(sample = "both")`. For example:
#'
#'   * `"cp_1 > 30"`:  the first change point is above 30.
#'   * `"Intercept_1 > Intercept_2"`: the intercept is greater in segment 1 than 2.
#'   * `"x_2 - x_1 <= 3"`: the difference between slope 1 and 2 is less
#'       than or equal to 3.
#'   * `"Intercept_1 > -2 & Intercept_1 < 2"`: Intercept_1 is between -2 and 2 (an interval hypothesis). This can be useful as a Region Of Practical Equivalence test (ROPE).
#'   * `"cp_1^2 < 30 | (log(x_1) + log(x_2)) > 5"`: be creative.
#'   * \code{"`cp_1_id[1]` > `cp_1_id[2]`"}: id1 is greater than id2, as estimated
#'       through the group-level change-point deviation for `id` at change point 1 (starting segment 2).
#'       Note that \code{``} are required when using `[i]`.
#'
#'   **Equality hypotheses** use the equal sign (=) and a Savage-Dickey density
#'   ratio: posterior density divided by prior density at the tested point equality
#' 
#'   \eqn{\theta = \theta_0}:
#'   \deqn{\text{BF}_{01} = \frac{p(\theta = \theta_0 \mid \text{data})}{p(\theta = \theta_0)}}
#' 
#'   where \eqn{\theta} is the evaluated parameter (or affine contrast), \eqn{\theta_0} is the hypothesized null value, \eqn{p(\theta = \theta_0 \mid \text{data})}{p(theta = theta_0 | data)} is the posterior density, and \eqn{p(\theta = \theta_0)}{p(theta = theta_0)} is the prior density. This is a Bayes factor for a nested point-null model against the fitted
#'   continuous model (\eqn{\text{BF}_{10} = 1 / \text{BF}_{01}}). Prior and posterior draws are required, using
#'   `mcp(sample = "both")`.
#' 
#'   The point-null model's nuisance prior is the fitted model's conditional 
#'   prior at the equality. Equality tests are limited to named scalar parameters 
#'   and affine contrasts. 
#' 
#'   Examples:
#'
#'   * `"cp_1 = 30"`: compare the point-null model `cp_1 = 30` with the
#'       continuous alternative.
#'   * `"Intercept_1 - Intercept_2 = 0"`: compare equal segment intercepts
#'       with the continuous alternative.
#' @return A tibble with a row per hypothesis and the following columns:
#'
#'   * `hypothesis` is the hypothesis; often re-arranged to test against zero.
#'   * `mean` is the posterior mean of the left-hand side of the hypothesis, or
#'       the prior mean when `prior = TRUE`.
#'   * `lower` is the lower bound of the central posterior interval of width `width`,
#'       or the corresponding prior interval when `prior = TRUE`.
#'   * `upper` is the upper bound of ditto.
#'   * `prob` is the posterior probability of a directional hypothesis, or the prior
#'       probability when `prior = TRUE`. It is `NA`
#'       for equality hypotheses, which compare models rather than an event
#'       within the fitted model.
#'   * `BF` Bayes Factor in favor  of the hypothesis.
#'       For "=" it is the Savage-Dickey density ratio.
#'       For directional hypotheses, it is the posterior odds divided by the
#'       prior odds. It is `Inf` (or `0`) when all (or no) posterior draws satisfy
#'       the hypothesis, i.e., the evidence exceeds what the number of draws can
#'       resolve; increase `iter` for a finite estimate. It is `NA` when `prior = TRUE`,
#'       when prior draws are not available (e.g., `sample = "post"`), or when all
#'       or no prior draws satisfy the hypothesis.
#'
#' @examples
#' # demo_fit contains both posterior and prior draws
#' # A directional hypothesis returns its posterior probability and Bayes factor
#' hypothesis(demo_fit, "cp_1 > 30")
#'
#' # Combine directional statements for an interval (a ROPE-style hypothesis)
#' hypothesis(demo_fit, "cp_1 > 20 & cp_1 < 30")
#'
#' # Evaluate several directional hypotheses at once
#' hypothesis(demo_fit, c("cp_1 > 30", "cp_2 > 70"))
#'
#' # Inspect the corresponding prior probability without a Bayes factor
#' hypothesis(demo_fit, "cp_1 > 30", prior = TRUE)
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
hypothesis = function(fit, ...) UseMethod("hypothesis")


#' @rdname hypothesis
#' @export
hypothesis.mcpfit = function(fit, hypotheses, width = 0.95, prior = FALSE, ...) {
  rlang::check_dots_empty()
  checkmate::assert_class(fit, "mcpfit")
  checkmate::assert_character(hypotheses)
  legacy_parameter_names(hypotheses, "hypothesis()")
  checkmate::assert_number(width, lower = 0, upper = 1)
  checkmate::assert_flag(prior)

  # Selected draws, plus prior draws for Bayes factors when the posterior is selected
  draws = posterior::as_draws_df(posterior_draws(fit, prior = prior))
  has_prior_draws = !prior && !is.null(mcmclist_draws(fit, prior = TRUE, error = FALSE))
  draws_prior = if (has_prior_draws) posterior::as_draws_df(posterior_draws(fit, prior = TRUE))

  # Loop through hypotheses and populate df_result
  df_result = tibble::tibble()
  for (expression in hypotheses) {
    ####################
    # PREPARE FOR TEST #
    ####################
    # Check input
    if (stringr::str_detect(expression, "!=|=="))
      stop("`!=` and `==` are not supported. Use `=` for equality hypotheses; the Bayes factor for the inequality is 1 / BF: ", expression, call. = FALSE)

    n_equals = stringr::str_count(expression, "(?<!(<|>))=")
    n_directional = stringr::str_count(expression, "<|<=|>|>=")

    if (n_equals > 1)
      stop("Only one equals-test (Savage-Dickey ratio) allowed in each hypothesis: ", expression)

    if (n_equals == 1 && n_directional > 0)
      stop("Equals cannot be combined with directional tests: ", expression)

    if (n_equals + n_directional == 0)
      stop("At least one operator must be present: <, >, =, <=, or >=: ", expression)

    if (stringr::str_detect(expression, "\\[|\\]") && !stringr::str_detect(expression, "`"))
      stop("Needs `` around group-level deviations, e.g., `cp_1_id[2]`. Got this: ", expression)

    if (n_equals == 1)
      validate_savage_dickey_expression(expression, posterior::variables(draws))


    # If this is a single expression (does not contain & or |), we can estimate
    # the test value by putting everything on the LHS and zero on the RHS.
    if (!stringr::str_detect(expression, "\\||&")) {
      # Determine which comparator is used here
      this_comparator = stringr::str_extract(expression, "<=|>=|<|>|=")

      # Re-arrange to LHS [comparator] 0.
      sides_split = strsplit(expression, "<=|>=|<|>|=", perl = TRUE)[[1]]
      sides_split = stringr::str_trim(sides_split)
      if (stringr::str_detect(sides_split[2], "\\+|\\-"))
        sides_split[2] = paste0("(", sides_split[2], ")")
      LHS = paste0(sides_split[1], " - ", sides_split[2])
      expression = paste0(LHS, " ", this_comparator, " 0")

      # Get effect estimate
      LHS_expr = rlang::parse_expr(LHS)
      effect = rlang::eval_tidy(LHS_expr, data = draws)

      tail_prob = (1 - width) / 2
      estimate = list(
        effect = mean(effect),
        .lower = stats::quantile(effect, tail_prob, names = FALSE),
        .upper = stats::quantile(effect, 1 - tail_prob, names = FALSE)
      )
    } else {
      estimate = list(effect = NA, .lower = NA, .upper = NA)
    }

    # SAVAGE-DICKEY: compute BF. Equality hypotheses always have an LHS (see validate_savage_dickey_expression()).
    if (n_equals == 1) {
      if (prior) {
        BF = NA_real_
      } else {
        if (!has_prior_draws)
          stop("Model contains '='. Both prior and posterior draws are needed to compute Savage-Dickey density ratios. Run mcp(..., sample = 'both').")

        # Check that that both prior contrast and posterior contrast are distributed (not constants)
        # Protects against "x_1 - x_1 = 1" and prior = list(x_1 = 3)
        prior_values = rlang::eval_tidy(LHS_expr, data = draws_prior)
        post_values = effect
        if (diff(range(prior_values)) == 0)
          stop("Cannot compute Savage-Dickey density ratio because the prior contrast is constant: ", expression, call. = FALSE)
        if (diff(range(post_values)) == 0)
          stop("Cannot compute Savage-Dickey density ratio because the posterior contrast is constant: ", expression, call. = FALSE)

        # Warning if testing against default priors
        hypothesis_pars = sub("\\[.*\\]$", "", all.vars(rlang::parse_expr(expression)))
        is_default_pars = fit$.internal$prior_table$parameter %in% hypothesis_pars & fit$.internal$prior_table$source == "default"
        default_pars = fit$.internal$prior_table$parameter[is_default_pars]
        if (length(default_pars) > 0) {
          warning(
            "Savage-Dickey Bayes factor was computed using default prior(s) for ",
            and_collapse(paste0("`", default_pars, "`")),
            ". Point Bayes factors are sensitive to the prior distribution; consider specifying informed priors.",
            call. = FALSE
          )
        }

        if (is_sparse_tail(prior_values, 0) || is_sparse_tail(post_values, 0)) {
          warning(
            "The tested value is in a sparse tail of the prior or posterior draws; the Savage-Dickey estimate may be unreliable.",
            call. = FALSE
          )
        }

        dens_prior = get_density(prior_values, 0)
        dens_post = get_density(post_values, 0)
        BF = dens_post / dens_prior
      }

      prob = NA_real_
    }

    # DIRECTIONAL: compute probability and BF
    if (n_directional != 0) {
      expr_parsed = rlang::parse_expr(expression)
      result = rlang::eval_tidy(expr_parsed, data = draws)
      prob = mean(result == TRUE)

      if (has_prior_draws) {
        prior_result = rlang::eval_tidy(expr_parsed, data = draws_prior)
        prior_prob = mean(prior_result == TRUE)

        # A Bayes factor is the update from prior odds to posterior odds.
        # It is undefined if the prior already implies or excludes the hypothesis.
        if (prior_prob %in% c(0, 1)) {
          warning("No prior draws ", if (prior_prob == 1) "violate" else "satisfy", " `", expression,
            "`, so its Bayes factor is undefined. Is it implied or excluded by the model?", call. = FALSE)
          BF = NA_real_
        } else {
          BF = (prob / (1 - prob)) / (prior_prob / (1 - prior_prob))
        }
      } else {
        BF = NA_real_
      }
    }

    # Add to list
    new_row = tibble::tibble(
      hypothesis = stringr::str_trim(expression),
      mean = estimate$effect,
      lower = estimate$.lower,
      upper = estimate$.upper,
      prob = prob,
      BF = BF
    )
    df_result = dplyr::bind_rows(df_result, new_row)
  }

  # Finally return
  df_result
}


# Validate the deliberately small expression language used by Savage-Dickey tests.
validate_savage_dickey_expression = function(expression, parameters) {
  sides = strsplit(expression, "=", fixed = TRUE)[[1]]
  sides = stringr::str_trim(sides)

  inspect = function(x) {
    if (is.numeric(x) && length(x) == 1)
      return(list(valid = TRUE, constant = TRUE))

    if (is.symbol(x)) {
      is_parameter = as.character(x) %in% parameters
      return(list(valid = is_parameter, constant = FALSE))
    }

    if (!is.call(x))
      return(list(valid = FALSE, constant = FALSE))

    operator = as.character(x[[1]])
    arguments = as.list(x)[-1]

    if (operator == "(" && length(arguments) == 1)
      return(inspect(arguments[[1]]))

    if (operator %in% c("+", "-") && length(arguments) == 1)
      return(inspect(arguments[[1]]))

    if (operator %in% c("+", "-") && length(arguments) == 2) {
      left = inspect(arguments[[1]])
      right = inspect(arguments[[2]])
      return(list(
        valid = left$valid && right$valid,
        constant = left$constant && right$constant
      ))
    }

    if (operator == "*" && length(arguments) == 2) {
      left = inspect(arguments[[1]])
      right = inspect(arguments[[2]])
      return(list(
        valid = left$valid && right$valid && (left$constant || right$constant),
        constant = left$constant && right$constant
      ))
    }

    list(valid = FALSE, constant = FALSE)
  }

  parsed = lapply(sides, rlang::parse_expr)
  checked = lapply(parsed, inspect)
  valid = length(sides) == 2 && all(vapply(checked, `[[`, logical(1), "valid"))
  has_parameter = any(!vapply(checked, `[[`, logical(1), "constant"))

  if (!valid || !has_parameter) {
    stop(
      "Savage-Dickey equality hypotheses must use a named scalar parameter or affine contrast. ",
      "Use only parameter names, numeric constants, +, -, parentheses, and multiplication by a numeric constant: ",
      expression,
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Compute the density at a specific point.
#'
#' Used in \link{hypothesis}
#'
#' @param x Numeric draws.
#' @param value What value to evaluate the density at
#' @return A float
#' @noRd
get_density = function(x, value) {
  if (diff(range(x)) == 0) {
    stop("Cannot estimate density for constant values.", call. = FALSE)
  }
  bandwidth = tryCatch(stats::bw.SJ(x), error = function(e) stats::bw.nrd0(x))
  mean(stats::dnorm(value, mean = x, sd = bandwidth))
}


# Is the tested value outside the central 98% of draws?
is_sparse_tail = function(x, value) {
  mean(x <= value) < 0.01 || mean(x >= value) < 0.01
}
