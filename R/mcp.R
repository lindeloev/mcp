#' Fit Multiple Linear Segments And Their Change Points
#'
#' Given a model (a list of segment formulas), `mcp` infers the posterior
#' distributions of the parameters of each segment as well as the change points
#' between segments. See details or [the mcp website](https://lindeloev.github.io/mcp/).
#'
#' @aliases mcp
#' @param data Table-like data in long format (data.frame, tibble, data.table, etc.)
#'   with syntactic column names.
#'   Missing values in the response variable are imputed using the posterior predictive.
#'   \code{\link{fitted.mcpfit}} or \code{\link{predict.mcpfit}} details how to see the imputed values.
#' @param model A list of formulas - one for each segment. The general format is
#'   `response ~ cp ~ predictors` (e.g., `y ~ 1 ~ 1 + x`), except the first segment
#'   has no change point and uses `response ~ predictors`. The response and
#'   change-point parts can be omitted (`cp ~ predictor` assumes the same
#'   response; `~ predictor` assumes an intercept-only change point). Non-$x$ terms persist
#'   into later segments until replaced or removed by an intercept reset (see details). See examples on the
#'   [mcp website](https://lindeloev.github.io/mcp/).
#'
#'   **1. Response (segment 1 only):**
#'   * `y ~ ...`: Standard continuous or count response (Gaussian, Poisson, Bernoulli).
#'   * `successes | trials(total) ~ ...`: Binomial response (`family = binomial()`).
#'   * `y | weights(w) ~ ...`: Observation log-likelihood weights (multiplies each observation's
#'     log-likelihood contribution by `w > 0`; affects posterior inference and `log_lik()`, but
#'     not predictions).
#'   * `y | trials(total) + weights(w) ~ ...`: Combine response auxiliaries using `+`.
#'
#'   **2. Change-point modeling (`cp`, segments 2+):**
#'   * `~ 1 ~ ...` (or omitted, e.g., `~ x`): Population-level change point (default).
#'   * `1 + (1 | id) ~ ...`: Group-level change-point deviations around the population change point.
#'     [Read more](https://lindeloev.github.io/mcp/articles/group_effects.html).
#'
#'   **3. Regression formula (all segments):** [Read more](https://lindeloev.github.io/mcp/articles/formulas.html)
#'   * `~ 1 + x`: Disjoined slope with a new segment intercept.
#'   * `~ 0 + x`: Joined slope (no intercept; continuous from previous segment).
#'   * `~ 1`: Plateau (intercept only, no slope).
#'   * `~ x:group + I(x^2) + exp(z)`: Extended terms, interactions, and R-side bases
#'     (`scale()`, `poly()`, `splines::ns()`). Bases are evaluated before sampling and reused for `newdata`.
#'   * `~ 1 + (1 | id)`: Group-level intercepts (or `(1 + x || id)` for independent slopes and intercepts).
#'     [Read more](https://lindeloev.github.io/mcp/articles/group_effects.html).
#'   * `~ sigma(1 + x)`: Distributional parameters on the link scale (e.g., log residual SD).
#'     [Read more](https://lindeloev.github.io/mcp/articles/dpar.html).
#'   * `~ ar(1) + ma(1)`: Autoregressive and moving-average time-series residuals on the link scale
#'     (accepts regression formulas, `series = id`, and `boundary`; use `ar(0)` or `ma(0)` to turn off in later segments).
#'     [Read more](https://lindeloev.github.io/mcp/articles/arma.html).
#'
#' @param prior Named list. Names are parameter names (`cp_i`, `Intercept_i`, `xvar_i`,
#'  `sigma_1`, etc.) and the values are either
#'
#'  * A distribution in `mcp`'s JAGS-string syntax (e.g.,
#'      `Intercept_1 = "dnorm(0, 1) T(0,)"`) indicating a
#'      conventional prior distribution. Data-calibrated, regularizing defaults
#'      are used where priors are not specified. These are designed for stable
#'      estimation and prediction, but should be justified before hypothesis testing.
#'      `mcp` uses conventional distribution scales rather than JAGS precision:
#'      SD for `dnorm()`, scale for `dt()`, `ddexp()`, and `dlogis()`, and
#'      log-SD for `dlnorm()`. See details.
#'  * A numerical value (e.g., `Intercept_1 = -2.1`) indicating a fixed value.
#'  * A model parameter name (e.g., `Intercept_2 = "Intercept_1"`), indicating that this parameter is shared -
#'      typically between segments. If two group-level deviations are shared this way,
#'      they will need to have the same grouping variable.
#'  * The default prior on change points is `dirichlet(1)` (uniform order statistics).
#'      For a single change point, this is the Beta(1, 1) / Uniform distribution over `[min(x), max(x)]`.
#'      For multiple change points, it corresponds to a flat Dirichlet distribution over segment lengths.
#'      You can also explicitly set `cp_i = "dirichlet(alpha)"` with the same positive `alpha` for all
#'      change points to regularize spacing (`alpha > 1` penalizes change points from occurring close together,
#'      while `alpha < 1` favors clustering). Under the hood, this is parameterized as an exact sequential
#'      stick-breaking Beta chain for fast and robust sampling.
#'      [Read more](https://lindeloev.github.io/mcp/articles/priors.html).
#' @param family A supported family: `gaussian()`, `binomial()`, `bernoulli()`,
#'   `poisson()`, or `negbinomial()`, with a supported link function; e.g.,
#'   `gaussian(link = "log")`.
#' @param par_x String (default: `NULL` which is auto-detect).
#' @param sample One of
#'   * `"post"`: Sample the posterior.
#'   * `"prior"`: Sample only the prior. Use `prior = TRUE` explicitly in
#'       plots, summaries, and other methods to select prior draws.
#'   * `"both"`: Sample both prior and posterior. Plots, summaries, etc.
#'       will default to using the posterior. Use `prior = TRUE` to select prior draws.
#'   * `"none"` or `FALSE`: Do not sample. Returns an mcpfit
#'       object without sample. This is useful if you only want to check
#'       prior strings (`fit$prior`), the JAGS model (`fit$jags_code`), etc.
#' @param cores Deprecated and ignored. Configure parallel processing with a
#'   [future][future::plan] plan instead, for example
#'   `future::plan(future::multisession, workers = 3)`. With the default future
#'   plan, chains are sampled sequentially. The argument remains available for
#'   backwards compatibility.
#' @param chains Positive integer. Number of chains to run.
#' @param iter Positive integer. Number of post-warmup draws from each chain.
#'   The total number of draws is `iter * chains`.
#' @param warmup Positive integer. Number of initial iterations per chain which
#'   are discarded before sampling. Set higher if needed for sampler adaptation;
#'   use diagnostics to assess convergence.
#' @param adapt Deprecated; use `warmup` instead.
#' @param inits A list of initial values for the parameters. This can be useful
#'   if a model fails to converge. Read more in \code{\link[rjags]{jags.model}}.
#'   Defaults to `NULL`, i.e., no inits.
#' @param jags_code String. Pass JAGS code to `mcp` to use directly. This is useful if
#'   you want to tweak the code in `fit$jags_code` and run it within the `mcp`
#'   framework. R-side simulation and prediction methods continue to use the
#'   mcp-default formulas (with warning), so they may no longer match custom JAGS code.
#' @param seed `NULL` or a positive integer. Seed for the JAGS random-number
#'   generators.
#' @param diagnostics Named list of diagnostic warning thresholds. Available
#'   elements are `rhat = 1.01`, `ess_bulk = 400`, `ess_tail = 400`,
#'   `ar = 0.10`, and `ma = 0.10`. An empty list uses these defaults; a partial
#'   list overrides only the supplied values. Set an element to `NULL` to disable
#'   that diagnostic, or use `FALSE` to disable all configurable diagnostic
#'   warnings. In [summary.mcpfit()], `NULL` inherits the settings used to fit
#'   the model, while a list or `FALSE` overrides the diagnostic footer.
#' @param quiet Logical. Suppress routine JAGS output and mcp sampling-status
#'   messages? Defaults to `FALSE`.
#'
#' @section The mcp model:
#' Consider the following model which you can find in `demo_fit` and `mcp_example("demo")`:
#'
#' ```r
#' model = list(
#'   response ~ 1,  # Plateau in the first segment (Intercept_1)
#'   ~ 0 + time,    # Joined slope (time_2) in segment 2 which starts at cp_1
#'   ~ 1 + time     # Disjoined slope (Intercept_3, time_3) at cp_2
#' )
#' ```
#'
#' \if{html}{\figure{mcp_demo.png}{options: width="500" alt="Fitted 3-segment mcp model with a plateau, joined slope, and disjoined slope"}}
#' \if{latex}{\figure{mcp_demo.png}{options: width=5in}}
#' \if{text}{\figure{mcp_demo.png}{[Fitted 3-segment mcp model with a plateau, joined slope, and disjoined slope]}}
#'
#' This model has \eqn{K=3} segments separated by \eqn{K-1=2} change points: \eqn{\tau_1} and \eqn{\tau_2}.
#'
#' More generally, an `mcp` model divides a continuous predictor \eqn{x} into \eqn{K} segments separated by
#' ordered change points \eqn{\tau_1 < \dots < \tau_{K-1}}. In each segment \eqn{k \in \{1, \dots, K\}},
#' the linear predictor \eqn{\eta_i} is evaluated directly from the segment-local distance \eqn{(x_i - \tau_{k-1})}:
#'
#' \deqn{\eta_i = \alpha_k + \beta_{k,1} (x_i - \tau_{k-1}) \quad (\text{with } \tau_0 = 0)}
#'
#' where the segment-start level \eqn{\alpha_k} is freely estimated for the first and disjoined segments, as in non-segmented regression, and determined by continuity for joined segments:
#'
#' \deqn{\alpha_k = \begin{cases}
#'   \beta_{k,0}, & \text{Disjoined segments } (\sim \texttt{1 + x}, \text{ including } k = 1) \\
#'   \alpha_{k-1} + \beta_{k-1,1} (\tau_{k-1} - \tau_{k-2}), & \text{Joined segments } (k \ge 2, \sim \texttt{0 + x})
#' \end{cases}}
#'
#' Here, \eqn{\beta_{k,0}} is the segment-start intercept, and \eqn{\beta_{k,1}} is the slope on \eqn{x}. In all segments, estimated slope and intercept parameters are absolute values (not changes relative to the preceding segment).
#'
#' If additional continuous covariates or categorical factors are included (e.g., `+ z + group`),
#' they enter additively on their original scale (\eqn{\dots + \sum \gamma_{k,j} z_{j,i}} for covariate \eqn{j}); only bare
#' change-point predictor terms \eqn{x} and \code{I(x^k)} are converted to segment-local coordinates.
#'
#' The inverse-linked parameter is \eqn{\mu_i = g^{-1}(\eta_i)} via link function \eqn{g(\mu_i) = \eta_i},
#' representing the expected response for most families (or the success probability for binomial models,
#' where the expected count is \eqn{n_i \mu_i}). Distributional parameters (\code{sigma()}, \code{shape()}, etc.)
#' and autoregressive terms (\code{ar()}, \code{ma()}) follow this exact same segmented structure on their respective link scales.
#' See more details on the `mcp` model in [mcp-package] and on the [mcp website](https://lindeloev.github.io/mcp/articles/formulas.html).
#'
#' @section Time-series residuals (link-scale observation-driven GARMA):
#' Autoregressive (`ar(p)`) and moving-average (`ma(q)`) terms define a finite conditional recurrence
#' on the link scale (generalized autoregressive moving-average, GARMA). They support Gaussian (`identity`),
#' binomial (`logit`), Bernoulli (`logit`), Poisson (`log`), and negative-binomial (`log`) families.
#' If \eqn{\eta^{\text{reg}}_t} is the ordinary regression predictor from the segment formulas and \eqn{\eta_t} is the
#' predictor including serial dependence, the recurrence decomposes into components:
#'
#' \deqn{\begin{aligned}
#'   \text{AR}_t &= \sum_{j=1}^{p} \phi_{j,t} \left[g(y^*_{t-j}) - \eta^{\text{reg}}_{t-j}\right] \\
#'   \text{MA}_t &= \sum_{k=1}^{q} \theta_{k,t} \left[g(y^*_{t-k}) - \eta_{t-k}\right] \\
#'   \eta_t &= \eta^{\text{reg}}_t + \text{AR}_t + \text{MA}_t
#' \end{aligned}}
#'
#' where \eqn{\phi_{j,t}} is the lag-\eqn{j} autoregressive (AR) coefficient at time \eqn{t},
#' \eqn{\theta_{k,t}} is the lag-\eqn{k} moving-average (MA) coefficient at time \eqn{t},
#' \eqn{g(\cdot)} is the link function, and \eqn{y^*_t} is the boundary-constrained observation with pseudo-count \eqn{b} (set via argument \code{boundary = 0.1} in \code{ar()} / \code{ma()}) to keep residuals finite on the link scale:
#' * **Gaussian:** \eqn{y^*_t = y_t}.
#' * **Poisson / Negative Binomial:** \eqn{y^*_t = \max(y_t, b)} to prevent \eqn{\log(0)}. Here \eqn{b} replaces zero counts with a small positive count.
#' * **Binomial / Bernoulli:** \eqn{y^*_t = \min(\max(y_t, b), n_t - b) / n_t}, where \eqn{y_t} is observed successes, \eqn{n_t} is the number of trials (\eqn{n_t = 1} for Bernoulli), and \eqn{b} clamps counts to the interval \eqn{[b, n_t - b]} before converting to a rate, preventing \eqn{\text{logit}(0)} and \eqn{\text{logit}(1)}.
#'
#' Implications:
#' * For an \eqn{N}-order component, the last \eqn{N} values *before* the segment onset are input to the first \eqn{\eta_t} in the segment.
#' * AR and MA components persist into later segments until replaced or turned off via \code{ar(0)} or \code{ma(0)}.
#' * AR coefficients are not jointly constrained to stationarity; nor MA coefficients to invertibility.
#' * See [the arma article](https://lindeloev.github.io/mcp/articles/arma.html) for more details.
#'
#' @section Notes on priors:
#'
#'   * *Ordered change point priors:* Default population-level `cp_i` priors are ordered and the ordering is imposed through the priors. For user-defined priors,
#'       `mcp` adds truncation (e.g., `T(cp_1, )`) only when the prior has neither
#'       explicit truncation nor an inherently bounded form such as `dunif()` or
#'       `dirichlet()`.
#'   * *Data-dependent terms:* If `mcp` encounters a data-dependent term like `min(time)`, `max(time)`, `median(response)`, or `mad(response)` in the prior string, they are resolved from the model data so a numerical value is passed to JAGS. The following terms are also allowed: `n_segments()` and `n_cp()`.
#'       The older constants `MINX`, `MAXX`, `MEANX`, `SDX`, `MINY`, `MAXY`, `MEANY`, `SDY`, and `N_CP` remain accepted with a deprecation warning.
#'   * *Group-level change points:* Group-specific locations follow a hierarchical
#'       normal distribution around their population change point, truncated so
#'       that realized locations remain in the observed range and ordered.
#'   * *Parameterization:* Prior strings use conventional scale parameterizations. `mcp` converts
#'       these to the parameterization required by JAGS when generating code:
#'       inverse variance for `dnorm()`, `dt()`, and `dlnorm()`, and inverse
#'       scale for `ddexp()` and `dlogis()`. Use
#'       `prior_summary(fit)` for resolved priors and
#'       `prior_summary(fit, verbose = TRUE)` for their rules and descriptions.
#' @references
#' * Lindeløv, J. K. (2020). mcp: An R Package for Regression With Multiple Change Points.
#'   *OSF Preprints*. [doi:10.31219/osf.io/fzqxv](https://doi.org/10.31219/osf.io/fzqxv)
#'   Introduces the `mcp` package, formula syntax, default priors, and workflow for regression
#'   with multiple change points across generalized linear and time-series models. Newer-than-2020 versions of
#'   the paper may be available at that link. Please cite the newest version.
#' * Carlin, B. P., Gelfand, A. E., & Smith, A. F. (1992). Hierarchical Bayesian Analysis of
#'   Changepoint Problems. *Applied Statistics*, 41(2), 389–405. [doi:10.2307/2347570](https://doi.org/10.2307/2347570)
#'   Introduced the Gibbs sampling approach for continuous change points in segmented regression
#'   and hierarchical models, providing the computational foundation used by BUGS and JAGS.
#' @return An \code{\link{mcpfit}} object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' \donttest{
#' # Define the segments using formulas. A change point is estimated between each formula.
#' model = list(
#'   response ~ 1,  # Plateau in the first segment (Intercept_1)
#'   ~ 0 + time,    # Joined slope (time_2) in segment 2 which starts at cp_1
#'   ~ 1 + time     # Disjoined slope (Intercept_3, time_3) at cp_2
#' )
#'
#' # Fit it and sample the prior too.
#' # future::plan(future::multisession, workers = 3)  # Uncomment for parallel sampling
#' data = mcp_example_data("demo")  # Simulated data example
#' demo_fit = mcp(model, data = data, sample = "both", seed = 42)
#'
#' # See parameter estimates
#' summary(demo_fit)
#'
#' # Visual inspection of the results
#' plot(demo_fit)  # Visualization of model fit/predictions
#' plot_pars(demo_fit)  # Parameter distributions
#' pp_check(demo_fit)  # Prior/Posterior predictive checks
#'
#' # Test a hypothesis
#' hypothesis(demo_fit, "cp_1 > 10")
#'
#' # Make predictions
#' head(fitted(demo_fit))
#' head(predict(demo_fit))
#' head(predict(demo_fit, newdata = data.frame(time = c(55.545, 80, 132))))
#'
#' # Compare to a one-intercept-only model (no change points) with default prior
#' model_null = list(response ~ 1)
#' fit_null = mcp(model_null, data = data, par_x = "time")  # fit another model here
#' demo_loo = loo(demo_fit)
#' null_loo = loo(fit_null)
#' loo::loo_compare(demo_loo, null_loo)
#'
#' # Inspect the prior. Useful for prior predictive checks.
#' summary(demo_fit, prior = TRUE)
#' plot(demo_fit, prior = TRUE)
#'
#' # Show all priors. Default priors are added where you don't provide any
#' prior_summary(demo_fit)
#'
#' # Set priors and re-run
#' prior = list(
#'   Intercept_1 = 15,
#'   time_2 = "dt(0, 2, 1) T(0, )",  # t-dist slope. Truncated to positive.
#'   cp_2 = "dunif(cp_1, 80)",       # change point to segment 2 > cp_1 and < 80.
#'   Intercept_3 = "Intercept_1"     # Shared intercept between segment 1 and 3
#' )
#'
#' fit3 = mcp(model, data = data, prior = prior, warmup = 2000, iter = 6000, seed = 42)
#'
#' # Show the JAGS model
#' demo_fit$jags_code
#' }
mcp = function(model,
               data,
               prior = list(),
               family = gaussian(),
               par_x = NULL,
               sample = "post",
               cores = NULL,
               chains = 3,
               iter = 3000,
               warmup = 1500,
               adapt = lifecycle::deprecated(),
               inits = NULL,
               jags_code = NULL,
               seed = NULL,
               diagnostics = list(),
               quiet = FALSE) {
  custom_jags_code = !is.null(jags_code)

  matched_call = match.call()

  ################
  # CHECK INPUTS #
  ################
  # Check model
  checkmate::assert_true(is.mcpmodel(model), .var.name = "model")
  assert_rel(model)

  # Check data and data-model correspondence
  if (missing(data) || is.null(data) || !is.data.frame(data))
    stop("`data` is required in mcp() since mcp v0.4.0. Passing data = NULL or omitting data is no longer supported.", call. = FALSE)

  checkmate::assert_data_frame(data)
  non_syntactic = names(data)[make.names(names(data)) != names(data)]
  if (length(non_syntactic) > 0)
    stop(
      "`data` has non-syntactic column name(s): ",
      paste0("`", non_syntactic, "`", collapse = ", "), ". Rename them before fitting.",
      call. = FALSE
    )
  data = data.frame(data)
  series = get_arma_series(model)
  assert_arma_series(data, series)

  checkmate::assert_string(par_x, null.ok = TRUE)
  par_x = get_par_x(model, data, par_x)
  rhs_vars = get_rhs_vars(model)
  model_vars = unique(c(get_model_vars(model), par_x, series))
  assert_model_data(data, par_x, rhs_vars)
  assert_data_cols(data, cols = model_vars, fail_funcs = c(is.infinite))
  data = data[, model_vars]  # Remove unused data
  response_var = get_model_vars(model)[1]
  if (response_var %in% names(data) && is.logical(data[[response_var]]))
    data[[response_var]] = as.numeric(data[[response_var]])

  # Plain named lists are the v0.4 JAGS-string prior format. Future classed
  # sampler-agnostic prior objects can dispatch before this legacy path.
  validate_prior_v1(prior)

  # Transform family to mcpfamily
  if (!is.family(family) && !is.mcpfamily(family))
    stop("`family` is not a valid family or mcpfamily. Should be gaussian(), binomial(), mcpfamily(gaussian(link = 'log')), etc.")

  family = mcpfamily(family)

  # More checking...
  checkmate::assert(
    checkmate::check_choice(sample, c("post", "prior", "both", "none")),
    checkmate::check_false(sample),
    .var.name = "sample"
  )
  if (!is.null(cores)) {
    if (!identical(cores, "all"))
      checkmate::assert_int(cores, lower = 1)

    cores_details = paste0(
      "`cores` is ignored. Parallel processing is now controlled by the active ",
      "future plan. For example, call ",
      "`future::plan(future::multisession, workers = 3)` before `mcp()`, and ",
      "`future::plan(future::sequential)` when the workers are no longer needed."
    )
    if (identical(cores, "all") || (is.numeric(cores) && cores > 1))
      cores_details = paste0(
        "Setting `cores` above one no longer enables parallel processing. ",
        cores_details
      )

    lifecycle::deprecate_warn(
      when = "0.4.0",
      what = "mcp(cores)",
      details = cores_details
    )
  }

  if (lifecycle::is_present(adapt)) {
    if (!missing(warmup))
      stop("Supply only one of `warmup` and deprecated `adapt`.", call. = FALSE)

    lifecycle::deprecate_soft(
      when = "0.4.0",
      what = "mcp(adapt)",
      with = "mcp(warmup)"
    )
    warmup = adapt
  }

  checkmate::assert_int(chains, lower = 1)
  checkmate::assert_int(iter, lower = 1)
  checkmate::assert_int(warmup, lower = 1)
  checkmate::assert_list(inits, null.ok = TRUE)
  checkmate::assert_int(seed, lower = 1, null.ok = TRUE)
  diagnostics = resolve_diagnostics(diagnostics)
  if (missing(quiet)) {
    quiet = isTRUE(getOption("mcp.quiet", FALSE)) || identical(Sys.getenv("IN_PKGDOWN"), "true")
  }
  checkmate::assert_flag(quiet)

  # jags_code
  if(!is.null(jags_code))
    if (!is.character(jags_code) || !stringr::str_detect(gsub(" ", "", jags_code), "model\\{"))
      stop("`jags_code` must be NULL or a string with a JAGS model, including 'model {...}'.")


  ##################
  # MODEL BUILDING #
  ##################
  # Build model metadata.
  segment_tables = get_segment_tables(model, data, family, par_x)
  segments = segment_tables$segments
  cps = segment_tables$cps
  assert_model_data(data, par_x, group_cols = stats::na.omit(cps$group_col))
  predictor_tables = get_predictor_tables(model, data, family, par_x)
  predictors = predictor_tables$predictors
  family = resolve_dpar_specs(family, predictors, model)
  group_effects = get_group_effects(cps, predictor_tables$group_effects)

  if (nrow(cps) == 0 && nrow(predictors) == 0 && nrow(group_effects) == 0)
    stop("The model does not contain any parameters to estimate.", call. = FALSE)

  # Make prior
  prior = get_prior(segments, cps, predictors, group_effects, family, prior, data)
  prior_table = attr(prior, "prior_table")
  prior_context = attr(prior, "prior_context")
  attr(prior, "prior_table") = NULL
  attr(prior, "prior_context") = NULL
  assert_fixed_sigma(prior_table, predictors, family)

  # Assemble model metadata used by fitted-model methods
  parameters = get_pars_table(predictors, cps, group_effects, family)
  data_columns = c(
    list(par_x = par_x, response = unique(segments$y), series = series),
    lapply(get_family_aux_columns(family, segments), function(column) {
      if (is.na(column)) NULL else column
    })
  )
  model_tables = list(
    data_columns = data_columns,
    segments = segments,
    cps = cps,
    predictors = predictors,
    group_effects = group_effects,
    parameters = parameters,
    design_specs = predictor_tables$design_specs
  )

  # Validate AR/MA configuration
  has_arma = any(predictors$dpar %in% c("ar", "ma"))
  if (has_arma) {
    if (is.null(family$garma))
      stop(
        "family = ", family$family, "(link = \"", family$link,
        "\") does not define the GARMA behavior required by ar() or ma()."
      )

    response_data = get_family_response_data(family, segments, data)
    assert_arma_boundaries(family, predictors$boundary[predictors$dpar %in% c("ar", "ma")], response_data)

    x_by_series = split(data[[par_x]], if (is.null(series)) 1 else data[[series]])
    x_unordered = any(vapply(
      x_by_series,
      function(x) is.unsorted(x) && is.unsorted(rev(x)),
      logical(1)
    ))
    if (x_unordered)
      message("'", par_x, "' is unordered. Please note that ar() and ma() apply in data-frame row order, not the values of '", par_x, "'.")
  }

  # Make formulas
  formula_jags = get_formula_jags(
    segments, predictors, group_effects, par_x, family,
    design_specs = predictor_tables$design_specs
  )
  formula_r = get_formula_r(formula_jags, predictors, group_effects, cps, par_x)

  # Make jags code if it is not provided by the user
  if (is.null(jags_code)) {
    ar_order = get_arma_order(predictors, "ar")
    ma_order = get_arma_order(predictors, "ma")
    jags_code = get_jags_code(
      prior, segments, group_effects, formula_jags, ar_order, ma_order, family, par_x,
      prior_table, prior_context, series = !is.null(series)
    )
  }


  ##########
  # SAMPLE #
  ##########
  jags_data = get_jags_data(
    data, family, segments, predictors, group_effects, jags_code, series,
    generated = !custom_jags_code,
    design_specs = predictor_tables$design_specs
  )

  # Monitor model parameters and, for generated JAGS code, latent responses
  all_pars = names(prior)
  missing_response_rows = which(is.na(data[[segments$y[1]]]))
  imputed_response_nodes = if (custom_jags_code || length(missing_response_rows) == 0) character() else
    paste0(segments$y[1], "[", missing_response_rows, "]")

  # Sample posterior
  if (sample %in% c("post", "both")) {
    mcmc_post = run_jags(
      jags_code = jags_code,
      jags_data = jags_data,
      pars = c(all_pars, imputed_response_nodes),
      sample = "post",
      n.chains = chains,
      n.iter = iter,
      n.adapt = warmup,
      inits = inits,
      seed = seed,
      quiet = quiet
    ) %>%
      recover_levels(data, group_effects)

    class(mcmc_post) = "mcmc.list"
    assert_ordered_cp_draws(mcmc_post, cps)

    # Handle missing data
    if (length(imputed_response_nodes) > 0) {
      mcmc_imputed = mcmc_post[, imputed_response_nodes, drop = FALSE]
      retained_parameter_nodes = setdiff(colnames(mcmc_post[[1]]), imputed_response_nodes)
      mcmc_post = mcmc_post[, retained_parameter_nodes, drop = FALSE]
    } else {
      mcmc_imputed = NULL
    }


    # Diagnostics check - also for single-chain fits
    fixed_pars = if (!is.null(prior_table$parameter) && !is.null(prior_table$kind))
      prior_table$parameter[prior_table$kind == "constant"] else character()
    warn_nonconvergence(mcmc_post, diagnostics, fixed_pars = fixed_pars)
  } else {
    mcmc_post = NULL
    mcmc_imputed = NULL
  }

  # Sample prior
  if (sample %in% c("prior", "both")) {
    # Set response = NA if we only sample prior
    jags_data_prior = jags_data
    if (!grepl("likelihood_phi_", jags_code, fixed = TRUE))
      jags_data_prior[[segments$y[1]]] = rep(NA, nrow(data))
    if (!is.null(jags_data_prior$response_observed_))
      jags_data_prior$response_observed_[] = 0L

    mcmc_prior = run_jags(
      jags_code = jags_code,
      jags_data = jags_data_prior,
      pars = all_pars,
      sample = "prior",
      n.chains = chains,
      n.iter = iter,
      n.adapt = warmup,
      inits = inits,
      seed = seed,
      quiet = quiet
    ) %>%
      recover_levels(data, group_effects)

    class(mcmc_prior) = "mcmc.list"
    assert_ordered_cp_draws(mcmc_prior, cps)
  } else {
    mcmc_prior = NULL
  }


  ##########
  # RETURN #
  ##########
  # Return normalized formulas without discarding user-defined environments.
  model = Map(stats::as.formula, segments$form, env = segments$form_env)
  class(model) = c("mcplist", "list")
  class(prior) = c("mcplist", "list")
  class(jags_code) = c("mcptext", "character")  # for nicer printing

  # Make mcpfit object
  mcpfit = list(
    # By user (same order as mcp argument)
    model = model,
    data = data,
    prior = prior,
    family = family,
    call = matched_call,

    # Results
    mcmc_post = mcmc_post,
    mcmc_prior = mcmc_prior,

    # Extracted model
    jags_code = jags_code,
    simulate = get_fitsimulate(cps, predictors, group_effects),

    # Pass info to *.mcpfit() functions.
    # Not meant to be used by the end user.
    .internal = list(
      model_tables = model_tables,
      formula_jags = formula_jags,
      formula_r = formula_r,
      prior_table = prior_table,
      prior_context = prior_context,
      prior_format = "jags_string_v1",
      diagnostics = diagnostics,
      custom_jags_code = custom_jags_code,
      imputed_response = mcmc_imputed,
      imputed_response_rows = if (is.null(mcmc_imputed)) integer() else missing_response_rows
    )
  )
  class(mcpfit) = "mcpfit"

  if (!is.null(mcmc_post) && has_arma)
    warn_arma_fit(mcpfit, diagnostics = diagnostics)

  # Return it
  mcpfit
}
