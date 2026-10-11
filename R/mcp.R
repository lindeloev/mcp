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
#' @param model A list of formulas, one for each segment, in the format
#'   `model = list(response ~ predictors, ~ cp(...) + predictors)`, where `cp()` is optional. See [mcp-formula] for how segments connect
#'   and the full formula syntax.
#' @param prior Named list, e.g., `prior = list(cp_1 = "dunif(0, 100)", Intercept_1 = "dnorm(0, 1)")`. 
#'   List names are parameter names (`cp_1`, `Intercept_1`, `x_2`,
#'   `sigma_1`, etc.) and values are a distribution string (e.g., `"dnorm(0, 1) T(0, )"`),
#'   a fixed value (e.g., `-2.1`), or the name of another parameter to share its value.
#'   Weakly informative defaults are used where priors are not specified;
#'   see them with [prior_summary()]. See [mcp-priors] for the full syntax.
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
#'   Defaults to `NULL` where all chains start at equispaced change points and the prior medians.
#'   User-provided `inits` take precedence.
#' @param jags_code String. Pass JAGS code to `mcp` to use directly. This is useful if
#'   you want to tweak the code in `fit$jags_code` and run it within the `mcp`
#'   framework. R-side simulation and prediction methods continue to use the
#'   mcp-default formulas (with warning), so they may no longer match custom JAGS code.
#' @param reparam Logical. Whether to enable reparameterization for faster mixing. Ignored
#'   when `jags_code` is supplied. Read more in
#'   [From formulas to sampler code](https://lindeloev.github.io/mcp/articles/implementation.html#reparameterizations).
#'   * `TRUE` (default) This changes only `fit$jags_code` and sampling efficiency, not the 
#'     model, the priors, or the posterior. 
#'   * `FALSE` runs the plain model. This mostly serve the pedagogical purpose of transparency.
#' @param seed `NULL` or a positive integer. Seed for the JAGS random-number
#'   generators. If `NULL` (default), the seed is drawn from R's random-number generator,
#'   so `set.seed()` before `mcp()` also makes sampling reproducible. 
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
#' @details
#' Here is the demo model, which ships already fitted as `demo_fit` (see also `mcp_example("demo")`):
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
#' Segment 2 continues from where the plateau left off, while segment 3 starts afresh with a new intercept.
#' See [mcp-formula] for the rules of how segments connect, the full formula syntax, and the
#' underlying model. See [mcp-priors] for how to specify priors.
#'
#' @references
#' * Lindeløv, J. K. (2020). mcp: An R Package for Regression With Multiple Change Points.
#'   *OSF Preprints*. \doi{10.31219/osf.io/fzqxv}
#'   Introduces the `mcp` package, formula syntax, default priors, and workflow for regression
#'   with multiple change points across generalized linear and time-series models. Newer-than-2020 versions of
#'   the paper may be available at that link. Please cite the newest version.
#' * Carlin, B. P., Gelfand, A. E., & Smith, A. F. (1992). Hierarchical Bayesian Analysis of
#'   Changepoint Problems. *Applied Statistics*, 41(2), 389–405. \doi{10.2307/2347570}
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
#'   cp_2 = "dunif(cp_1, 80)"        # change point to segment 3 > cp_1 and < 80.
#' )
#'
#' fit3 = mcp(model, data = data, prior = prior, seed = 42)
#'
#' # Share coefficients across segments using same() (e.g., reuse Intercept_1 in segment 3)
#' model_same = list(
#'   response ~ 1,
#'   ~ 0 + time,
#'   ~ same(1, as = 1) + time  # Reuse Intercept_1 instead of estimating Intercept_3
#' )
#' fit_same = mcp(model_same, data = data, sample = FALSE)
#'
#' # Show the JAGS model
#' demo_fit$jags_code
#' }
mcp = function(
  model,
  data,
  prior = list(),
  family = gaussian(),
  par_x = NULL,
  sample = "post",
  chains = 3,
  iter = 3000,
  warmup = 1000,
  inits = NULL,
  jags_code = NULL,
  reparam = TRUE,
  seed = NULL,
  diagnostics = list(),
  quiet = FALSE,
  cores = NULL,
  adapt = lifecycle::deprecated()
) {
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
  assert_reserved_output_namespace(names(data), context = "data")
  assert_reserved_offset_namespace(names(data), context = "data")
  non_syntactic = names(data)[make.names(names(data)) != names(data)]
  if (length(non_syntactic) > 0)
    stop(
      "`data` has non-syntactic column name(s): ",
      paste0("`", non_syntactic, "`", collapse = ", "), ". Rename them before fitting.",
      call. = FALSE
    )
  data = data.frame(data)
  vars = get_formula_vars(model)
  series = logical0_to_null(vars$series)
  assert_arma_series(data, series)

  checkmate::assert_string(par_x, null.ok = TRUE)
  par_x = get_par_x(vars, data, par_x)
  model_vars = unique(c(unlist(lapply(model, all.vars)), par_x, series))  # In formula order
  assert_model_data(data, par_x, unlist(vars[c("predictor", "group", "offset", "series")]))
  assert_data_cols(data, cols = model_vars, fail_funcs = c(is.infinite))
  data = data[, model_vars]  # Remove unused data
  for (col in vars$response) {
    if (is.logical(data[[col]]))
      data[[col]] = as.numeric(data[[col]])
    if (!is.numeric(data[[col]]))
      stop("The response column '", col, "' must be numeric or logical. Got ", class(data[[col]])[1], ".", call. = FALSE)
  }

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
  checkmate::assert_flag(reparam)
  checkmate::assert_flag(quiet)

  # jags_code
  if (custom_jags_code && (!is.character(jags_code) || !stringr::str_detect(gsub(" ", "", jags_code), "model\\{")))
    stop("`jags_code` must be NULL or a string with a JAGS model, containing 'model {...}'.")


  ##################
  # MODEL BUILDING #
  ##################
  # Build model metadata.
  segment_tables = get_segment_tables(model, data, family, par_x)
  segments = segment_tables$segments
  cps = segment_tables$cps
  # Pass the same canonical segment formulas to all model-building consumers
  model = Map(stats::as.formula, segments$form, env = segments$form_env)
  segments$form_env = NULL
  assert_model_data(data, par_x, group_cols = stats::na.omit(cps$group_col))
  predictor_tables = get_predictor_tables(model, data, family, par_x)
  predictor_occurrences = predictor_tables$predictors
  group_occurrences = get_group_effects(cps, predictor_tables$group_effects)
  predictor_definitions = definition_rows(predictor_occurrences)
  group_definitions = definition_rows(group_occurrences)
  family = resolve_dpar_specs(family, predictor_tables$components)

  if (nrow(cps) == 0 && nrow(predictor_definitions) == 0 && nrow(group_definitions) == 0)
    stop("The model does not contain any parameters to estimate.", call. = FALSE)

  # Check that priors are set on the original parameter, not ones replaced by same()
  shared = dplyr::bind_rows(
    dplyr::transmute(
      dplyr::filter(predictor_occurrences, .data$segment != .data$definition_segment),
      replaced = stringr::str_replace(.data$code_name, "_[0-9]+$", paste0("_", .data$segment)),
      original = .data$code_name
    ),
    dplyr::transmute(
      dplyr::filter(group_occurrences, .data$segment != .data$definition_segment),
      replaced = stringr::str_replace(.data$sd_name, paste0("_", .data$definition_segment, "_"), paste0("_", .data$segment, "_")),
      original = .data$sd_name
    )
  ) %>% dplyr::distinct()
  bad_prior = shared[shared$replaced %in% names(prior), , drop = FALSE]
  if (nrow(bad_prior) > 0) {
    guidance = if (nrow(bad_prior) == 1) {
      paste0("`", bad_prior$original, "` instead of `", bad_prior$replaced, "`")
    } else {
      paste0("`", bad_prior$original, "` (not `", bad_prior$replaced, "`)")
    }
    stop("`same()` reuses earlier parameters instead of creating new ones. Set prior on ",
      and_collapse(guidance), ".", call. = FALSE)
  }

  # Make prior
  prior = get_prior(segments, cps, predictor_definitions, group_definitions, family, prior, data, predictor_tables$design_specs, reparam)
  prior_table = attr(prior, "prior_table")
  attr(prior, "prior_table") = NULL
  assert_fixed_dpars(prior_table, predictor_definitions, family)

  # Assemble model metadata used by fitted-model methods
  parameters = get_pars_table(predictor_definitions, cps, group_definitions, family)
  data_columns = c(
    list(par_x = par_x, response = unique(segments$y), series = series),
    lapply(get_family_aux_columns(family, segments), function(column) {
      if (is.na(column)) NULL else column
    })
  )
  # Data-sized design values are only needed to build the model. Fitted-model
  # methods evaluate the stored design specifications on (new)data instead.
  model_tables = list(
    data_columns = data_columns,
    segments = segments,
    cps = cps,
    predictors = dplyr::select(predictor_occurrences, -"matrix_data"),
    group_effects = dplyr::select(group_occurrences, -"matrix_data"),
    parameters = parameters,
    design_specs = lapply(predictor_tables$design_specs, function(spec) {
      spec$offset_data = NULL
      spec
    })
  )

  # Validate AR/MA configuration
  has_arma = any(predictor_definitions$dpar %in% c("ar", "ma"))
  if (has_arma) {
    if (is.null(family$garma))
      stop(
        "family = ", family$family, "(link = \"", family$link,
        "\") does not define the GARMA behavior required by ar() or ma()."
      )

    response_data = get_family_response_data(family, segments, data)
    assert_arma_thresholds(family, predictor_definitions$threshold[predictor_definitions$dpar %in% c("ar", "ma")], response_data)

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
    segments, predictor_occurrences, group_occurrences, par_x, family,
    design_specs = predictor_tables$design_specs
  )

  # Make jags code if it is not provided by the user
  if (!custom_jags_code) {
    ar_order = get_arma_order(predictor_definitions, "ar")
    ma_order = get_arma_order(predictor_definitions, "ma")
    jags_code = get_jags_code(
      prior, segments, group_definitions, formula_jags, ar_order, ma_order, family, par_x,
      prior_table, series = !is.null(series)
    )
  }


  ##########
  # SAMPLE #
  ##########
  jags_data = get_jags_data(
    data, family, segments, predictor_occurrences, group_occurrences, jags_code, series,
    generated = !custom_jags_code,
    design_specs = predictor_tables$design_specs
  )

  # Monitor model parameters and, for generated JAGS code, latent responses
  all_pars = names(prior)
  missing_response_rows = which(is.na(data[[segments$y[1]]]))
  imputed_response_nodes = if (custom_jags_code || length(missing_response_rows) == 0) character() else
    paste0(segments$y[1], "[", missing_response_rows, "]")

  # The JAGS glm module helps some families always and others only with group-level predictor effects.
  glm_module = family$backends$jags$glm_module
  use_glm = identical(glm_module, "always") ||
    (identical(glm_module, "group_effects") && any(group_definitions$part != "cp"))

  sample_draws = function(jags_data, pars, sample) {
    draws = run_jags(
      jags_code = jags_code,
      jags_data = jags_data,
      pars = pars,
      sample = sample,
      n.chains = chains,
      n.iter = iter,
      n.adapt = warmup,
      inits = inits,
      seed = seed,
      quiet = quiet,
      use_glm = use_glm
    ) %>%
      recover_levels(data, group_definitions)

    class(draws) = "mcmc.list"
    assert_ordered_cp_draws(draws, cps)
    draws
  }

  # Sample posterior
  if (sample %in% c("post", "both")) {
    mcmc_post = sample_draws(jags_data, c(all_pars, imputed_response_nodes), "post")

    # Handle missing data
    if (length(imputed_response_nodes) > 0) {
      mcmc_imputed = mcmc_post[, imputed_response_nodes, drop = FALSE]
      retained_parameter_nodes = setdiff(colnames(mcmc_post[[1]]), imputed_response_nodes)
      mcmc_post = mcmc_post[, retained_parameter_nodes, drop = FALSE]
    } else {
      mcmc_imputed = NULL
    }


    # Diagnostics check - also for single-chain fits
    fixed_pars = prior_table$parameter[prior_table$kind == "constant"]
    warn_nonconvergence(mcmc_post, diagnostics, fixed_pars = fixed_pars)
  } else {
    mcmc_post = NULL
    mcmc_imputed = NULL
  }

  # Sample prior
  if (sample %in% c("prior", "both")) {
    # Set stochastic response nodes to NA for prior sampling. 
    # Except for weighted discrete likelihoods which use the response in a zeros trick.
    jags_data_prior = jags_data
    response_node = paste0(segments$y[1], "[i_] ~")
    response_is_stochastic = if (custom_jags_code)
      !grepl("likelihood_phi_", jags_code, fixed = TRUE) else
      grepl(response_node, jags_code, fixed = TRUE)
    if (response_is_stochastic)
      jags_data_prior[[segments$y[1]]] = rep(NA, nrow(data))
    if (!is.null(jags_data_prior$response_observed_))
      jags_data_prior$response_observed_[] = 0L

    mcmc_prior = sample_draws(jags_data_prior, all_pars, "prior")
  } else {
    mcmc_prior = NULL
  }


  ##########
  # RETURN #
  ##########
  # Return canonical formulas attached to .GlobalEnv to avoid environment bloat
  model = lapply(segments$form, stats::as.formula, env = globalenv())
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
    simulate = get_fitsimulate(cps, predictor_definitions, group_definitions),

    # Pass info to *.mcpfit() functions.
    # Not meant to be used by the end user.
    .internal = list(
      model_tables = model_tables,
      formula_jags = formula_jags,
      prior_table = prior_table,
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
