# ABOUT: These functions call the JAGS sampler and set up inputs in a way
# that is appropriate for that, including data (get_jags_data()).
# ------------------------------------------------------------------------


# Sample the posterior and/or prior using JAGS
# Holds all JAGS-specific sampling steps for `mcp()`.
sample_jags = function(model_tables, prior_table, family, data, sample, chains, iter, warmup, inits, jags_code, custom_jags_code, seed, quiet) {
  data_columns = model_tables$data_columns
  parameters = model_tables$parameters
  jags_data = get_jags_data(data, family, model_tables, jags_code, generated = !custom_jags_code)

  # Monitor model parameters and, for generated JAGS code, latent responses and the deviance (see check_jags_likelihood())
  all_pars = prior_table$parameter
  missing_response_rows = which(is.na(data[[data_columns$response]]))
  imputed_response_nodes = if (custom_jags_code || length(missing_response_rows) == 0) character() else
    paste0(data_columns$response, "[", missing_response_rows, "]")
  deviance_node = if (!custom_jags_code) "deviance"

  sample_draws = function(pars, sample) {
    draws = run_jags(
      jags_code = jags_code,
      jags_data = if (sample == "prior") get_jags_prior_data(jags_data, jags_code, data_columns) else jags_data,
      pars = pars,
      n.chains = chains,
      n.iter = iter,
      n.adapt = warmup,
      inits = get_chain_inits(model_tables, prior_table, family, data, sample, chains, inits, seed, if (custom_jags_code) jags_code),
      quiet = quiet,
      use_glm = use_jags_glm(family, definition_rows(model_tables$group_coefs))
    ) |>
      recover_levels(data, parameters[parameters$scope == "group", ])

    class(draws) = "mcmc.list"

    # JAGS sorts nodes alphabetically. Order them like mcp_pars() instead, keeping group levels in factor-level order.
    draws[, order(match(sub("\\[.*", "", coda::varnames(draws)), parameters$parameter)), drop = FALSE]
  }

  # Sample posterior and split off the missing responses and the deviance
  mcmc_post = NULL
  imputed_response = NULL
  deviance = NULL
  if (sample %in% c("post", "both")) {
    mcmc_post = sample_draws(c(all_pars, imputed_response_nodes, deviance_node), "post")
    if (length(imputed_response_nodes) > 0)
      imputed_response = mcmc_post[, imputed_response_nodes, drop = FALSE]
    if (!is.null(deviance_node))
      deviance = mcmc_post[, deviance_node, drop = FALSE]
    mcmc_post = mcmc_post[, setdiff(colnames(mcmc_post[[1]]), c(imputed_response_nodes, deviance_node)), drop = FALSE]
  }

  # Sample prior
  mcmc_prior = if (sample %in% c("prior", "both")) sample_draws(all_pars, "prior")

  list(mcmc_post = mcmc_post, mcmc_prior = mcmc_prior, imputed_response = imputed_response, deviance = deviance)
}


# Initial values for each chain: the user's `inits` with JAGS RNG seeds (see get_jags_inits()) and defaults for the
# other parameters (see get_default_inits()). Used for sampling and by as_jags().
# - jags_code: Custom JAGS code or NULL for the code that mcp writes.
get_chain_inits = function(model_tables, prior_table, family, data, sample, chains, inits, seed, jags_code = NULL) {
  x_range = range(data[[model_tables$data_columns$par_x]], na.rm = TRUE)
  cp_deviance = if (is.null(jags_code) && sample == "post") get_cp_deviance(model_tables, family, data)
  defaults = get_default_inits(prior_table, model_tables$parameters, x_range, inits, jags_code, cp_deviance)
  lapply(get_jags_inits(inits, seed, chains, sample), function(chain_inits) {
    c(chain_inits, defaults[setdiff(names(defaults), names(chain_inits))])
  })
}


#' Run parallel MCMC sampling using JAGS.
#'
#' @inheritParams mcp
#' @inheritParams rjags::jags.model
#' @inheritParams rjags::coda.samples
#' @param jags_code A string. JAGS model, usually returned by `get_jags_code()`.
#' @param pars Character vector of parameters to save/monitor. `"deviance"` loads the JAGS dic module.
#' @param inits List with initial values for each chain, e.g., from `get_chain_inits()`.
#' @param use_glm Logical. Load the JAGS glm module, which block-samples
#'   coefficients that enter the linear predictor jointly.
#' @return `mcmc.list`
#' @noRd
run_jags = function(jags_code,
                    jags_data,
                    pars,
                    n.chains,
                    n.iter,
                    n.adapt,
                    inits,
                    quiet,
                    use_glm = FALSE
) {

  # Define the sampling function in this environment.
  # Can be used sequentially or in parallel.
  do_sampling = function(inits, n.chains) {
    sample_chains = function() {
      # Samplers are assigned at compile time, so load the modules here to reach
      # parallel workers too, and unload them again to leave the user's session unchanged.
      modules = setdiff(c(if (use_glm) "glm", if ("deviance" %in% pars) "dic"), rjags::list.modules())
      for (module in modules)
        rjags::load.module(module, quiet = TRUE)
      on.exit(for (module in modules) rjags::unload.module(module, quiet = TRUE), add = TRUE)

      # Compile model silently. Adapt separately to keep the progress bar.
      jags_connection = textConnection(jags_code)
      on.exit(close(jags_connection), add = TRUE)
      jags_model = rjags::jags.model(
        file = jags_connection,
        data = jags_data,
        inits = inits,
        n.chains = n.chains,
        n.adapt = 0,
        quiet = TRUE
      )
      progress_bar = if (quiet) "none" else getOption("jags.pb")
      if (!rjags::adapt(jags_model, n.iter = n.adapt, end.adaptation = TRUE, progress.bar = progress_bar))
        warning("Adaptation incomplete")

      # Ensure requested warmup iterations are discarded even if JAGS skips adaptation
      burnin = n.adapt - jags_model$iter()
      if (burnin > 0) {
        stats::update(
          object = jags_model,
          n.iter = burnin,
          progress.bar = progress_bar
        )
      }

      # Sample and return
      rjags::coda.samples(
        model = jags_model,
        variable.names = pars,
        n.iter = n.iter,
        quiet = quiet
      )
    }

    # rjags still prints some adaptation notices with quiet = TRUE.
    if (quiet) {
      draws = NULL
      capture.output({ draws = sample_chains() })
    } else {
      draws = sample_chains()
    }
    draws
  }

  # Use JAGS directly under a sequential future plan. This compiles the model
  # only once for all chains and preserves JAGS's progress output.
  n_workers = future::nbrOfWorkers()
  timer = proc.time()
  if (n_workers == 1) {
    draws = do_sampling(
      inits = inits,
      n.chains = n.chains
    )
  } else {
    # Submit one chain per future. The user's future plan controls the backend.
    if (!quiet)
      message("Parallel sampling in progress...")
    draws = future.apply::future_lapply(
      inits,
      n.chains = 1,
      FUN = do_sampling,
      future.seed = TRUE
    )

    # Get result as mcmc.list
    draws = unlist(draws, recursive = FALSE)
    class(draws) = "mcmc.list"
  }

  # Sampling finished
  passed = proc.time() - timer
  if (!quiet)
    message("Finished sampling in ", round(passed["elapsed"], 1), " seconds\n")

  draws
}


# The JAGS glm module helps some families always and others only with group-level effects in the predictor part.
use_jags_glm = function(family, group_coefs) {
  glm_module = family$backends$jags$glm_module
  identical(glm_module, "always") || (identical(glm_module, "group_effects") && nrow(group_coefs) > 0)
}


# JAGS data for sampling the prior. The response is set to NA, except for likelihoods without a
# response node, e.g., the zeros trick, which `response_observed_ = 0` turns off.
get_jags_prior_data = function(jags_data, jags_code, data_columns) {
  if (has_response_node(jags_code, data_columns))
    jags_data[[data_columns$response]][] = NA
  if (!is.null(jags_data$response_observed_))
    jags_data$response_observed_[] = 0L
  jags_data
}


# Add deterministic JAGS RNG initial values to the `inits` passed to `mcp()`.
# Returns one list of initial values per chain. Prior sampling (`sample = "prior"`) uses other chain seeds.
get_jags_inits = function(inits, seed, n.chains, sample) {
  if (is.list(inits[[1]])) {
    if (length(inits) != n.chains)
      stop("`inits` must have length equal to `chains` (", n.chains, ") when supplying chain-specific initialization lists.")
  } else {
    inits = rep(list(inits), n.chains)
  }

  # Space the chain seeds by the golden-ratio fraction of the seed range, so that nearby seeds,
  # e.g., 1 and 2, never share chains (for seeds less than 10^7 apart).
  chain = seq_len(n.chains) - 1 + if (sample == "prior") n.chains else 0
  rng_seed = as.integer((as.double(seed) - 1 + chain * 1327217884) %% .Machine$integer.max + 1)

  lapply(seq_len(n.chains), function(i) {
    chain_inits = inits[[i]]
    chain_inits[c(".RNG.name", ".RNG.seed", ".RNG.state")] = NULL
    c(chain_inits, list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = rng_seed[i]))
  })
}


# Check names of user-data nodes passed to JAGS against generated JAGS nodes.
# - generated: Whether `jags_code` was generated by mcp (`FALSE` for custom `jags_code`).
# - parameters: The parameters table from `get_pars_table()`.
assert_jags_namespace = function(data_names, family, n_segments, parameters, jags_code, generated = TRUE) {
  group_cols = unique(stats::na.omit(parameters$group_col))
  dpars = family$dpar_specs$dpar

  # Parameters and the ends of the observed range (cp_0 and cp_K) share one namespace.
  parameter_names = c(paste0("cp_", c(0, n_segments)), parameters$parameter)

  # Data helpers are always present; model helpers occur only in generated code.
  offset_nodes = if (!is.null(jags_code)) unique(stringr::str_extract_all(jags_code, "offset_[A-Za-z0-9_]+_")[[1]])
  helper_names = unique(c(
    "rhs_matrix_", "series_id_", "response_observed_",
    "likelihood_weight_", "likelihood_phi_", "likelihood_zero_", offset_nodes,
    paste0("n_unique_", group_cols), paste0(group_cols, "_"), paste0(group_cols, "_index_"),
    if (generated) c(
      "i_", "garma_threshold_", "garma_y_", "garma_link_y_",
      "resid_abs_", "resid_ma_", "resid_garma_", "nb_prob_",
      "cp_order_", "cp_weight_",
      paste0("x_local_", seq_len(n_segments), "_"),
      paste0("link_", dpars, "_"), paste0(dpars, "_"),
      paste0(parameters$parameter, "_inverse")
    )
  ))

  # Catch parameter-to-parameter and parameter-to-helper collisions before priors.
  duplicate_parameters = unique(parameter_names[duplicated(parameter_names)])
  parameter_helpers = intersect(parameter_names, helper_names)
  generated_collisions = unique(c(duplicate_parameters, parameter_helpers))
  if (length(generated_collisions) > 0)
    stop(
      "Generated parameter name(s) collide in the JAGS namespace: ",
      and_collapse(paste0("'", generated_collisions, "'")),
      ". Rename the corresponding predictor column."
    )

  # Finally, ensure user data cannot overwrite any generated JAGS node.
  collisions = intersect(data_names, c(parameter_names, helper_names))
  reserved_prefix = grep("^\\.?mcp_", data_names, value = TRUE)
  bad = unique(c(collisions, reserved_prefix))
  if (length(bad) > 0)
    stop(
      "Data column name(s) collide with mcp's generated JAGS namespace: ",
      and_collapse(paste0("'", bad, "'")), ". Rename them before fitting."
    )
  invisible(TRUE)
}


#' Adds helper variables for use in `run_jags`
#'
#' Returns the relevant data columns as a list and add elements with unique
#' grouping-factor levels.
#'
#' @inheritParams run_jags
#' @param data A tibble
#' @param model_tables Model tables, as in `get_fit_model_tables(fit)`.
#' @noRd
get_jags_data = function(data, family, model_tables, jags_code, generated = TRUE) {
  data_columns = model_tables$data_columns
  parameters = model_tables$parameters
  n_segments = nrow(model_tables$cps) + 1
  group_cols = unique(stats::na.omit(parameters$group_col))
  series = data_columns$series

  # Start with "raw" data
  aux_columns = get_family_aux_columns(family, data_columns)
  cols_data = unique(stats::na.omit(c(data_columns$response, data_columns$par_x, unname(aux_columns))))
  assert_jags_namespace(
    unique(c(names(data), group_cols, cols_data, series)), family, n_segments, parameters, jags_code, generated
  )
  jags_data = as.list(data[, cols_data])

  if (!is.null(series))
    jags_data$series_id_ = as.integer(factor(data[[series]], levels = unique(data[[series]])))

  for (group_col in group_cols) {
    # Add metadata for the grouping-factor levels.
    n_unique_name = paste0("n_unique_", group_col)
    group_levels = levels(factor(data[[group_col]]))
    jags_data[[n_unique_name]] = length(group_levels)

    # Index grouping levels in a separate node so the grouping column can also
    # be par_x, the response, etc. Levels are recovered using recover_levels()
    jags_data[[paste0(group_col, "_index_")]] = as.numeric(factor(data[[group_col]], levels = group_levels))
  }

  # Predictor design matrix and offsets, evaluated like predictions. Keep the JAGS data name for custom-code compatibility.
  designs = evaluate_designs(model_tables, data)
  jags_data$rhs_matrix_ = designs$matrix
  jags_data[names(designs$offsets)] = designs$offsets

  # Compatibility for custom JAGS code written with the released data
  # constants. Generated code and newly resolved priors no longer need these.
  context = prior_context(data, data_columns, n_segments, model_tables$design_specs)
  jags_data = legacy_prior_jags_data(jags_data, jags_code, context)

  # Zeros trick (see context$zeros_trick() in get_jags_family_context()), also in custom code copied from as_jags(fit)$code
  if (grepl("likelihood_zero_[i_] ~", jags_code, fixed = TRUE)) {
    if (anyNA(data[[data_columns$response]]) && !has_response_node(jags_code, data_columns))
      stop(
        "The likelihood of family = ", family$family, "() uses the response as data, e.g., to apply weights, ",
        "so it cannot impute missing responses (NA in '", data_columns$response, "'). Remove these rows before fitting.",
        call. = FALSE
      )
    jags_data$response_observed_ = as.integer(!is.na(data[[data_columns$response]]))
    jags_data$likelihood_zero_ = numeric(nrow(data))
  }

  # Observed nodes that restrict the change point priors, also in custom code copied from as_jags(fit)$code
  if (grepl("cp_order_ ~", jags_code, fixed = TRUE))
    jags_data$cp_order_ = 1
  if (grepl("cp_weight_ ~", jags_code, fixed = TRUE))
    jags_data$cp_weight_ = 0

  # Return
  jags_data
}


# Whether the likelihood has a node for the response, e.g., `y[i_] ~ dnorm(...)`. Otherwise, e.g., with the zeros trick,
# the response is data, so JAGS can neither impute missing responses nor sample it from the prior.
has_response_node = function(jags_code, data_columns) {
  grepl(paste0(" ", data_columns$response, "[i_] ~"), jags_code, fixed = TRUE)
}


# Recover the original levels (numeric or string) of group-level effects in an mcmc.list.
# JAGS uses 1, 2, 3, ... for indexing group-level effects.
# - group_pars: Group-level rows of the parameters table, i.e., `mcp_pars(fit, scope = "group")`.
recover_levels = function(draws, data, group_pars) {
  # JAGS names levels by index, e.g., "x_1_id[2]"; use the original levels, e.g., "x_1_id[Omar]"
  index_names = character()
  level_names = character()
  for (i in seq_len(nrow(group_pars))) {
    group_levels = levels(factor(data[[group_pars$group_col[i]]]))
    index_names = c(index_names, paste0(group_pars$parameter[i], "[", seq_along(group_levels), "]"))
    level_names = c(level_names, paste0(group_pars$parameter[i], "[", group_levels, "]"))
  }

  lapply(draws, function(x) {
    matched = match(colnames(x), index_names)
    colnames(x)[!is.na(matched)] = level_names[matched[!is.na(matched)]]
    x
  })
}


# How much the JAGS likelihood and the R-side log_lik() disagree in `ndraws` posterior draws of `fit`. 
# Should ideally be zero, i.e., mcp's JAGS code, log_lik(), and loo() describe the same model.
# The JAGS deviance may differ from -2 * sum(log_lik()) by a constant (e.g., from the zeros trick in
# context$zeros_trick()), so the check is on their difference across draws. NA where log_lik() does not apply.
# - deviance: mcmc.list with the JAGS deviance of each draw in `fit`.
get_likelihood_mismatch = function(fit, deviance, ndraws = 20) {
  columns = mcp_columns(fit)
  if (is_arma(fit) && anyNA(fit$data[[columns$response]]))
    return(NA_real_)  # log_lik() refuses missing AR/MA histories

  ndraws = min(ndraws, count_draws(fit))
  fit = thin_post_draws(fit, ndraws)
  jags_deviance = as.numeric(as.matrix(thin_chains(deviance, ndraws))) - get_cp_weight_deviance(fit)
  r_deviance = -2 * rowSums(log_lik(fit, summary = FALSE, draws_format = "matrix"), na.rm = TRUE)
  difference = jags_deviance - r_deviance
  max(abs(difference - stats::median(difference))) / (1 + mean(abs(r_deviance)))
}


# Warn if mcp's JAGS code and R-side methods describe different likelihoods (see get_likelihood_mismatch()),
# e.g., when `jags$likelihood` and `r$log_lik` of a custom family disagree.
check_jags_likelihood = function(fit, deviance) {
  if (isTRUE(get_likelihood_mismatch(fit, deviance) > 1e-6))
    warning(
      "The JAGS likelihood and the R-side `log_lik()` of family = ", fit$family$family, "() disagree, so R-side methods ",
      "such as `loo()`, `predict()`, and `plot()` may be wrong. For families made with `mcpfamily()`, check that ",
      "`jags$likelihood` and `r$log_lik` describe the same distribution.", call. = FALSE
    )
}


# The JAGS deviance of the dirichlet weight on the segment lengths (`cp_weight_` in get_jags_code()) in each draw of `fit`
get_cp_weight_deviance = function(fit) {
  model_tables = get_fit_model_tables(fit)
  alpha = get_dirichlet_alpha(fit$.internal$prior_table, model_tables$cps)
  if (alpha == 1)
    return(0)

  # cp_weight_ ~ dpois(-(alpha - 1) * sum(log(segment length / range of x))) adds 2 * lambda to the deviance
  x_range = range(fit$data[[model_tables$data_columns$par_x]])
  cps = cbind(x_range[1], as.matrix(mcmclist_draws(fit))[, model_tables$cps$name, drop = FALSE], x_range[2])
  -2 * (alpha - 1) * rowSums(log((cps[, -1, drop = FALSE] - cps[, -ncol(cps), drop = FALSE]) / diff(x_range)))
}
