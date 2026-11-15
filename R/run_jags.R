# ABOUT: These functions call the JAGS sampler and set up inputs in a way
# that is appropriate for that, including data (get_jags_data()).
# ------------------------------------------------------------------------


# Sample the posterior and/or prior using JAGS
# Holds all JAGS-specific sampling steps for `mcp()`.
sample_jags = function(model_tables, prior_table, family, data, sample, chains, iter, warmup, inits, jags_code, custom_jags_code, seed, quiet) {
  data_columns = model_tables$data_columns
  parameters = model_tables$parameters
  jags_data = get_jags_data(data, family, model_tables, jags_code, generated = !custom_jags_code)

  # Monitor model parameters and, for generated JAGS code, latent responses
  all_pars = prior_table$parameter
  missing_response_rows = which(is.na(data[[data_columns$response]]))
  imputed_response_nodes = if (custom_jags_code || length(missing_response_rows) == 0) character() else
    paste0(data_columns$response, "[", missing_response_rows, "]")

  sample_draws = function(pars, sample) {
    draws = run_jags(
      jags_code = jags_code,
      jags_data = if (sample == "prior") get_jags_prior_data(jags_data, jags_code, data_columns) else jags_data,
      pars = pars,
      sample = sample,
      n.chains = chains,
      n.iter = iter,
      n.adapt = warmup,
      inits = inits,
      default_inits = get_default_inits(prior_table, parameters, range(data[[data_columns$par_x]], na.rm = TRUE), inits, if (custom_jags_code) jags_code,
                                        cp_deviance = if (!custom_jags_code && sample == "post") get_cp_deviance(model_tables, family, data)),
      seed = seed,
      quiet = quiet,
      use_glm = use_jags_glm(family, definition_rows(model_tables$group_coefs))
    ) |>
      recover_levels(data, parameters[parameters$scope == "group", ])

    class(draws) = "mcmc.list"

    # JAGS sorts nodes alphabetically. Order them like mcp_pars() instead, keeping group levels in factor-level order.
    draws[, order(match(sub("\\[.*", "", coda::varnames(draws)), parameters$parameter)), drop = FALSE]
  }

  # Sample posterior and split off the missing responses
  mcmc_post = NULL
  imputed_response = NULL
  if (sample %in% c("post", "both")) {
    mcmc_post = sample_draws(c(all_pars, imputed_response_nodes), "post")
    if (length(imputed_response_nodes) > 0) {
      imputed_response = mcmc_post[, imputed_response_nodes, drop = FALSE]
      mcmc_post = mcmc_post[, setdiff(colnames(mcmc_post[[1]]), imputed_response_nodes), drop = FALSE]
    }
  }

  # Sample prior
  mcmc_prior = if (sample %in% c("prior", "both")) sample_draws(all_pars, "prior")

  list(mcmc_post = mcmc_post, mcmc_prior = mcmc_prior, imputed_response = imputed_response)
}


#' Run parallel MCMC sampling using JAGS.
#'
#' @inheritParams mcp
#' @inheritParams rjags::jags.model
#' @inheritParams rjags::coda.samples
#' @param jags_code A string. JAGS model, usually returned by `get_jags_code()`.
#' @param pars Character vector of parameters to save/monitor.
#' @param default_inits Named list. Initial values for parameters not in `inits`, e.g., from `get_cp_inits()`.
#' @param use_glm Logical. Load the JAGS glm module, which block-samples
#'   coefficients that enter the linear predictor jointly.
#' @return `mcmc.list`
#' @noRd
run_jags = function(jags_code,
                    jags_data,
                    pars,
                    sample,
                    n.chains,
                    n.iter,
                    n.adapt,
                    inits,
                    default_inits = list(),
                    seed,
                    quiet,
                    use_glm = FALSE
) {

  # Define the sampling function in this environment.
  # Can be used sequentially or in parallel.
  do_sampling = function(inits, n.chains) {
    sample_chains = function() {
      # Samplers are assigned at compile time, so load the glm module here to reach
      # parallel workers too, and unload it again to leave the user's session unchanged.
      if (use_glm && !("glm" %in% rjags::list.modules())) {
        rjags::load.module("glm", quiet = TRUE)
        on.exit(rjags::unload.module("glm", quiet = TRUE), add = TRUE)
      }

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
  if (is.null(seed))
    seed = sample.int(.Machine$integer.max, 1)
  inits = get_jags_inits(inits, seed, n.chains, sample)
  inits = lapply(inits, function(chain_inits) c(chain_inits, default_inits[setdiff(names(default_inits), names(chain_inits))]))

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

  # Zeros trick (see context$zeros_trick() in get_jags_family_context())
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

  # Observed nodes that restrict the change point priors, also in custom code copied from fit$jags_code
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
