# ABOUT: These functions call the JAGS sampler and set up inputs in a way
# that is appropriate for that, including data (get_jags_data()).
# ------------------------------------------------------------------------


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
    sample_jags = function() {
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
      capture.output({ draws = sample_jags() })
    } else {
      draws = sample_jags()
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


# JAGS data for sampling the prior. The response is set to NA, except for weighted discrete
# likelihoods, which use the response in a zeros trick that `response_observed_ = 0` turns off.
get_jags_prior_data = function(jags_data, family, data_columns) {
  if (is.null(data_columns$weights) || !isTRUE(family$response$is_discrete))
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


# Check that `inits` only set parameters that are sampled directly (see add_nodes())
assert_inits = function(inits, pars) {
  invalid = setdiff(unlist(lapply(if (is.list(inits[[1]])) inits else list(inits), names)), pars)
  if (length(invalid) > 0)
    stop("`inits` can only set these parameters of this model: `", paste0(pars, collapse = "`, `"), "`. Got `", paste0(unique(invalid), collapse = "`, `"),
         "`. Parameters with constant priors and parameters that are reparameterized for faster sampling (see `mcp_reparam` in ?mcp-package) cannot be initialized.", call. = FALSE)
}


# Initial values where users supplied none: ordered change points (see get_cp_inits()) and prior medians of
# other population-level parameters, e.g., SDs away from zero. Custom JAGS code may define parameters differently.
get_default_inits = function(prior_table, parameters, x_range, inits, custom_jags_code = NULL) {
  is_cp = parameters$scope == "population" & parameters$role == "change_point"
  cp_rows = prior_table[match(parameters$parameter[is_cp], prior_table$parameter), ]
  is_median = !is.na(prior_table$node) & prior_table$parameter %in% parameters$parameter[parameters$scope == "population" & !is_cp]
  node_priors = stats::setNames(prior_table$node_prior[is_median], prior_table$node[is_median])
  defaults = c(get_cp_inits(cp_rows, x_range, if (!is.list(inits[[1]])) inits), get_median_inits(node_priors))  # Ignores chain-specific inits
  if (is.null(custom_jags_code)) defaults else
    defaults[vapply(names(defaults), function(name) grepl(paste0("\\n\\s*", name, " ~"), custom_jags_code), logical(1))]
}


# Initial values at the prior medians of nodes. Covers priors whose arguments and truncation bounds are numbers.
# - prior: Named list of parsed priors, e.g., `prior_table$node_prior` named by `prior_table$node`.
get_median_inits = function(prior) {
  as_number = function(x) tryCatch(as.numeric(eval(str2lang(x), baseenv())), error = function(e) NA_real_)

  inits = list()
  for (name in names(prior)) {
    parsed = prior[[name]]
    args = vapply(parsed$args, as_number, numeric(1), USE.NAMES = FALSE)
    bounds = c(-Inf, Inf)
    is_bounded = nzchar(parsed$trunc)
    bounds[is_bounded] = vapply(parsed$trunc[is_bounded], as_number, numeric(1), USE.NAMES = FALSE)
    if (!anyNA(args) && !anyNA(bounds))
      inits[[name]] = prior_quantile(parsed$dist, args, bounds, 0.5)
  }
  inits
}


# Initial values for the change points with a prior distribution, except those in `user_inits`:
# equispaced in the observed range (`x_range`) with room for minimum segment lengths, then moved the
# least possible to respect their bounds and order. JAGS fails to initialize if they violate cp_order_.
# - cp_rows: Rows of the prior table for the population-level change points, in order.
get_cp_inits = function(cp_rows, x_range, user_inits = list()) {
  n_cp = nrow(cp_rows)
  if (n_cp == 0)
    return(list())

  # Lower and upper bounds (dunif() or T()) and minimum segment length (see add_nodes()).
  # Constant change points and user-supplied initial values are bounded by themselves.
  env = stats::setNames(as.list(x_range), paste0("cp_", c(0, n_cp + 1)))
  evaluate = function(x, empty) if (is.null(x) || !nzchar(x)) empty else eval(str2lang(x), env, baseenv())
  is_free = !is.na(cp_rows$node)
  lower = upper = min_lengths = numeric(n_cp)
  for (i in seq_len(n_cp)) {
    fixed = c(user_inits[[cp_rows$parameter[i]]], if (cp_rows$kind[i] == "constant") as.numeric(cp_rows$value[i]), NA)[1]
    bounds = if (is_free[i]) get_prior_bounds(cp_rows$node_prior[[i]])
    min_lengths[i] = evaluate(cp_rows$min_length[i], 0)
    lower[i] = if (!is.na(fixed)) fixed else evaluate(bounds[1], -Inf)
    upper[i] = if (!is.na(fixed)) fixed else evaluate(bounds[2], Inf)
  }

  # Keep 1% of the range of x from the bounds and neighboring change points, but within narrow priors
  margin = pmin(0.01 * diff(x_range), (upper - lower) / 3)
  spacing = min_lengths + c(0, pmax(margin[-1], margin[-n_cp]))
  lower = lower + margin
  upper = upper - margin

  # Earliest and latest positions of ordered change points within these bounds
  for (i in seq_len(n_cp)[-1])
    lower[i] = max(lower[i], lower[i - 1] + spacing[i])
  for (i in rev(seq_len(n_cp)[-1]))
    upper[i - 1] = min(upper[i - 1], upper[i] - spacing[i])
  if (any(lower > upper))
    stop("The change point priors (and initial values) leave no room for ordered change points. Check the lower and upper bounds (in dunif() or T()), minimum segment lengths, and initial values of the change points.", call. = FALSE)

  target = x_range[1] + cumsum(min_lengths) + max((diff(x_range) - sum(min_lengths)) / (n_cp + 1), 0) * seq_len(n_cp)
  inits = stats::setNames(as.list(pmin(pmax(target, lower), upper)), cp_rows$parameter)
  inits[is_free & cp_rows$parameter %notin% names(user_inits)]
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
#' @param n_segments The number of segments.
#' @param data_columns Data-column roles, as returned by `mcp_columns()`.
#' @param coefs The `coefs` table from `get_coef_tables()`.
#' @param group_coefs The `group_coefs` table from `get_coef_tables()`.
#' @param parameters The parameters table from `get_pars_table()`.
#' @param design_specs Named fitted design specifications.
#' @noRd
get_jags_data = function(data, family, n_segments, data_columns, coefs, group_coefs, parameters, jags_code,
                          generated = TRUE, design_specs = list()) {
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

  # Predictor design matrix. Keep the JAGS data name for custom-code compatibility.
  jags_data$rhs_matrix_ = get_predictor_matrix(coefs, group_coefs)

  # Add any offset vectors to JAGS data
  for (design_spec in design_specs) {
    if (isTRUE(design_spec$has_offset)) {
      jags_data[[design_spec$offset_name]] = design_spec$offset_data
    }
  }

  # Compatibility for custom JAGS code written with the released data
  # constants. Generated code and newly resolved priors no longer need these.
  context = prior_context(data, data_columns, n_segments, design_specs)
  jags_data = legacy_prior_jags_data(jags_data, jags_code, context)

  if (generated && !is.null(data_columns$weights)) {  # Zeros trick (see jags_weighted_likelihood())
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
