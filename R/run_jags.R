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
#' @param cp_inits Named list. Default initial values for change points, from `get_cp_inits()`.
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
                    cp_inits = list(),
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
      jm = rjags::jags.model(
        file = jags_connection,
        data = jags_data,
        inits = inits,
        n.chains = n.chains,
        n.adapt = 0,
        quiet = TRUE
      )
      progress_bar = if (quiet) "none" else getOption("jags.pb")
      if (!rjags::adapt(jm, n.iter = n.adapt, end.adaptation = TRUE, progress.bar = progress_bar))
        warning("Adaptation incomplete")

      # Ensure requested warmup iterations are discarded even if JAGS skips adaptation
      burnin = n.adapt - jm$iter()
      if (burnin > 0) {
        stats::update(
          object = jm,
          n.iter = burnin,
          progress.bar = progress_bar
        )
      }

      # Sample and return
      rjags::coda.samples(
        model = jm,
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

  # Default inits where users supplied none: ordered change points (see get_cp_inits()) and
  # prior medians for other nodes, e.g., SDs away from zero.
  # Custom JAGS code may define change points differently, so only initialize sampled nodes.
  is_sampled = vapply(names(cp_inits), function(name) grepl(paste0("\\n\\s*", name, " ~"), jags_code), logical(1))
  defaults = c(cp_inits[is_sampled], get_jags_median_inits(jags_code))
  defaults = defaults[!duplicated(names(defaults))]
  inits = lapply(inits, function(ch) c(ch, defaults[setdiff(names(defaults), names(ch))]))

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


# The JAGS glm module helps some families always and others only with group-level predictor effects.
use_jags_glm = function(family, group_effects) {
  glm_module = family$backends$jags$glm_module
  identical(glm_module, "always") || (identical(glm_module, "group_effects") && any(group_effects$part != "cp"))
}


# JAGS data for sampling the prior. Stochastic response nodes are set to NA, except for weighted
# discrete likelihoods, which use the response in a zeros trick that `response_observed_ = 0` turns off.
get_jags_prior_data = function(jags_data, jags_code, response, custom_jags_code) {
  response_is_stochastic = if (custom_jags_code)
    !grepl("likelihood_phi_", jags_code, fixed = TRUE) else
    grepl(paste0(response, "[i_] ~"), jags_code, fixed = TRUE)
  if (response_is_stochastic)
    jags_data[[response]][] = NA
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
    ch = inits[[i]]
    ch[c(".RNG.name", ".RNG.seed", ".RNG.state")] = NULL
    c(ch, list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = rng_seed[i]))
  })
}


# Check that `inits` only set parameters that are sampled directly
assert_inits = function(inits, pars, jags_code) {
  is_sampled = vapply(pars, function(name) grepl(paste0("\\n\\s*", name, "(\\[[^]]*\\])? ~"), jags_code), logical(1))
  invalid = setdiff(unlist(lapply(if (is.list(inits[[1]])) inits else list(inits), names)), pars[is_sampled])
  if (length(invalid) > 0)
    stop("`inits` can only set these parameters of this model: `", paste0(pars[is_sampled], collapse = "`, `"), "`. Got `", paste0(unique(invalid), collapse = "`, `"),
         "`. Fixed parameters and parameters that are reparameterized for faster sampling (see `mcp_reparam` in ?mcp-package) cannot be initialized.", call. = FALSE)
}


# Initial values at the prior medians. Covers scalar stochastic nodes in JAGS code whose prior
# arguments and truncation bounds are numbers. JAGS initializes the remaining nodes.
get_jags_median_inits = function(jags_code) {
  nodes = stringr::str_match_all(jags_code, "\\n\\s*([A-Za-z0-9_.]+) ~ ([^\\n#]+)")[[1]]
  as_number = function(x) tryCatch(as.numeric(eval(str2lang(x), baseenv())), error = function(e) NA_real_)

  inits = list()
  for (i in seq_len(nrow(nodes))) {
    parts = split_prior_truncation(nodes[i, 3])
    call = parse_prior_call(parts$distribution)
    if (is.null(call))
      next

    # Numeric arguments with JAGS precision converted back to scale
    a = vapply(call$args, as_number, numeric(1), USE.NAMES = FALSE)
    exponent = jags_scale_exponent(call$name)
    if (!is.null(exponent))
      a[2] = a[2]^(-1 / exponent)

    bounds = c(-Inf, Inf)
    if (!is.null(parts$truncation)) {
      trunc_args = parse_prior_call(parts$truncation)$args
      bounds[nzchar(trunc_args)] = vapply(trunc_args[nzchar(trunc_args)], as_number, numeric(1), USE.NAMES = FALSE)
    }
    if (!anyNA(a) && !anyNA(bounds))
      inits[[nodes[i, 2]]] = prior_quantile(call$name, a, bounds, 0.5)
  }
  inits
}


# Initial values for the change points with a prior distribution, except those in `user_inits`.
# Equispaced in the observed range of x (`x_range`) with room for minimum segment lengths, then moved the least
# possible to be within the lower and upper bounds of their priors and ordered. Fixed change points and the
# user's initial values stay where they are. Sampler fails at init if impossible values are chosen because
# it is incompatible with the 1 ~ dbern(...) ordering node.
# - prior: Named list of resolved priors, as in `fit$prior`.
get_cp_inits = function(prior, x_range, user_inits = list()) {
  cps = prior[grepl("^cp_[0-9]+$", names(prior))]
  n = length(cps)
  if (n == 0)
    return(list())

  # Lower and upper bounds (dunif() or T()) and minimum segment length.
  # Fixed change points and user-supplied initial values are bounded by themselves.
  env = stats::setNames(as.list(x_range), paste0("cp_", c(0, n + 1)))
  evaluate = function(x, empty) if (is.null(x) || !nzchar(x)) empty else eval(str2lang(x), env, baseenv())
  is_free = vapply(cps, is_distribution, logical(1))
  lower = upper = min_lengths = numeric(n)
  for (i in seq_len(n)) {
    fixed = c(user_inits[[names(cps)[i]]], suppressWarnings(as.numeric(cps[[i]])))[1]
    bounds = if (is_free[i] && !startsWith(cps[[i]], "dirichlet")) split_cp_bounds(cps[[i]])$bounds
    min_length = if (!is.null(bounds)) get_min_segment_length(bounds, i, names(cps))
    if (!is.null(min_length))
      bounds[1] = "cp_0"
    min_lengths[i] = evaluate(min_length, 0)
    lower[i] = if (!is.na(fixed)) fixed else evaluate(bounds[1], -Inf)
    upper[i] = if (!is.na(fixed)) fixed else evaluate(bounds[2], Inf)
  }

  # Keep 1% of the range of x from the bounds and neighboring change points, but within narrow priors
  margin = pmin(0.01 * diff(x_range), (upper - lower) / 3)
  spacing = min_lengths + c(0, pmax(margin[-1], margin[-n]))
  lower = lower + margin
  upper = upper - margin

  # Earliest and latest positions of ordered change points within these bounds
  for (i in seq_len(n)[-1])
    lower[i] = max(lower[i], lower[i - 1] + spacing[i])
  for (i in rev(seq_len(n)[-1]))
    upper[i - 1] = min(upper[i - 1], upper[i] - spacing[i])
  if (any(lower > upper))
    stop("The change point priors (and initial values) leave no room for ordered change points. Check the lower and upper bounds (in dunif() or T()), minimum segment lengths, and initial values of the change points.", call. = FALSE)

  target = x_range[1] + cumsum(min_lengths) + max((diff(x_range) - sum(min_lengths)) / (n + 1), 0) * seq_len(n)
  inits = stats::setNames(as.list(pmin(pmax(target, lower), upper)), names(cps))
  inits[is_free & names(cps) %notin% names(user_inits)]
}


# Check names of user-data nodes passed to JAGS against generated JAGS nodes.
# - generated: Whether `jags_code` was generated by mcp (`FALSE` for custom `jags_code`).
assert_jags_namespace = function(data_names, family, segments, predictors,
                                  group_effects, jags_code, generated = TRUE) {
  group_cols = unique(stats::na.omit(group_effects$group_col))
  dpars = family$dpar_specs$dpar

  # Public parameters and fixed data range endpoints share one namespace.
  parameter_names = stats::na.omit(c(
    paste0("cp_", 0:nrow(segments)),
    predictors$name,
    group_effects$name, group_effects$sd_name
  ))

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
      paste0("x_local_", seq_len(nrow(segments)), "_"),
      paste0("link_", dpars, "_"), paste0(dpars, "_"),
      paste0(c(predictors$name, group_effects$name), "_inverse")
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
#' @param segments A segment table returned by `get_segment_tables()`.
#' @param predictors Returned by `get_predictors()`.
#' @param group_effects Returned by `get_group_effects()`.
#' @param design_specs Named fitted design specifications.
#' @noRd
get_jags_data = function(data, family, segments, predictors, group_effects, jags_code,
                          series = NULL, generated = TRUE, design_specs = list()) {
  group_cols = unique(stats::na.omit(group_effects$group_col))

  # Start with "raw" data
  aux_columns = get_family_aux_columns(family, segments)
  cols_data = unique(stats::na.omit(c(segments$y, segments$x, unname(aux_columns))))
  assert_jags_namespace(
    unique(c(names(data), group_cols, cols_data, series)), family, segments,
    definition_rows(predictors), definition_rows(group_effects), jags_code, generated
  )
  jags_data = as.list(data[, cols_data])

  if (!is.null(series))
    jags_data$series_id_ = as.integer(factor(data[[series]], levels = unique(data[[series]])))

  for (col in group_cols) {
    # Add metadata for the grouping-factor levels.
    tmp = paste0("n_unique_", col)
    group_levels = levels(factor(data[[col]]))
    jags_data[[tmp]] = length(group_levels)

    # Index grouping levels in a separate node so the grouping column can also
    # be par_x, the response, etc. Levels are recovered using recover_levels()
    jags_data[[paste0(col, "_index_")]] = as.numeric(factor(data[[col]], levels = group_levels))
  }

  # Predictor design matrix. Keep the JAGS data name for custom-code compatibility.
  jags_data$rhs_matrix_ = get_predictor_matrix(predictors, group_effects)

  # Add any offset vectors to JAGS data
  for (spec in design_specs) {
    if (isTRUE(spec$has_offset)) {
      jags_data[[spec$offset_name]] = spec$offset_data
    }
  }

  # Compatibility for custom JAGS code written with the released data
  # constants. Generated code and newly resolved priors no longer need these.
  context = prior_context(data, segments, design_specs)
  jags_data = legacy_prior_jags_data(jags_data, jags_code, context)

  if (generated && grepl("response_observed_", jags_code, fixed = TRUE))
    jags_data$response_observed_ = as.integer(!is.na(data[[segments$y[1]]]))
  if (generated && grepl("likelihood_zero_", jags_code, fixed = TRUE))
    jags_data$likelihood_zero_ = numeric(nrow(data))

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
recover_levels = function(draws, data, group_effects) {
  # JAGS names levels by index, e.g., "x_1_id[2]"; use the original levels, e.g., "x_1_id[Omar]"
  index_names = character()
  level_names = character()
  for (i in seq_len(nrow(group_effects))) {
    group_levels = levels(factor(data[[group_effects$group_col[i]]]))
    index_names = c(index_names, paste0(group_effects$name[i], "[", seq_along(group_levels), "]"))
    level_names = c(level_names, paste0(group_effects$name[i], "[", group_levels, "]"))
  }

  lapply(draws, function(x) {
    matched = match(colnames(x), index_names)
    colnames(x)[!is.na(matched)] = level_names[matched[!is.na(matched)]]
    x
  })
}
