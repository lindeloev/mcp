# ABOUT: These functions call the JAGS sampler and set up inputs in a way
# that is appropriate for that, including data (get_jags_data()).
# ------------------------------------------------------------------------


#' Run parallel MCMC sampling using JAGS.
#'
#'
#' @aliases run_jags
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @inheritParams rjags::jags.model
#' @inheritParams rjags::coda.samples
#' @param jags_code A string. JAGS model, usually returned by `get_jags_code()`.
#' @param pars Character vector of parameters to save/monitor.
#' @return `mcmc.list``
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
run_jags = function(data,
                    jags_code,
                    jags_data,
                    pars,
                    sample,
                    n.chains,
                    n.iter,
                    n.adapt,
                    inits,
                    seed,
                    quiet
) {

  # Prevent failure of all mcp methods when length(pars) <= 2.
  # This always happens when there is only one parameter, so we just
  # save samples from the dummy change points.
  if (length(pars) <= 2)
    pars = c(pars, "cp_0", "cp_1")

  # Define the sampling function in this environment.
  # Can be used sequentially or in parallel.
  do_sampling = function(inits, n.chains) {
    sample_jags = function() {
      # Compile model
      jags_connection = textConnection(jags_code)
      on.exit(close(jags_connection), add = TRUE)
      jm = rjags::jags.model(
        file = jags_connection,
        data = jags_data,
        inits = inits,
        n.chains = n.chains,
        n.adapt = n.adapt,
        quiet = quiet
      )

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
      samples = NULL
      save_samples = function() samples <<- sample_jags()
      capture.output(save_samples())
    } else {
      samples = sample_jags()
    }
    samples
  }

  # Use JAGS directly under a sequential future plan. This compiles the model
  # only once for all chains and preserves JAGS's progress output.
  n_workers = future::nbrOfWorkers()
  timer = proc.time()
  if (n_workers == 1) {
    inits = get_jags_inits(inits, seed, n.chains, sample)
    samples = try(do_sampling(
      inits = inits,
      n.chains = n.chains
    ))
  } else {
    # Submit one chain per future. The user's future plan controls the backend.
    if (!quiet)
      message("Parallel sampling in progress...")
    if (is.null(seed))
      seed = sample.int(.Machine$integer.max - 2 * n.chains, 1)
    inits = get_jags_inits(inits, seed, n.chains, sample)
    samples = future.apply::future_lapply(
      inits,
      n.chains = 1,
      FUN = do_sampling,
      future.seed = TRUE
    )

    # Get result as mcmc.list
    samples = unlist(samples, recursive = FALSE)
    class(samples) = "mcmc.list"
  }

  # Sampling finished
  passed = proc.time() - timer
  if (!quiet)
    message("Finished sampling in ", round(passed["elapsed"], 1), " seconds\n")

  # Recover the levels of group-level effects if it succeeded
  if (coda::is.mcmc.list(samples)) {
    return(samples)
  } else {
    # If it didn't succeed, quit gracefully.
    warning("--------------\nJAGS failed with the above error. Returning an `mcpfit` without samples. Inspect fit$prior and fit$jags_code to identify the problem. Read about typical problems and fixes here: https://lindeloev.github.io/mcp/articles/tips.html.")
    return(NULL)
  }
}


#' Add deterministic JAGS RNG initial values
#'
#' @keywords internal
#' @noRd
#' @param inits Initial values passed to `mcp()`.
#' @param seed A positive integer or `NULL`.
#' @param n.chains Number of chains.
#' @param sample One of `"post"` or `"prior"`.
#' @return Initial values in list-of-lists form when `seed` is supplied.
get_jags_inits = function(inits, seed, n.chains, sample) {
  if (is.null(seed))
    return(inits)

  if (length(inits) > 0 && all(vapply(inits, is.list, logical(1))))
    stop(
      "When `seed` is supplied, `inits` must be a single named list ",
      "shared by all chains, not a list of chain-specific lists."
    )

  if (is.null(inits))
    inits = list()
  inits[c(".RNG.name", ".RNG.seed", ".RNG.state")] = NULL

  rng_seed = seq_len(n.chains)
  if (sample == "prior")
    rng_seed = rng_seed + n.chains
  rng_seed = as.integer(
    ((as.double(seed) - 1 + rng_seed - 1) %% .Machine$integer.max) + 1
  )

  lapply(
    rng_seed,
    function(seed) {
      c(inits, list(
        .RNG.name = "base::Wichmann-Hill",
        .RNG.seed = seed
      ))
    }
  )
}


#' Adds helper variables for use in `run_jags`
#'
#' Returns the relevant data columns as a list and add elements with unique
#' grouping-factor levels.
#'
#' @aliases get_jags_data
#' @keywords internal
#' @noRd
#' @inheritParams run_jags
#' @param data A tibble
#' @param segments A segment table returned by `get_segment_tables()`.
#' @param predictors Returned by `get_predictors()`.
#' @param group_effects Returned by `get_group_effects()`.
get_jags_data = function(data, family, segments, predictors, group_effects, jags_code) {
  group_cols = unique(stats::na.omit(group_effects$group_col))

  # Start with "raw" data
  aux_columns = get_family_aux_columns(family, segments)
  cols_data = unique(stats::na.omit(c(segments$y, segments$x, unname(aux_columns))))
  jags_data = as.list(data[, c(group_cols, cols_data)])

  for (col in group_cols) {
    # Add metadata for the grouping-factor levels.
    tmp = paste0("n_unique_", col)
    jags_data[[tmp]] = length(unique(data[, col]))

    # Make grouping columns numeric in order of appearance.
    # They will be recovered using the recover_levels()
    jags_data[[col]] = as.numeric(factor(jags_data[[col]], levels = unique(jags_data[[col]])))
  }

  # Predictor design matrix. Keep the JAGS data name for custom-code compatibility.
  jags_data$rhs_matrix_ = get_predictor_matrix(predictors, group_effects)

  # Add named data constants for change point prior bounds (see jagsify_constants())
  jags_constants = attr(jags_code, "jags_constants")
  if (!is.null(jags_constants))
    jags_data[names(jags_constants)] = jags_constants

  # Compatibility for custom JAGS code written with the released data
  # constants. Generated code and newly resolved priors no longer need these.
  context = prior_context(data, segments)
  jags_data = add_legacy_prior_jags_data(jags_data, jags_code, context)

  # Return
  jags_data
}


#' Recover the levels of group-level effects in mcmc.list
#'
#' JAGS uses 1, 2, 3, ... for indexing group-level effects.
#' This function adds back the original levels, whether numeric or string
#'
#' @aliases recover_levels
#' @keywords internal
#' @noRd
#' @param samples An mcmc.list with group-level columns starting in `mcmc_col`.
#' @param data A tibble or data.frame
#' @param group_effects Returned by `get_group_effects()`.
recover_levels = function(samples, data, group_effects) {
  for (i in seq_len(nrow(group_effects))) {
    effect = group_effects[i, ]
    # Get vectors of old ("from") and replacement column names in samples
    from = colnames(samples[[1]])[stringr::str_starts(colnames(samples[[1]]), paste0(effect$name, '\\['))]
    to = sprintf(paste0(effect$name, '[%s]'), unique(data[, effect$group_col]))

    # Recode column names on each list (chain) using lapply
    names(to) = from
    samples = lapply(samples, function(x) {
      colnames(x) = dplyr::recode(colnames(x), !!!to)
      x
    })
  }

  samples
}
