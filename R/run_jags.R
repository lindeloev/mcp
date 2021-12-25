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
                    cores,
                    sample,
                    n.chains,
                    n.iter,
                    n.adapt,
                    inits
) {

  # Prevent failure of all mcp methods when length(pars) <= 2.
  # This always happens when there is only one parameter, so we just
  # save samples from the dummy change points.
  if (length(pars) <= 2)
    pars = c(pars, "cp_0", "cp_1")

  # Set number of cores from "all" or mc.cores if `cores` is not specified.
  # Max at 2 for CRAN etc.
  opts_cores = options()$mc.cores
  if (is.numeric(opts_cores) && cores == 1)
    cores = opts_cores
  if (cores == "all") {
    cores = future::availableCores() - 1
    n.chains = cores
  }
  if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") == "TRUE") {
    if (cores > 1)
      cores = 2
  }

  # Define the sampling function in this environment.
  # Can be used sequentially or in parallel.
  do_sampling = function(seed, n.chains, quiet) {
    # Optionally seed JAGS. Typically for parallel processing to avoid risk of identical seeds.
    if (!is.null(seed))
      inits = c(inits, list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed))

    # Compile model
    jm = rjags::jags.model(
      file = textConnection(jags_code),
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

  # Time for sampling!
  if (cores == 1) {
    # SERIAL
    timer = proc.time()
    samples = try(do_sampling(
      seed = NULL,
      n.chains = n.chains,
      quiet = FALSE
    ))

  } else if (cores == "all" || cores > 1) {
    # PARALLEL using the future package and one chain per worker
    if (future::nbrOfWorkers() == 1)
      message("Setting up parallel workers...")
    future::plan(future::multisession, .skip = TRUE)

    message("Parallel sampling in progress...")
    timer = proc.time()
    samples = future.apply::future_lapply(
      sample(1:1000, n.chains),  # Random seed to JAGS
      n.chains = 1,
      quiet = TRUE,
      FUN = do_sampling,
      future.seed = TRUE
    )

    # Get result as mcmc.list
    samples = unlist(samples, recursive = FALSE)
    class(samples) = "mcmc.list"
  }

  # Sampling finished
  passed = proc.time() - timer
  message("Finished sampling in ", round(passed["elapsed"], 1), " seconds\n")

  # Recover the levels of varying effects if it succeeded
  if (coda::is.mcmc.list(samples)) {
    return(samples)
  } else {
    # If it didn't succeed, quit gracefully.
    warning("--------------\nJAGS failed with the above error. Returning an `mcpfit` without samples. Inspect fit$prior and fit$jags_code to identify the problem. Read about typical problems and fixes here: https://lindeloev.github.io/mcp/articles/tips.html.")
    return(NULL)
  }
}


#' Adds helper variables for use in `run_jags`
#'
#' Returns the relevant data columns as a list and add elements with unique
#' varying group levels.
#'
#' @aliases get_jags_data
#' @keywords internal
#' @noRd
#' @inheritParams run_jags
#' @param data A tibble
#' @param ST A segment table (tibble), returned by `get_segment_table`.
#' @param rhs_table Returned by `get_rhs()`
get_jags_data = function(data, family, ST, rhs_table, jags_code, sample) {
  cols_varying = unique(stats::na.omit(ST$cp_group_col))

  # Start with "raw" data
  cols_data = unique(stats::na.omit(c(ST$y, ST$x, ST$trials, ST$weights)))
  jags_data = as.list(data[, c(cols_varying, cols_data)])

  for (col in cols_varying) {
    # Add meta-data (now many varying group levels)
    tmp = paste0("n_unique_", col)
    jags_data[[tmp]] = length(unique(data[, col]))

    # Make varying columns numeric in order of appearance
    # They will be recovered using the recover_levels()
    jags_data[[col]] = as.numeric(factor(jags_data[[col]], levels = unique(jags_data[[col]])))
  }

  # RHS data matrix
  jags_data$rhs_matrix_ = get_rhs_matrix(rhs_table)

  # Add e.g. MINX = min(data$x) for all variables where they are used.
  # Searches whether it is in jags_code. If yes, add to jags_data
  # TO DO: there must be a more concise way than this...
  funcs = c("min", "max", "sd", "mean", "median", "mad")

  # For x
  for (func in funcs) {
    constant_name = toupper(paste0(func, "x"))  # No link function for x
    if (stringr::str_detect(jags_code, constant_name)) {
      column = ST$x[1]  # from x/y to data column name
      linkdata = data[, column]  # No link function on x stuff
      jags_data[[constant_name]] = get(func)(linkdata, na.rm = TRUE)  # Call func
    }
  }

  # For y
  for (func in funcs) {
    constant_name =  toupper(paste0(func, "LINK", "y"))  # Link function applies to y
    if (stringr::str_detect(jags_code, constant_name)) {
      column = ST$y[1]  # from x/y to data column name
      linkdata = family$linkfun(data[, column])  # No link function on x stuff
      jags_data[[constant_name]] = get(func)(linkdata, na.rm = TRUE)  # Call func
    }
  }

  # For default prior
  if (stringr::str_detect(jags_code, "N_CP"))
    jags_data$N_CP = nrow(ST) - 1

  # Set response = NA if we only sample prior
  if (sample == "prior")
    jags_data[[ST$y[1]]] = rep(NA, nrow(data))

  # Return
  jags_data
}


#' Recover the levels of varying effects in mcmc.list
#'
#' Jags uses 1, 2, 3, ..., etc. for indexing of varying effects.
#' This function adds back the original levels, whether numeric or string
#'
#' @aliases recover_levels
#' @keywords internal
#' @noRd
#' @param samples An mcmc.list with varying columns starting in `mcmc_col`.
#' @param data A tibble or data.frame
recover_levels = function(samples, data, ST) {
  for (i in seq_len(nrow(ST))) {
    S = ST[i, ]
    if (!is.na(S$cp_group_col)) {
      # Get vectors of old ("from") and replacement column names in samples
      from = colnames(samples[[1]])[stringr::str_starts(colnames(samples[[1]]), paste0(S$cp_group, '\\['))]  # Current column names
      to = sprintf(paste0(S$cp_group, '[%s]'), unique(data[, S$cp_group_col]))  # Desired column names

      # Recode column names on each list (chain) using lapply
      names(to) = from
      samples = lapply(samples, function(x) {
        colnames(x) = dplyr::recode(colnames(x), !!!to)
        x
      })
    }
  }

  samples
}
