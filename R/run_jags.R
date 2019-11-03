#' Run parallel MCMC sampling using JAGS.
#'
#'
#' @aliases run_jags
#' @inheritParams mcp
#' @inheritParams dclone::jags.fit
#' @param jags_code A string. JAGS model, usually returned by \code{make_jagscode()}.
#' @param ST A segment table (tibble), returned by \code{get_segment_table}.
#'   Only really used when the model contains varying effects.
#' @param model_file A temporary file. Makes parallel sampling possible
#'   samples that are discarded between n.adapt and n.iter (to improve convergence)..
#' @param ... Parameters for \code{jags.parfit} which channels them to \code{jags.fit}.
#' @return \code{mcmc.list}
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' \dontrun{
#' run_jags(data, model, params)
#'}

run_jags = function(data,
                    jags_code,
                    params,
                    ST,
                    cores,
                    sample = sample,
                    model_file = "tmp_jags_code.txt",
                    ...  # Otherwise run with default JAGS settings
) {

  # Prevent failure of all mcp methods when length(params) <= 2 (one parameter +
  # loglik_).This always happens when there is only one segment, so we just save samples
  # from the dummy change points.
  if (length(params) <= 2)
    params = c(params, "cp_0", "cp_1")

  # Start timer
  timer = proc.time()

  if (cores == 1) {
    # SERIAL
    samples = try(dclone::jags.fit(
      data = get_jags_data(data, ST, jags_code, sample),
      params = params,
      model = textConnection(jags_code),
      ...
    ))
  } else if (cores == "all" | cores > 1) {
    # PARALLEL
    # Write model to disk
    sink(model_file)
    cat(jags_code)
    sink()  # stops sinking :-)

    # Start parallel cluster
    if (cores == "all") {
      cores = parallel::detectCores() - 1
    }
    cl = parallel::makePSOCKcluster(cores)

    # Do the sampling. Yield mcmc.list
    samples = try(dclone::jags.parfit(
      cl = cl,
      data = get_jags_data(data, ST, jags_code, sample),
      params = params,
      model = model_file,
      ...
    ))

    # Stop the cluster, delete the file
    parallel::stopCluster(cl)
    file.remove(model_file)
  }

  if(coda::is.mcmc.list(samples)) {  # If sampling succeeded
    # Recover the levels of varying effects
    for (i in seq_len(nrow(ST))) {
      S = ST[i, ]
      if (!is.na(S$cp_group_col)) {
        samples = recover_levels(samples, data, S$cp_group, S$cp_group_col)
      }
    }
    # Return
    print(proc.time() - timer)  # Return time
    return(samples)
  } else {
    warning("--------------\nJAGS failed with the above error. This is often caused by priors which allow for impossible values, such as negative standard deviations, probabilities outside [0, 1], etc.\n\nAnother class of error is directed cycles: when using a parameter to set a prior for another parameter, it may end up ultimately predicting itself.\n\nReturning an `mcpfit` without samples. Inspect fit$prior and fit$jags_code to identify the problem.")
    return(NULL)
  }
}


#' Adds helper variables for use in \code{run_jags}
#'
#' Returns the relevant data columns as a list and add elements with unique
#' varying group levels.
#'
#' @param data A tibble
#' @param ST A segment table (tibble), returned by \code{get_segment_table}.

get_jags_data = function(data, ST, jags_code, sample) {
  cols_varying = unique(stats::na.omit(ST$cp_group_col))

  # Start with "raw" data
  cols_data = unique(stats::na.omit(c(ST$y, ST$x, ST$trials)))
  jags_data = as.list(data[, c(cols_varying, cols_data)])

  for (col in cols_varying) {
    # Add meta-data (now many varying group levels)
    tmp = paste0("n_unique_", col)
    jags_data[[tmp]] = length(unique(dplyr::pull(data, col)))

    # Make varying columns numeic in order of appearance
    # They will be recovered using the recover_levels()
    jags_data[[col]] = as.numeric(factor(jags_data[[col]], levels = unique(jags_data[[col]])))
  }


  # Add e.g. MINX = min(data$x) for all variables where they are used.
  # Searches whether it is in jags_code. If yes, add to jags_data
  # TO DO: there must be a prettier way to do this.
  funcs = c("min", "max", "sd", "mean")
  xy_vars = c("x", "y")
  for (func in funcs) {
    for(xy_var in xy_vars) {
      constant_name = toupper(paste0(func, xy_var))
      if (stringr::str_detect(jags_code, constant_name)) {
        func_eval = eval(parse(text = func))  # as real function
        jags_data[[constant_name]] = func_eval(dplyr::pull(data, ST[, xy_var][[1]][1]))
      }
    }
  }

  # Set response = NA if we only sample prior
  if(sample == "prior")
    jags_data[[ST$y[1]]] = rep(NA, nrow(data))

  # Return
  jags_data
}



#' Recover the levels of varying effects in mcmc.list
#'
#' Jags uses 1, 2, 3, ..., etc. for indexing of varying effects.
#' This function adds back the original levels, whether numeric or string
#'
#' @param samples An mcmc.list with varying columns starting in \code{mcmc_col}.
#' @param data A tibble or data.frame with the cols in \code{data_col}.
#' @param mcmc_col A vector of strings.
#' @param data_col A vector of strings. Has to be same length as \code{mcmc_col}
#'
recover_levels = function(samples, data, mcmc_col, data_col) {
  # Get vectors of old ("from") and replacement column names in samples
  from = colnames(samples[[1]])[stringr::str_starts(colnames(samples[[1]]), paste0(mcmc_col, '\\['))]  # Current column names
  to = sprintf(paste0(mcmc_col, '[%s]'), unique(dplyr::pull(data, data_col)))  # Desired column names

  # Recode column names on each list (chain) using lapply
  names(to) = from
  lapply(samples, function(x) {
    colnames(x) = dplyr::recode(colnames(x), !!!to)
    x
  })
}
