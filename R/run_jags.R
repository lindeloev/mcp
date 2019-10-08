library(dclone)

#' Run parallel MCMC sampling using JAGS.
#'
#' @param data A data.frame or tibble containing the variables expressed in `model` in long format.
#' @param model A JAGS model, usually returned by `make_jagscode()`.
#' @param model_file A temporary file. Makes parallel sampling possible
#' @param ... Parameters for `jags.parfit` which channels them to `jags.fit`.
#' @keywords mcmc, jags, mct
#' @export
#' @examples
#' run_jags(data, model, params)
run_jags = function(data,
                    model,
                    params,
                    model_file = "tmp_jags_model.txt",
                    
                    # JAGS arguments
                    n.chains = parallel::detectCores() - 1,  # Use all cores but one
                    n.iter = 2000,  # This model has a high number of effective samples
                    n.adapt = 2000,  # This model is quick to adapt
                    n.update = 2000,  # This model is quick to adapt
                    ...  # Otherwise run with default JAGS settings
) {
  
  # Write model to disk
  sink(model_file)
  cat(model)
  sink()  # stops sinking :-)
  
  # Start parallel cluster and timer
  cl = makePSOCKcluster(n.chains)
  timer = proc.time()
  
  # Do the sampling. Yield mcmc.list
  samples = jags.parfit(
    cl = cl,
    data = data,
    params = params,
    model = model_file,
    n.chains = n.chains,
    n.iter = n.iter,
    n.adapt = n.adapt,
    n.update = n.update,
    ...
  )
  mcmc = as.mcmc.list(samples)
  
  # Stop the cluster, delete the file, and show time
  stopCluster(cl)
  file.remove(model_file)
  print(proc.time() - timer)  # Print time
  
  # Return
  mcmc
}