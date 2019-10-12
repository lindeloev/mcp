#' Run parallel MCMC sampling using JAGS.
#'
#'
#' @aliases run_jags
#' @inheritParams dclone::jags.fit
#' @param jags_code A string. JAGS model, usually returned by \code{make_jagscode()}.
#' @param model_file A temporary file. Makes parallel sampling possible
#' @param n.chains Number of chains to run. Defaults to all-but-one cores.
#' @param n.iter Number of post-warmup samples to draw.
#' @param n.adapt Number of iterations to adapt sampler values.
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
                    model_file = "tmp_jags_code.txt",

                    # JAGS arguments
                    n.chains = parallel::detectCores() - 1,  # Use all cores but one
                    n.iter = 3000,  # Number of iterations post-warmup.
                    n.adapt = 1500,  # Takes some time to adapt
                    n.update = 1500,  # Same as n.adapt
                    ...  # Otherwise run with default JAGS settings
) {

  # Write model to disk
  sink(model_file)
  cat(jags_code)
  sink()  # stops sinking :-)

  # Start parallel cluster and timer
  cl = parallel::makePSOCKcluster(n.chains)
  timer = proc.time()

  # Do the sampling. Yield mcmc.list
  samples = dclone::jags.parfit(
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

  # Stop the cluster, delete the file, and show time
  parallel::stopCluster(cl)
  file.remove(model_file)
  print(proc.time() - timer)  # Print time

  # Return
  samples
}
