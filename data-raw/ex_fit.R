# Define the model
segments = list(
  response ~ 1,
  ~ 0 + time,
  ~ 1 + time
)

# Fit it. The `ex_demo` dataset is included in mcp
options(mc.cores = 3)
ex_fit = mcp(segments, data = ex_demo, adapt = 3000, iter = 1000, sample = "both")
ex_fit$mcmc_loglik = NULL  # Make the object small

# Save to mcp
usethis::use_data(ex_fit, overwrite = TRUE)
