ex = mcp_example("demo")
demo_fit = mcp(ex$model, data = ex$data, adapt = 3000, iter = 1000, sample = "both")
demo_fit$mcmc_loglik = NULL  # Make the object small

# Save to mcp
usethis::use_data(demo_fit, overwrite = TRUE)
