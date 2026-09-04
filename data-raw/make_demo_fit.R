devtools::load_all()
empty = mcp_example("demo", sample = "none")
set.seed(78)
demo_fit = mcp(
    empty$model, data = empty$data, 
        warmup = 2000, iter = 6000, chains = 2, sample = "both", seed = 78
)

# Keep evenly spaced draws for a compact bundled example with reliable diagnostics.
thin_draws = function(draws, n = 500) {
    coda::mcmc.list(lapply(draws, function(chain) {
        indices = round((seq_len(n) - 0.5) * nrow(chain) / n)
        coda::mcmc(chain[indices, , drop = FALSE])
    }))
}
demo_fit[["mcmc_post"]] = thin_draws(.subset2(demo_fit, "mcmc_post"))
demo_fit[["mcmc_prior"]] = thin_draws(.subset2(demo_fit, "mcmc_prior"))

# Save to mcp
usethis::use_data(demo_fit, overwrite = TRUE, compress = "xz")
