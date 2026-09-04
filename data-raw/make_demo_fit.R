make_demo_fit = function() {
    # Set reproducible state
    previous_plan = future::plan(future::sequential)
    on.exit(future::plan(previous_plan), add = TRUE)

    # Load latest version of mcp
    devtools::load_all()

    # Demo fit!
    empty = mcp_example("demo", sample = "none")
    for (i in seq_along(empty$model)) environment(empty$model[[i]]) = globalenv()  # Do not include environment
    demo_fit = mcp(
        empty$model, data = empty$data,
        warmup = 2000, iter = 6000, chains = 2,
        sample = "both", seed = 78
    )

    # Thinning for max ESS/draws to keep mcp package compact
    thin_draws = function(draws, n = 500) {
        coda::mcmc.list(lapply(draws, function(chain) {
            indices = as.integer(floor((seq_len(n) - 1) * nrow(chain) / n) + 1)
            coda::mcmc(chain[indices, , drop = FALSE])
        }))
    }
    demo_fit[["mcmc_post"]] = thin_draws(.subset2(demo_fit, "mcmc_post"))
    demo_fit[["mcmc_prior"]] = thin_draws(.subset2(demo_fit, "mcmc_prior"))

    # Normalize closure environments and strip srcrefs for compact, reproducible packaging
    demo_fit = rapply(demo_fit, function(f) {
        f = utils::removeSource(f)
        if (is.environment(environment(f))) {
            environment(f) = if (identical(formals(f), formals(stats::gaussian()$variance))) asNamespace("stats") else asNamespace("mcp")
        }
        f
    }, classes = "function", how = "replace")

    # Save to mcp
    usethis::use_data(demo_fit, overwrite = TRUE, compress = "xz")
    invisible(demo_fit)
}
