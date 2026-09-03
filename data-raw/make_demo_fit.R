devtools::load_all()
empty = mcp_example("demo", sample = "none")
set.seed(78)
demo_fit = mcp(
    empty$model, data = empty$data, 
    warmup = 3000, iter = 1000, chains = 2, sample = "both",
    diagnostics = FALSE, seed = 78
)

# Save to mcp
usethis::use_data(demo_fit, overwrite = TRUE)
