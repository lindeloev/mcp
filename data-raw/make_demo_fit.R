devtools::load_all()
empty = mcp_example("demo", sample = FALSE)
set.seed(40)
demo_fit = mcp(
    empty$model, data = empty$data, 
    warmup = 3000, iter = 1000, chains = 2, sample = "both",
    diagnostics = FALSE, seed = 40
)

# Save to mcp
usethis::use_data(demo_fit, overwrite = TRUE)
