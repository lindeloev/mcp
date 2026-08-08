devtools::load_all()
ex = mcp_example("demo", sample = FALSE)
set.seed(40)
demo_fit = mcp(ex$model, data = ex$data, adapt = 3000, iter = 1000, chains = 2, sample = "both", seed = 40)

# Save to mcp
usethis::use_data(demo_fit, overwrite = TRUE)
