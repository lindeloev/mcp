ex = mcp_example("demo", sample = FALSE)
demo_fit = mcp(ex$model, data = ex$data, adapt = 3000, iter = 1000, sample = "both")

# Save to mcp
usethis::use_data(demo_fit, overwrite = TRUE)
