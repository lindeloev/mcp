# Define model
segments = list(
  y ~ 1,
  ~ 1
)
empty = mcp::mcp(segments, sample = FALSE, par_x = "x")

# Simulate data
set.seed(40)
ex_plateaus = tibble::tibble(
  x = runif(100, 0, 100),
  y = empty$func_y(x, cp_1 = 50, int_1 = 10, int_2 = 20, sigma = 8)
)

# Save to mcp
usethis::use_data(ex_plateaus, overwrite = TRUE)
