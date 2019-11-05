# Define the model
segments = list(
  y ~ 1 + x,
  1 ~ rel(1) + rel(x),
  rel(1) ~ 0 + x
)

# Simulate data
empty = mcp::mcp(segments, sample = FALSE)
ex_rel_prior = tibble::tibble(
  x = 1:100,
  y = empty$func_y(x, 25, 40, 1, -3, 0.5, 30, 40, 7)
)

# Save to mcp
usethis::use_data(ex_rel_prior, overwrite = TRUE)
