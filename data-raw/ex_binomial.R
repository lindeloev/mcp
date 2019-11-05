# Define the model
segments = list(
  y | trials(N) ~ 1,  # constant rate
  1 ~ 0 + x,  # joined changing rate
  1 ~ 1 + x  # disjoined changing rate
)

# Simulate data
empty = mcp::mcp(segments, family = binomial(), sample = FALSE)
ex_binomial = tibble::tibble(
  x = 1:100,
  N = sample(10, length(x), replace=TRUE),
  y = empty$func_y(x, N, 2, 0, -0.2, 0.05, 25, 65))

# Save to mcp
usethis::use_data(ex_binomial, overwrite = TRUE)
