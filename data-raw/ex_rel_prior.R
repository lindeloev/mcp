# Define the model
segments = list(
  y ~ 1 + x,
  ~ rel(1) + rel(x),
  rel(1) ~ 0 + x
)

# Simulate data
empty = mcp::mcp(segments, sample = FALSE)
set.seed(40)
ex_rel_prior = tibble::tibble(
  x = 1:100,
  y = empty$simulate(
    x,
    cp_1 = 25,
    cp_2 = 40,
    int_1 = 25,
    int_2 = -10,
    sigma = 7,
    x_1 = 1,
    x_2 = -3,
    x_3 = 0.5)
)

# Save to mcp
ex_rel_prior = data.frame(ex_rel_prior)
usethis::use_data(ex_rel_prior, overwrite = TRUE)
