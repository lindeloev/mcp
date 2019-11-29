# Define the model
segments = list(
  y | trials(N) ~ 1,  # constant rate
  ~ 0 + x,  # joined changing rate
  ~ 1 + x  # disjoined changing rate
)

# Simulate data
empty = mcp::mcp(segments, family = binomial(), sample = FALSE)
ex_binomial = tibble::tibble(
  x = 1:100,
  N = sample(10, length(x), replace=TRUE),
  y = empty$simulate(
    x,
    N,
    cp_1 = 35,
    cp_2 = 70,
    int_1 = 2,
    int_3 = 0,
    x_2 = -0.25,
    x_3 = 0.05)
  )

# Save to mcp
ex_binomial = data.frame(ex_binomial)
usethis::use_data(ex_binomial, overwrite = TRUE)
