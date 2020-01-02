# Define the model
segments = list(
  y ~ 1 + sin(x),
  ~ 1 + cos(x) + x
)

# Simulate the data
empty = mcp::mcp(segments, sample = FALSE)
set.seed(40)
ex_trig = tibble::tibble(
  x = seq(0, 35, by = 0.2),
  y = empty$simulate(
    x,
    cp_1 = 17,
    int_1 = 10,
    x_1_sin = 10,
    x_2_cos = 8,
    int_2 = 10,
    x_2 = 3,
    sigma = 3
  )
)

# Save to mcp
ex_trig = data.frame(ex_trig)
usethis::use_data(ex_trig, overwrite = TRUE)
