# Define the model
segments = list(
  y ~ 1 + sin(x),
  1 ~ 0 + cos(x) + x
)

# Simulate the data
empty = mcp(segments, sample = FALSE)
ex_trig = tibble::tibble(
  x = seq(0, 35, by = 0.15),
  y = empty$func_y(x, cos_x_2 = 8, cp_1 = 16.5, int_1 = 10, sigma = 3, sin_x_1 = 10, x_2 = 3)
)

# Save to mcp
usethis::use_data(ex_trig, overwrite = TRUE)
