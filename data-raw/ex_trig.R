# Define the model
segments_trig = list(
  y ~ 1 + sin(x),
  1 ~ 0 + cos(x) + x
)

# Simulate the data
empty_trig = mcp(segments_trig, sample = FALSE)
ex_trig = tibble::tibble(
  x = seq(0, 35, by = 0.15),
  y = empty_trig$func_y(x, 3, 10, 10, 8, 3, 16.5)
)

# Save to mcp
usethis::use_data(ex_trig, overwrite = TRUE)






