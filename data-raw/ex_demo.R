# Define model
segments = list(
  response ~ 1,
  1 ~ 0 + time,
  1 ~ 1 + time
)
empty = mcp::mcp(segments, sample = FALSE)

# Simulate data
set.seed(40)
ex_demo = tibble::tibble(
  time = runif(100, 0, 100),
  response = empty$func_y(time, 10, 0, 0.5, - 0.2, 30, 70, 4)
)

# Save it to package
usethis::use_data(ex_demo, overwrite = TRUE)
