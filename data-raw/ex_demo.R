# Define model
segments = list(
  response ~ 1,
  ~ 0 + time,
  ~ 1 + time
)

# Simulate data
empty = mcp::mcp(segments, sample = FALSE)

set.seed(40)
ex_demo = tibble::tibble(
  time = runif(100, 0, 100),
  response = empty$simulate(
    time,
    cp_1 = 30,
    cp_2 = 70,
    int_1 = 10,
    int_3 = 0,
    sigma = 4,
    time_2 = 0.5,
    time_3 = -0.2)
)

# Save it to package
ex_demo = data.frame(ex_demo)
usethis::use_data(ex_demo, overwrite = TRUE)
