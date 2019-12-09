# Define model
segments = list(
  y ~ 1,
  ~ 0 + sigma(1 + x),
  ~ 0 + x
)

# Simulate data
empty = mcp::mcp(segments, sample = FALSE)

set.seed(30)
ex_variance = tibble::tibble(
  x = 1:100,
  y = empty$simulate(
    x,
    cp_1 = 25,
    cp_2 = 75,
    int_1 = 20,
    x_3 = 2,
    sigma_1 = 7,
    sigma_2 = 25,
    sigma_x_2 = -0.45)
  )


# Save it to package
ex_variance = data.frame(ex_variance)
usethis::use_data(ex_variance, overwrite = TRUE)
