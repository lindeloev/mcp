# Define model
model = list(
  y ~ 1,
  ~ 1
)
empty = mcp::mcp(model, sample = FALSE, par_x = "x")

# Simulate data
set.seed(40)
ex_plateaus = tibble::tibble(
  x = runif(100, 0, 100),
  y = empty$simulate(
    x,
    cp_1 = 50,
    int_1 = 10,
    int_2 = 20,
    sigma = 8)
)

# Save to mcp
ex_plateaus = data.frame(ex_plateaus)
usethis::use_data(ex_plateaus, overwrite = TRUE)
