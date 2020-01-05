# Define the model
model = list(
  y ~ 1,
  ~ 0 + x + I(x^2)
)

# Simulate the data
empty = mcp::mcp(model, sample = FALSE)
set.seed(42)
ex_quadratic = tibble::tibble(
  x = seq(0, 40, by = 0.5),
  y = empty$simulate(
    x,
    cp_1 = 15,
    int_1 = 10,
    sigma = 30,
    x_2 = -30,
    x_2_E2 = 1.5)
)

# Save to mcp
ex_quadratic = data.frame(ex_quadratic)
usethis::use_data(ex_quadratic, overwrite = TRUE)
