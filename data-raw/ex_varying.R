# Define the model
model = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope
)

# Simulate data
empty = mcp::mcp(model, sample=FALSE)
ex_varying = tibble::tibble(id = c("John", "Benny", "Rose", "Cath", "Bill", "Erin")) %>%
  tidyr::expand_grid(x = seq(1, 100, by=4)) %>%
  dplyr::mutate(
    id_numeric = as.numeric(as.factor(id)),
    y = empty$simulate(
      x,
      cp_1 = 40,
      cp_1_id = 7*(id_numeric - mean(id_numeric)),
      int_1 = 15,
      x_1 = 3,
      x_2 = -2,
      sigma = 25
    )
  )

# Save to mcp
ex_varying = data.frame(ex_varying)
usethis::use_data(ex_varying, overwrite = TRUE)
