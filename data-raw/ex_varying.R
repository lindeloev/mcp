# Define the model
segments = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope
)

# Simulate data
empty = mcp::mcp(segments, sample=FALSE)
ex_varying = tibble::tibble(id = c("John", "Benny", "Rose", "Cath", "Bill", "Erin")) %>%
  tidyr::expand_grid(x = seq(1, 100, by=4)) %>%
  dplyr::mutate(
    id_numeric = as.numeric(as.factor(id)),
    y = empty$func_y(x, 15, 3, -2, 40, 25, 7*(id_numeric - mean(id_numeric)))
  )

# Save to mcp
usethis::use_data(ex_varying, overwrite = TRUE)
