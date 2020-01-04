# The model
model = list(
  price ~ 1 + ar(2),
  ~ 0 + time + ar(1)
)

# Simulate fitted data
empty = mcp::mcp(model, sample = FALSE)
set.seed(42)
ex_ar = tibble::tibble(
  time = 1:200,
  price = empty$simulate(
    time,
    cp_1 = 120,
    int_1 = 20,
    time_2 = 0.5,
    sigma_1 = 5,
    ar1_1 = 0.7,
    ar2_1 = 0.2,
    ar1_2 = -0.4)
)

# Save it
ex_ar = data.frame(ex_ar)
usethis::use_data(ex_ar, overwrite = TRUE)
