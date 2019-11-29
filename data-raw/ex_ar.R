# The model
segments = list(
  price ~ 1 + ar(2),
  ~ 0 + time
)

# Simulate fitted data
empty = mcp::mcp(segments, sample = FALSE)
ex_ar = tibble::tibble(
  time = 1:300,
  price = empty$simulate(
    time,
    cp_1 = 180,
    int_1 = 20,
    time_2 = 0.5,
    sigma_1 = NA,
    type = "fitted")
)

# Add residuals
set.seed(42)
residuals = stats::arima.sim(list(ar = c(0.5, 0.3)), nrow(ex_ar), sd = 5)
ex_ar$price = ex_ar$price + c(residuals)

# Save it
ex_ar = data.frame(ex_ar)
usethis::use_data(ex_ar, overwrite = TRUE)
