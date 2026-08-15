models_gauss = list(
  # Log-link Gaussian with a CP in the mean.
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(
      cp_1 = 93,
      Intercept_1 = log(9),
      x_1 = 0.003,
      Intercept_2 = log(14),
      x_2 = -0.004,
      sigma_1 = 2
    ),
    family = gaussian(link = "log"),
    chains = 2,
    warmup = 500,
    iter = 750,
    min_ess = 50)
)

apply_test_fit("Gaussian recovery", models_gauss)
