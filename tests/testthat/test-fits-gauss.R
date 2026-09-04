models_gauss = list(
  # Log-link Gaussian with a CP in the mean and weights.
  list(
    y | weights(w) ~ 1 + x,
    ~ 1 + x,
    simulated = list(
      cp_1 = 93,
      Intercept_1 = log(9),
      x_1 = 0.003,
      Intercept_2 = log(14),
      x_2 = -0.004,
      sigma_1 = 2
    ),
    newdata = data.frame(
      x = seq(1, 200, length.out = 400),
      w = rep(c(1, 2), 200)
    ),
    family = gaussian(link = "log"),
    chains = 2,
    warmup = 1000,
    iter = 2000,
    min_ess = 50,
    seed = 4)
)

apply_test_fit("Gaussian recovery", models_gauss)
