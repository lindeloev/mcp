models_families_gaussian = list(
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
    family = gaussian(link = "log")
  )
)

apply_test_fit("Gaussian link recovery", models_families_gaussian)
