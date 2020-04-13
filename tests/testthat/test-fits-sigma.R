models_sigma = list(
  # Simple sigma
  list(y ~ 1 + sigma(1 + x),
       simulated = list(
         int_1 = 30,
         sigma_1 = 5,
         sigma_x_1 = 0.1
       )),

  # Larger sigma
  list(y ~ 1 + sigma(1),
       ~ 0 + x + sigma(1 + x),
       ~ 0,
       simulated = list(
         cp_1 = 80,
         cp_2 = 140,
         int_1 = -20,
         sigma_1 = 2,
         x_2 = 0.5,
         sigma_2 = 5,
         sigma_x_2 = 0.2
       ))
)

apply_test_fit(models_sigma, "Sigma (gauss) fit")
