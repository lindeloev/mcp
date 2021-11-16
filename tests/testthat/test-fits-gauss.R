models_gauss = list(
  # Simple
  list(y ~ 1,
       ~ 1,
       simulated = list(
         Intercept_1 = 10,
         Intercept_2 = 20,
         sigma_1 = 5,
         cp_1 = 100)),

  # A lot of terms
  list(y ~ 1 + x,
       ~ 0 + x,
       ~ 1 + sigma(1),
       simulated = list(
         cp_1 = 70,
         cp_2 = 70,
         Intercept_1 = 10,
         Intercept_3 = 0,
         x_1 = 0.5,
         x_2 = -1,
         sigma_1 = 3,
         sigma_3 = 6))
)

apply_test_fit("Gaussian fit", models_gauss)
