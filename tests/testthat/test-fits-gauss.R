models_gauss = list(
  # Simple
  list(y ~ 1,
       ~ 1,
       simulated = list(
         int_1 = 10,
         int_2 = 20,
         sigma_1 = 5,
         cp_1 = 100)),

  # A lot of terms
  list(y ~ 1 + x + sin(x),
       ~ rel(1) + rel(x),
       ~ 0,
       simulated = list(
         cp_1 = 70,
         cp_2 = 140,
         int_1 = 10,
         x_1 = 0.5,
         x_2 = -1,
         x_1_sin = 5,
         sigma_1 = 3,
         int_2 = -50))
)

apply_test_fit(models_gauss, "Gaussian fit")
