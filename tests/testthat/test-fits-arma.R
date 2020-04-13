models_arma = list(
  # Simple AR
  list(y ~ 1 + ar(1),
       simulated = list(
         int_1 = 30,
         ar1_1 = 0.7,
         sigma_1 = 10
       )),

  # Larger AR
  list(y ~ 1 + ar(2),
       ~ 0 + x + ar(1, 1 + x),
       ~ 0,
       simulated = list(
         cp_1 = 80,
         cp_2 = 140,
         int_1 = -20,
         sigma_1 = 5,
         ar1_1 = 0.7,
         ar2_1 = -0.4,
         x_2 = 0.5,
         ar1_2 = 0.5,
         ar1_x_2 = -0.02
       ))
)

apply_test_fit(models_arma, "ARMA (gauss) fit")
