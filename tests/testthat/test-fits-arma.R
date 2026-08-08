models_arma = list(
  # Segment-specific AR coefficients and an AR predictor effect.
  list(y ~ 1 + ar(2),
       ~ 0 + x + ar(1, 1 + x),
       ~ 0,
       simulated = list(
         cp_1 = 80,
         cp_2 = 140,
         Intercept_1 = -20,
         sigma_1 = 5,
         ar1_1 = 0.7,
         ar2_1 = -0.4,
         x_2 = 0.5,
         ar1_2 = 0.5,
         ar1_x_2 = -0.005
       ),
       chains = 2,
       adapt = 2000,
       iter = 2000,
       min_ess = 50)
)

apply_test_fit("ARMA (gauss) fit", models_arma)
