models_ar = list(
  # Segment-specific AR coefficients and an AR predictor effect.
  list(y ~ 1 + ar(2, series = id),
       ~ 0 + x + ar(1, 1 + x),
       ~ 0,
       simulated = list(
         cp_1 = 79.5,
         cp_2 = 139.5,
         Intercept_1 = -20,
         sigma_1 = 5,
         ar1_1 = 0.7,
         ar2_1 = -0.4,
         x_2 = 0.5,
         ar1_2 = 0.5,
         ar1_x_2 = -0.002
       ),
       newdata = data.frame(
         id = rep(c("a", "b"), each = 200),
         x = rep(seq(1, 200, length.out = 200), 2),
         y = 0
       ),
       chains = 2,
       warmup = 2000,
       iter = 2000,
       min_ess = 50)
)

apply_test_fit("AR (Gaussian) fit", models_ar)
