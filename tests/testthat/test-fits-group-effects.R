group_ids = sprintf("id_%02d", 1:20)
group_fit_data = data.frame(
  x = seq(1, 200, length.out = 400),
  id = rep(group_ids, length.out = 400),
  y = 0
)
models_group_effects = list(
  # Group-specific residual SDs on the log scale.
  list(
    y ~ 1 + sigma(1 + (1 || id)),
    newdata = group_fit_data,
    simulated = list(
      Intercept_1 = 10,
      sigma_1 = log(2),
      sigma_1_id_sd = 0.35
    ),
    chains = 2,
    warmup = 1000,
    iter = 1500,
    min_ess = 50
  )
)

apply_test_fit("Predictor group-effect fit", models_group_effects)
