group_ids = sprintf("id_%02d", 1:20)
group_fit_data = data.frame(
  x = seq(1, 200, length.out = 400),
  id = rep(group_ids, length.out = 400),
  condition = factor(rep(c("A", "B"), each = 20, times = 10)),
  y = 0
)
group_index = match(group_fit_data$id, group_ids)
mean_deviations = 3 * stats::qnorm(((seq_along(group_ids) - 0.5) / length(group_ids)))
sigma_deviations = 0.35 * stats::qnorm(((seq_along(group_ids) - 0.5) / length(group_ids)))
factor_a_deviations = 2 * stats::qnorm(((seq_along(group_ids) - 0.5) / length(group_ids)))
factor_b_deviations = 1.5 * stats::qnorm(
  ((c(11:20, 1:10) - 0.5) / length(group_ids))
)

models_group_effects = list(
  # A mean-intercept deviation that carries into the second segment.
  list(
    y ~ 1 + (1 | id),
    ~ 0 + x,
    newdata = group_fit_data,
    simulated = list(
      cp_1 = 100,
      Intercept_1 = 10,
      x_2 = 0.15,
      sigma_1 = 2,
      Intercept_1_id = mean_deviations[group_index]
    ),
    hyperparameters = list(Intercept_1_id_sd = 3)
  ),

  # Group-specific residual SDs on the log scale.
  list(
    y ~ 1 + sigma(1 + (1 || id)),
    newdata = group_fit_data,
    simulated = list(
      Intercept_1 = 10,
      sigma_1 = log(2),
      sigma_1_id = sigma_deviations[group_index]
    ),
    hyperparameters = list(sigma_1_id_sd = 0.35)
  ),

  # Independent group coefficients for both levels of a no-intercept factor.
  list(
    y ~ 0 + (0 + condition || id),
    newdata = group_fit_data,
    simulated = list(
      sigma_1 = 1.5,
      conditionA_1_id = factor_a_deviations[group_index],
      conditionB_1_id = factor_b_deviations[group_index]
    ),
    hyperparameters = list(
      conditionA_1_id_sd = 2,
      conditionB_1_id_sd = 1.5
    )
  )
)

apply_test_fit("Predictor group-effect fit", models_group_effects)
