models_families_binomial = list(
  list(
    y | trials(N) ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = -0.4, x_1 = 0.006, Intercept_2 = 0.6, x_2 = -0.004),
    newdata = data.frame(x = seq(1, 200, length.out = 400), y = 0, N = 10L),
    family = binomial(link = "logit")
  ),
  list(
    y | trials(N) ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = -0.3, x_1 = 0.005, Intercept_2 = 0.4, x_2 = -0.003),
    newdata = data.frame(x = seq(1, 200, length.out = 400), y = 0, N = 10L),
    family = binomial(link = "probit")
  ),
  list(
    y | trials(N) ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = 0.4, x_1 = 0.001, Intercept_2 = 0.6, x_2 = -0.001),
    newdata = data.frame(x = seq(1, 200, length.out = 400), y = 0, N = 10L),
    family = binomial(link = "identity")
  )
)

apply_test_fit("Binomial family recovery", models_families_binomial)
