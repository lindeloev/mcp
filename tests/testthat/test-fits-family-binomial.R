models_families_binomial = list(
  list(
    y | trials(N) ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = -0.3, x_1 = 0.005, Intercept_2 = 0.8, x_2 = -0.003),
    newdata = data.frame(x = seq(1, 200, length.out = 400), y = 0, N = 10L),
    family = binomial(link = "probit"),
    chains = 2,
    adapt = 500,
    iter = 750,
    min_ess = 50
  ),
  list(
    y | trials(N) ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = 0.4, x_1 = 0.001, Intercept_2 = 0.7, x_2 = -0.001),
    newdata = data.frame(x = seq(1, 200, length.out = 400), y = 0, N = 10L),
    family = binomial(link = "identity"),
    chains = 2,
    adapt = 500,
    iter = 1000,
    min_ess = 50
  )
)

apply_test_fit("Binomial family recovery", models_families_binomial)
