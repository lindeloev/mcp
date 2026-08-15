models_families_poisson = list(
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = log(3), x_1 = 0.004, Intercept_2 = log(7), x_2 = -0.003),
    family = poisson(link = "log"),
    chains = 2,
    warmup = 500,
    iter = 750,
    min_ess = 50
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = 3, x_1 = 0.01, Intercept_2 = 6, x_2 = -0.01),
    family = poisson(link = "identity"),
    chains = 2,
    warmup = 500,
    iter = 500,
    min_ess = 50
  )
)

apply_test_fit("Poisson family recovery", models_families_poisson)
