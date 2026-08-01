models_families_poisson = list(
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = log(3), x_1 = 0.004, Intercept_2 = log(5), x_2 = -0.003),
    family = poisson(link = "log")
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = 3, x_1 = 0.01, Intercept_2 = 5, x_2 = -0.01),
    family = poisson(link = "identity")
  )
)

apply_test_fit("Poisson family recovery", models_families_poisson)
