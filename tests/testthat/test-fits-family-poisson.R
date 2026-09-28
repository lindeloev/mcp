models_families_poisson = list(
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = log(3), x_1 = 0.004, Intercept_2 = log(7), x_2 = -0.003),
    family = poisson(link = "log"),
    chains = 2,
    warmup = 200,
    iter = 700,
    min_ess = 50
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = 2, x_1 = 0.01, Intercept_2 = 7, x_2 = -0.01),
    family = poisson(link = "identity"),
    chains = 2,
    warmup = 1000,
    iter = 1000,
    min_ess = 50,
    seed = 1
  ),
  list(
    y ~ 1 + x + offset(log(pop)),
    ~ 1 + x + offset(log(pop)),
    simulated = list(cp_1 = 93, Intercept_1 = log(0.03), x_1 = 0.004, Intercept_2 = log(0.07), x_2 = -0.003),
    newdata = data.frame(x = seq(1, 200, length.out = 400), pop = rep(c(10, 100), 200)),
    family = poisson(link = "log"),
    chains = 3,
    warmup = 1500,
    iter = 3000,
    min_ess = 50,
    seed = 1
  )
)

apply_test_fit("Poisson family recovery", models_families_poisson)
