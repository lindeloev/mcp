models_families_bernoulli = list(
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = -0.4, x_1 = 0.006, Intercept_2 = 2, x_2 = -0.004),
    family = bernoulli(link = "logit"),
    chains = 2,
    warmup = 1000,
    iter = 2000,
    min_ess = 50
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = -0.3, x_1 = 0.005, Intercept_2 = 1.4, x_2 = -0.003),
    family = bernoulli(link = "probit"),
    chains = 2,
    warmup = 500,
    iter = 1000,
    min_ess = 50
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = 0.4, x_1 = 0.001, Intercept_2 = 0.8, x_2 = -0.001),
    family = bernoulli(link = "identity"),
    chains = 2,
    warmup = 500,
    iter = 1000,
    min_ess = 50
  )
)

apply_test_fit("Bernoulli family recovery", models_families_bernoulli)
