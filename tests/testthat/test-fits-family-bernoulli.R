models_families_bernoulli = list(
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 100, Intercept_1 = -0.4, x_1 = 0.006, Intercept_2 = 2.2, x_2 = -0.004),  # A unimodal posterior (cp_1 = 93 had a small second mode)
    family = bernoulli(link = "logit"),
    chains = 2,
    warmup = 500,
    iter = 1000,
    min_ess = 50,
    seed = 42
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 100, Intercept_1 = -0.3, x_1 = 0.005, Intercept_2 = 1.4, x_2 = -0.003),  # cp_1 = 93 mixed poorly for some seeds
    family = bernoulli(link = "probit"),
    chains = 2,
    warmup = 500,
    iter = 1000,
    min_ess = 50
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 85, Intercept_1 = 0.2, x_1 = 0.001, Intercept_2 = 0.8, x_2 = -0.001),  # cp_1 = 93 and 100 mixed poorly for some seeds
    family = bernoulli(link = "identity"),
    chains = 3,
    warmup = 500,
    iter = 1000,
    min_ess = 50
  )
)

apply_test_fit("Bernoulli family recovery", models_families_bernoulli)
