models_families_negbinomial = list(
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = log(4), x_1 = 0.003, Intercept_2 = log(8), x_2 = -0.002, shape_1 = log(10)),
    family = negbinomial(),
    chains = 2,
    warmup = 500,
    iter = 750,
    min_ess = 50
  )
)

apply_test_fit("Negative binomial family recovery", models_families_negbinomial)
