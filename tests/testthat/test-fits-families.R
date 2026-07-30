models_families = list(
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(
      cp_1 = 93,
      Intercept_1 = log(9),
      x_1 = 0.003,
      Intercept_2 = log(14),
      x_2 = -0.004,
      sigma_1 = 2
    ),
    family = gaussian(link = "log")
  ),
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
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = -0.4, x_1 = 0.006, Intercept_2 = 0.6, x_2 = -0.004),
    family = bernoulli(link = "logit")
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = -0.3, x_1 = 0.005, Intercept_2 = 0.4, x_2 = -0.003),
    family = bernoulli(link = "probit")
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = 0.4, x_1 = 0.001, Intercept_2 = 0.6, x_2 = -0.001),
    family = bernoulli(link = "identity")
  ),
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
  ),
  list(
    y ~ 1 + x,
    ~ 1 + x,
    simulated = list(cp_1 = 93, Intercept_1 = log(4), x_1 = 0.003, Intercept_2 = log(6), x_2 = -0.002, shape_1 = log(5)),
    family = negbinomial()
  )
)

apply_test_fit("Family and link recovery", models_families)
