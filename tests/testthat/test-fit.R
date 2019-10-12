########################################################################
# TEST THAT TRUE PARAMETERS ARE WITHIN 50% HDI OF ESTIMATED PARAMETERS #
########################################################################
# All relevant segments expressions
segments = list(
  y ~ 1 + x,
  1 ~ rel(1) + rel(x),
  1 ~ 0,
  1 ~ x
)

# Simulation parameters
sim_x = runif(200, 0, 100)
func_args = list(
  x = sim_x,
  sigma = 5,
  int_1 = 10,
  x_1 = -2,
  int_2 = 20,
  x_2 = 1,
  x_4 = 3,
  cp_1 = 30,
  cp_2 = 50,
  cp_3 = 80
)

# Make it a data.frame to use for joining with summary.mcpfit
df = func_args
df[["x"]] = NULL
df = data.frame(
  name = names(df),
  theory = as.numeric(df)
)
df$name = as.character(df$name)

# Simulate data
fit_empty = mcp(segments, sample=F)
data = data.frame(
  x = sim_x,
  y = do.call(fit_empty$func_y, func_args)
)

# Fit model to simulated data. A pretty long run to ensure convergence
# and small MCMC error
fit = mcp(segments, data, n.adapt=2500, n.update=2500, n.iter=3000)

# Check: expect all estimates to be within 50% HDI
X = summary(fit, width = 0.98) %>%
  left_join(df, by="name") %>%
  mutate(score = theory > .lower & theory < .upper)

test_that("fit approximate default priors", {
  expect_true(
    all(X$score),
    info = round(X, 3))
})
