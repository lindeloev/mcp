########################################################################
# TEST THAT TRUE PARAMETERS ARE WITHIN 50% HDI OF ESTIMATED PARAMETERS #
########################################################################
library(dplyr)

# All relevant segments expressions
segments = list(
  y ~ 1 + x,
  1 ~ rel(1) + rel(x),
  rel(1) ~ rel(1) + rel(x),
  rel(1) ~ 0
)

# Simulation parameters
sim_x = runif(200, 0, 100)
func_args = list(
  x = sim_x,
  sigma = 5,
  int_1 = 10,
  int_2 = 20,
  int_3 = -40,
  x_1 = -2,
  x_2 = 3,
  x_3 = -2,
  cp_1 = 30,
  cp_2 = 20,
  cp_3 = 30
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
set.seed(42)
data = data.frame(
  x = sim_x,
  y = do.call(fit_empty$func_y, func_args)
)

# Fit model to simulated data. A pretty long run to ensure convergence
# and small MCMC error
fit = mcp(segments, data, n.adapt=2500, n.update=2500, n.iter=3000)

# Check: expect all estimates to be within 98% HDI
results_table = summary(fit, width = 0.95) %>%
  left_join(df, by="name") %>%
  mutate(score = theory > .lower & theory < .upper) %>

test_that("fit approximate default priors", {
  expect_true(
    all(results_table$score),
    info = mutate_if(results_table, is.numeric, round, digits=1))
})
