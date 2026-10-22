#' Run example models
#'
#' @aliases mcp_example
#' @param name Name of the example. One of:
#'  * `"ar"`: One change point in autoregressive residuals (the `ar1` dpar)
#'  * `"binomial"`: Binomial with two change points. Much like `"demo"` on a logit scale.
#'  * `"demo"`: Two change points between intercepts and joined/disjoined slopes.
#'  * `"group_mu"`: Group-level predictor deviations (random intercepts/slopes) across a change point.
#'  * `"group_cp"`: Group-level change-point deviations (random effects).
#'  * `"intercepts"`: An intercept-only change point.
#'  * `"missing"`: Missing data imputation (NAs in response variable y).
#'  * `"multiple"`: Multiple regression with categorical predictors and interactions.
#'  * `"sigma"`: A change in "sigma" dpar, including a slope on sigma.
#'  * `"void"`: Missing data around change points.
#' @param plot Logical. Plot the fitted example? Requires sample != "none".
#' @inheritParams mcp
#' @return An `mcpfit`, enriched with an `$example_code` field. It contains the code to
#'   reproduce the data and the fit.
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' \donttest{
#' fit = mcp_example("multiple")
#' print(fit$example_code) # See how the data was simulated
#'
#' # Without sampling
#' empty = mcp_example("binomial", sample = "none", plot = FALSE)
#' print(empty)
#' print(empty$example_code)
#'
#' # Now sample this model
#' fit2 = mcp(empty$model, empty$data, family = empty$family, seed = 42)
#' plot(fit2)
#' }
mcp_example = function(name, sample = "post", plot = TRUE) {
  checkmate::assert_string(name)
  checkmate::assert_flag(plot)
  data = data.frame() # To make R CMD Check happy.
  fit = NULL # To make R CMD Check happy.
  plot = plot && !identical(sample, FALSE) && !identical(sample, "none")

  examples = list(
    ar = "# Define model
model = list(
  price ~ 1 + ar(2),
  ~ 0 + ar(1)
)

# Simulate data
set.seed(42)
data = data.frame(
  time = 1:150,
  price = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, par_x = 'time', sample = FALSE)
data$price = empty$simulate(empty, data,
  cp_1 = 90,
  Intercept_1 = 20,
  sigma_1 = 5,
  ar1_1 = 0.7,
  ar2_1 = 0.2,
  ar1_2 = -0.3
)

# Run sampling
fit = mcp(model, data, par_x = 'time', sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  gg1 = plot(fit) + ggplot2::labs(title = 'fit = mcp_example(\"ar\"); plot(fit)')
  set.seed(43)
  gg2 = plot_dpar(fit, \"ar1\") + ggplot2::labs(title = 'plot_dpar(fit, \"ar1\")')

  # Using patchwork: library(patchwork)
  print(gg1 / gg2)
}",


    binomial = "# Define model
model = list(
  y | trials(N) ~ 1,  # constant success probability
  ~ 0 + x,  # joined changing success probability
  ~ 1 + x  # disjoined changing success probability
)

# Simulate data
set.seed(42)
data = data.frame(
  x = 1:100,
  N = base::sample(15:25, 100, replace=TRUE),
  y = 0.  # Numeric placeholder that is valid for every sampled trial count.
)
empty = mcp(model, data, family = binomial(), sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 30,
  cp_2 = 70,
  Intercept_1 = 1.5,
  Intercept_3 = -1,
  x_2 = -0.15,
  x_3 = 0.05
)

# Run sampling
fit = mcp(model, data, family = binomial(), sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, interval = \"fitted\") + ggplot2::labs(title = 'fit = mcp_example(\"binomial\"); plot(fit, interval = \"fitted\")'))
}",


    demo = "# Define model
model = list(
  response ~ 1,
  ~ 0 + time,
  ~ 1 + time
)

# Simulate data
set.seed(78)
data = data.frame(
  time = runif(100, 0, 100),
  response = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, sample = FALSE)
data$response = empty$simulate(empty, data,
  cp_1 = 30,
  cp_2 = 70,
  Intercept_1 = 10,
  time_2 = 0.5,
  Intercept_3 = 15,
  time_3 = -0.3,
  sigma_1 = 3.5
)

# Run sampling
fit = mcp(model, data, sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(78)
  print(plot(fit, interval = \"fitted\") + ggplot2::labs(title = 'fit = mcp_example(\"demo\"); plot(fit, interval = \"fitted\")'))
}",


    group_mu = "# Define model
model = list(
  y ~ 1 + state + (state || id),
  ~ 1 + state + same((state || id))  # identical group effects in both segments
)

# Simulate balanced data with 9 levels of the grouping factor
set.seed(42)
data = tidyr::expand_grid(
  id = sprintf('id_%02d', 1:9),
  state = factor(c('A', 'B')),
  x = seq(0, 100, length.out = 9)
)
data$y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.

empty = mcp(model, data, par_x = 'x', sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 50,
  Intercept_1 = 10,
  stateB_1 = 4,
  Intercept_2 = 13,
  stateB_2 = -2,
  Intercept_1_id_sd = 2,
  stateB_1_id_sd = 2,
  sigma_1 = 1.5
)

# Run sampling
fit = mcp(model, data, par_x = 'x', sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, facet_by = 'id') + 
      ggplot2::labs(title = 'fit = mcp_example(\"group_mu\"); plot(fit, facet_by = \"id\")'))
}",
    intercepts = "# Define model
model = list(
  y ~ 1,
  ~ 1
)

# Simulate data
set.seed(42)
data = data.frame(
  x = runif(60, 0, 100),
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, sample = FALSE, par_x = 'x')
data$y = empty$simulate(empty, data,
  cp_1 = 50,
  Intercept_1 = 10,
  Intercept_2 = 20,
  sigma_1 = 7
)

# Run sampling
fit = mcp(model, data, par_x = 'x', sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, interval = \"fitted\") + ggplot2::labs(title = 'fit = mcp_example(\"intercepts\"); plot(fit, interval = \"fitted\")'))
}",
    missing = "# Define model
model = list(
  y ~ 1 + x + state,
  ~ 0 + x + same(state)  # reuse the state coefficients
)

# Simulate complete data
set.seed(42)
data = data.frame(
  x = 1:100,
  state = factor(rep(c('A', 'B'), 50)),
  y = 2.  # Numeric placeholder replaced by simulation below.
)
empty = mcp(model, data, par_x = 'x', sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 55,
  Intercept_1 = 10,
  x_1 = 0.25,
  stateB_1 = 22,
  x_2 = -0.4,
  sigma_1 = 4
)

# Remove scattered responses and a short run away from the change point
missing_rows = c(8, 19, 27:31, 68, 84, 96)
data$y[missing_rows] = NA

# Run sampling; JAGS imputes the missing responses jointly with the parameters
# Predictions at these rows follow the same distribution as the imputations (no AR/MA):
# predict(fit, newdata = fit$data[is.na(fit$data$y), ], summary = FALSE)
fit = mcp(model, data, par_x = 'x', sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, interval = \"fitted\") + 
      ggplot2::labs(title = 'fit = mcp_example(\"missing\"); plot(fit, interval = \"fitted\")'))
}",
    multiple = "# Define model
model = list(
  y ~ 1 + x:state + z,
  ~ 1 + x + state,
  ~ 0 + I(x^2) + same(state)
)

# Simulate data
set.seed(42)
data = data.frame(
  x = 1:132,
  state = rep(c('A', 'B', 'C', 'D'), 33),
  z = rnorm(132, mean = 1:132, sd = 25),
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, sample = FALSE, par_x = 'x')
data$y = empty$simulate(empty, data,
  cp_1 = 70,
  cp_2 = 100,

  Intercept_1 = 10,
  z_1 = 0.2,
  x.stateA_1 = -0.30,
  x.stateB_1 = -0.10,
  x.stateC_1 = 0.25,
  x.stateD_1 = 0.50,

  Intercept_2 = 15,
  x_2 = -0.8,
  stateB_2 = 10,
  stateC_2 = 20,
  stateD_2 = 30,

  xE2_3 = 0.02,

  sigma_1 = 5
)

# Run sampling
fit = mcp(model, data, par_x = 'x', sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit) + ggplot2::labs(title = 'fit = mcp_example(\"multiple\"); plot(fit)'))
}",
    sigma = "# Define model
model = list(
  y ~ 1,
  ~ 0 + sigma(1 + x),
  ~ 0 + x
)

# Simulate data
set.seed(42)
data = data.frame(
  x = 1:100,
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp::mcp(model, data, sample = FALSE)
data$y = empty$simulate(empty, data,
    cp_1 = 24.5,
    cp_2 = 75,
    Intercept_1 = 25,
    x_3 = 1,
    sigma_Intercept_1 = log(4),
    sigma_Intercept_2 = log(25),
    sigma_x_2 = log(3.5 / 25) / (75 - 24.5)
  )

# Run sampling
fit = mcp(model, data, sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  gg1 = plot(fit, interval = \"predict\") + ggplot2::labs(title = 'fit = mcp_example(\"sigma\"); plot(fit, interval = \"predict\")')
  set.seed(43)
  gg2 = plot_dpar(fit, 'sigma') + ggplot2::labs(title = 'plot_dpar(fit, \"sigma\")')

  # Using patchwork: library(patchwork)
  print(gg1 / gg2)
}",
    group_cp = "# Define model
model = list(
  y ~ 1 + x,              # intercept + slope
    ~ cp((1|id)) + 0 + x  # joined slope with group-level change points
)

# Simulate data
set.seed(42)
data = tibble::tibble(id = c('John', 'Omar', 'Rose', 'Cath', 'Ni', 'Erin', 'Frank', 'Mark', 'Slava'))
data = tidyr::expand_grid(data, x = seq(1, 100, by=20))
data$y = 2.# or whatever signals 'numeric'. Will be replaced by simulation below.

empty = mcp(model, data, sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 40,
  cp_1_id_sd = 20,
  Intercept_1 = 15,
  x_1 = 3,
  x_2 = -2,
  sigma_1 = 10
)

# Run sampling
fit = mcp(model, data, sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, facet_by = 'id') + ggplot2::labs(title = 'fit = mcp_example(\"group_cp\"); plot(fit, facet_by = \"id\")'))
}",

    void = "# demo_fit, with data removed near change points
data = demo_fit$data
data$response[data$time > 25 & data$time < 42] = NA
data$response[data$time > 70 & data$time < 80] = NA

# Run sampling
fit = mcp(demo_fit$model, data, sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  gg1 = plot(demo_fit, interval = \"fitted\") + ggplot2::labs(title = 'plot(demo_fit, interval = \"fitted\")')
  set.seed(43)
  gg2 = plot(fit, interval = \"fitted\") + ggplot2::labs(title = 'fit = mcp_example(\"void\"); plot(fit, interval = \"fitted\")')

  # Using patchwork: library(patchwork)
  print(gg1 / gg2)
}"
)

  # Run the code
  name = rlang::arg_match0(name, names(examples))
  eval(str2expression(examples[[name]]))

  # Prepend the arguments so the code runs standalone
  fit$example_code = paste0("# mcp_example() arguments\nsample = ", deparse(sample), "\nplot = ", plot, "\n\n", examples[[name]])
  class(fit$example_code) = c("mcptext", "character")
  fit
}


#' @aliases mcp_example_data
#' @export
#' @describeIn mcp_example Conveniently get simulated data only.
mcp_example_data = function(name) {
  mcp_example(name, sample = "none", plot = FALSE)$data
}
