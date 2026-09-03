#' Run example models
#'
#' @aliases mcp_example
#' @param name Name of the example. One of:
#'  * `"demo"`: Two change points between intercepts and joined/disjoined slopes.
#'  * `"intercepts"`: An intercept-only change point.
#'  * `"multiple"`: Multiple regression with categorical predictors and interactions.
#'  * `"binomial"`: Binomial with two change points. Much like `"demo"` on a logit scale.
#'  * `"group_mu"`: Group-level predictor deviations (random intercepts/slopes) across a change point.
#'  * `"group_cp"`: Group-level change-point deviations (random effects).
#'  * `"missing"`: Missing data imputation (NAs in response variable y).
#'  * `"quadratic"`: A change point to a quadratic segment where there is no data.
#'  * `"ar"`: One change point in autoregressive residuals (the `ar1` dpar)
#'  * `"sigma"`: A change in "sigma" dpar, including a slope on sigma.
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
#' fit2 = mcp(empty$model, empty$data, family = empty$family)
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
  ~ 0 + time + ar(1)
)

# Simulate data
set.seed(42)
data = data.frame(
  time = 1:120,
  price = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, sample = FALSE)
data$price = empty$simulate(empty, data,
  cp_1 = 74.5,
  Intercept_1 = 20,
  time_2 = 0.5,
  sigma_1 = 5,
  ar1_1 = 0.4,
  ar2_1 = 0.15,
  ar1_2 = -0.3
)

# Run sampling
fit = mcp(model, data, sample = sample, warmup = 2000, iter = 5000, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  gg1 = plot(fit) + ggplot2::labs(title = \"plot(fit)\")
  set.seed(43)
  gg2 = plot_dpar(fit, \"ar1\") + ggplot2::labs(title = 'plot_dpar(fit, \"ar1\")')
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
  N = base::sample(15, 25, replace=TRUE),
  y = 0.  # Numeric placeholder that is valid for every sampled trial count.
)
empty = mcp(model, data, family = binomial(), sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 29.5,
  cp_2 = 69.5,
  Intercept_1 = 1.5,
  Intercept_3 = -1,
  x_2 = -0.15,
  x_3 = 0.05
)

# Run sampling
fit = mcp(model, data, family = binomial(), iter = 4000, sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, q_fit = TRUE) + ggplot2::labs(title = 'plot(fit, q_fit = TRUE)'))
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
  Intercept_3 = 20,
  time_3 = -0.3,
  sigma_1 = 4
)

# Run sampling
fit = mcp(model, data, sample = sample, seed = 78)

# Illustrative plot
if (plot) {
  set.seed(78)
  print(plot(fit, q_fit = TRUE) + ggplot2::labs(title = 'plot(fit, q_fit = TRUE)'))
}",


    group_mu = "# Define model
model = list(
  y ~ 1 + condition + (condition || id),
  ~ 1 + condition  # group effects carry into this segment
)

# Simulate balanced data with 9 levels of the grouping factor
set.seed(200)
data = tidyr::expand_grid(
  id = sprintf('id_%02d', 1:9),
  condition = factor(c('A', 'B')),
  x = seq(0, 100, length.out = 9)
)
data$y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.

empty = mcp(model, data, par_x = 'x', sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 50,
  Intercept_1 = 10,
  conditionB_1 = 4,
  Intercept_2 = 13,
  conditionB_2 = -2,
  Intercept_1_id_sd = 2,
  conditionB_1_id_sd = 2,
  sigma_1 = 1.5
)

# Run sampling
fit = mcp(model, data, par_x = 'x', sample = sample, iter = 15000, seed = 200)

# Illustrative plot
if (plot) {
  set.seed(200)
  print(plot(fit, facet_by = 'id', color_by = 'condition') +
      ggplot2::labs(title = 'plot(fit, facet_by = \"id\", color_by = \"condition\")'))
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
  sigma_1 = 8
)

# Run sampling
fit = mcp(model, data, par_x = 'x', sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, q_fit = TRUE) + ggplot2::labs(title = 'plot(fit, q_fit = TRUE)'))
}",
    missing = "# Define model
model = list(
  y ~ 1 + x + condition,
  ~ 0 + x  # condition effect carries into the second segment
)

# Simulate complete data
set.seed(42)
data = data.frame(
  x = 1:100,
  condition = factor(rep(c('A', 'B'), 50)),
  y = 2.  # Numeric placeholder replaced by simulation below.
)
empty = mcp(model, data, par_x = 'x', sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 54.5,
  Intercept_1 = 10,
  x_1 = 0.25,
  conditionB_1 = 22,
  x_2 = -0.4,
  sigma_1 = 4
)

# Remove scattered responses and a short run away from the change point
missing_rows = c(8, 19, 27:31, 68, 84, 96)
data$y[missing_rows] = NA

# Run sampling; JAGS retains posterior draws for the missing responses
# See them using fitted(fit, newdata = fit$data[is.na(fit$data$y), ], summary = FALSE)
fit = mcp(model, data, par_x = 'x', iter = 5000, sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, q_fit = TRUE, color_by = 'condition') +
      ggplot2::labs(title = 'plot(fit, q_fit = TRUE, color_by = \"condition\")'))
}",
    multiple = "# Define model
model = list(
  y ~ 1 + x:group + z,
  ~ 1 + x + group,
  ~ 0 + I(x^2)
)

# Simulate data
set.seed(42)
data = data.frame(
  x = 1:120,
  group = rep(c('A', 'B', 'C', 'D'), 30),
  z = rnorm(120, mean = 1:120, sd = 25),
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, sample = FALSE, par_x = 'x')
data$y = empty$simulate(empty, data,
  cp_1 = 69.5,
  cp_2 = 99.5,

  Intercept_1 = 10,
  z_1 = 0.2,
  xgroupA_1 = -0.75,
  xgroupB_1 = -0.25,
  xgroupC_1 = 0.25,
  xgroupD_1 = 0.75,

  Intercept_2 = 10,
  x_2 = -1,
  groupB_2 = 15,
  groupC_2 = 30,
  groupD_2 = 45,

  xE2_3 = 0.2,

  sigma_1 = 5
)

# Run sampling
fit = mcp(model, data, par_x = 'x', iter = 10000, sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, color_by = 'group') + ggplot2::labs(title = 'plot(fit, color_by = \"group\")'))
}",
    quadratic = "# Define model
model = list(
  y ~ 1,
  ~ 0 + x + I(x^2)
)

# Simulate data
set.seed(40)
data = data.frame(
  x = c(seq(0, 10, by = 0.3), seq(20, 40, by = 0.3)),
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp::mcp(model, data, sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 15,
  Intercept_1 = 10,
  x_2 = -30,
  xE2_2 = 1.5,
  sigma_1 = 20
)

# Run sampling
fit = mcp(model, data, warmup = 2000, iter = 10000, sample = sample, diagnostics = FALSE, seed = 6)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, q_fit  = TRUE, q_predict = TRUE) + ggplot2::labs(title = 'plot(fit, q_fit  = TRUE, q_predict = TRUE)'))
}",
    sigma = "# Define model
model = list(
  y ~ 1,
  ~ 0 + sigma(1 + x),
  ~ 0 + x
)


# Simulate data
set.seed(40)
data = data.frame(
  x = 1:100,
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp::mcp(model, data, sample = FALSE)
data$y = empty$simulate(empty, data,
    cp_1 = 24.5,
    cp_2 = 75,
    Intercept_1 = 20,
    x_3 = 1,
    sigma_1 = log(7),
    sigma_2 = log(25),
    sigma_x_2 = log(2.5 / 25) / (75 - 24.5)
  )

# Run sampling
fit = mcp(model, data, iter = 3000, sample = sample, seed = 40)

# Illustrative plot
if (plot) {
  set.seed(40)
  gg1 = plot(fit, q_predict = TRUE) + ggplot2::labs(title = 'plot(fit, q_predict = TRUE)')
  set.seed(41)
  gg2 = plot_dpar(fit, 'sigma') + ggplot2::labs(title = 'plot_dpar(fit, \"sigma\")')
  print(gg1 / gg2)
}",
    group_cp = "# Define model
model = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope
)

# Simulate data
set.seed(42)
data = tibble::tibble(id = c('John', 'Omar', 'Rose', 'Cath', 'Ni', 'Erin', 'Frank', 'Mark', 'Slava'))
data = tidyr::expand_grid(data, x = seq(1, 100, by=20))
data$y = 2.# or whatever signals 'numeric'. Will be replaced by simulation below.

empty = mcp(model, data, sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 40,
  cp_1_sd = 20,
  Intercept_1 = 15,
  x_1 = 3,
  x_2 = -2,
  sigma_1 = 10
)

# Run sampling
fit = mcp(model, data, iter = 4500, sample = sample, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, facet_by = 'id') + ggplot2::labs(title = 'plot(fit, facet_by = \"id\")'))
}"
  )

  # Run the code
  name = rlang::arg_match0(name, names(examples))
  eval(str2expression(examples[[name]]))

  fit$example_code = examples[[name]]
  class(fit$example_code) = c("mcptext", "character")
  fit
}


#' @aliases mcp_example_data
#' @export
#' @describeIn mcp_example Conveniently get simulated data only.
mcp_example_data = function(name) {
  mcp_example(name, sample = "none", plot = FALSE)$data
}
