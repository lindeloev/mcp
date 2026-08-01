#' Get example models and data
#'
#' @aliases mcp_example
#' @param name Name of the example. One of:
#'  * `"demo"`: Two change points between intercepts and joined/disjoined slopes.
#'  * `"ar"`: One change point in autoregressive residuals.
#'  * `"binomial"`: Binomial with two change points. Much like `"demo"` on a logit scale.
#'  * `"group"`: Group-level (random) intercepts and factor effects across a change point.
#'  * `"intercepts"`: An intercept-only change point.
#'  * `"multiple"`: Multiple regression with categorical predictors and interactions.
#'  * `"quadratic"`: A change point to a quadratic segment.
#'  * `"varying"`: Varying / hierarchical change points.
#'  * `"variance"`: A change in variance, including a variance slope.
#' @param plot Logical. Plot the fitted example? No plot is produced when
#'   `sample = FALSE`.
#' @inheritParams mcp
#' @return An `mcpfit`, enriched with a `$call` field. It contains the code to
#'   reproduce the data and the fit.
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' \donttest{
#' fit = mcp_example("multiple")
#' print(fit$call) # See how the data was simulated
#'
#' # Without sampling
#' empty = mcp_example("binomial", sample = FALSE, plot = FALSE)
#' print(empty)
#' print(empty$call)
#' print(empty$data)
#' #'
#' #' # Now sample this model
#' fit2 = mcp(empty$model, empty$data, family = empty$family)
#' plot(fit2)
#' }
mcp_example = function(name, sample = "post", warn = FALSE, plot = TRUE) {
  checkmate::assert_string(name)
  checkmate::assert_flag(warn)
  checkmate::assert_flag(plot)
  data = data.frame() # To make R CMD Check happy.
  fit = NULL # To make R CMD Check happy.
  plot = plot && !(sample %in% c(FALSE, "none"))

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
  cp_1 = 75,
  Intercept_1 = 20,
  time_2 = 0.5,
  sigma_1 = 5,
  ar1_1 = 0.6,
  ar2_1 = 0.3,
  ar1_2 = -0.5
)

# Run sampling
fit = mcp(model, data, sample = sample, adapt = 2000, iter = 5000, warn = warn, seed = 42)

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
  y | trials(N) ~ 1,  # constant rate
  ~ 0 + x,  # joined changing rate
  ~ 1 + x  # disjoined changing rate
)

# Simulate data
set.seed(42)
data = data.frame(
  x = 1:100,
  N = base::sample(10, 100, replace=TRUE),
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
fit = mcp(model, data, family = binomial(), sample = sample, warn = warn, seed = 42)

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
set.seed(42)
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
  Intercept_3 = 0,
  time_3 = -0.2,
  sigma_1 = 4
)

# Run sampling
fit = mcp(model, data, iter = 4000, sample = sample, warn = warn, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit) + ggplot2::labs(title = 'plot(fit)'))
}",
    group = "# Define model
model = list(
  y ~ 1 + condition + (condition || participant),
  ~ 1 + condition  # group effects carry into this segment
)

# Simulate balanced data with 12 levels of the grouping factor
set.seed(200)
data = tidyr::expand_grid(
  participant = sprintf('participant_%02d', 1:12),
  condition = factor(c('A', 'B')),
  x = seq(0, 100, length.out = 12)
)
data$y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.

empty = mcp(model, data, par_x = 'x', sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 50,
  Intercept_1 = 10,
  conditionB_1 = 4,
  Intercept_2 = 15,
  conditionB_2 = -2,
  Intercept_1_participant_sd = 2,
  conditionB_1_participant_sd = 1.5,
  sigma_1 = 1.5
)

# Run sampling
fit = mcp(
  model, data, par_x = 'x', sample = sample, warn = warn,
  adapt = 2000, iter = 20000, seed = 200
)

# Illustrative plot
if (plot) {
  set.seed(200)
  print(plot(fit, facet_by = 'participant', color_by = 'condition') +
      ggplot2::labs(title = 'plot(fit, facet_by = \"participant\", color_by = \"condition\")'))
}",
    intercepts = "# Define model
model = list(
  y ~ 1,
  ~ 1
)

# Simulate data
set.seed(42)
data = data.frame(
  x = runif(100, 0, 100),
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
fit = mcp(model, data, par_x = 'x', sample = sample, warn = warn, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit) + ggplot2::labs(title = 'plot(fit)'))
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
  cp_1 = 70,
  cp_2 = 100,

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
fit = mcp(model, data, par_x = 'x', sample = sample, warn = warn, seed = 42)

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
set.seed(42)
data = data.frame(
  x = seq(0, 40, by = 0.5),
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp::mcp(model, data, sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 15,
  Intercept_1 = 10,
  x_2 = -30,
  xE2_2 = 1.5,
  sigma_1 = 30
)

# Run sampling
fit = mcp(model, data, sample = sample, warn = warn, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit) + ggplot2::labs(title = 'plot(fit)'))
}",
    variance = "# Define model
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
    x_3 = 2,
    sigma_1 = log(7),
    sigma_2 = log(25),
    sigma_x_2 = log(2.5 / 25) / (75 - 24.5)
  )

# Run sampling
fit = mcp(model, data, iter = 4000, adapt = 3000, sample = sample, warn = warn, seed = 40)

# Illustrative plot
if (plot) {
  set.seed(40)
  gg1 = plot(fit, q_predict = TRUE) + ggplot2::labs(title = 'plot(fit, q_predict = TRUE)')
  set.seed(41)
  gg2 = plot_dpar(fit, 'sigma') + ggplot2::labs(title = 'plot_dpar(fit, \"sigma\")')
  print(gg1 / gg2)
}",
    varying = "# Define model
model = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope
)

# Simulate data
set.seed(42)
data = tibble::tibble(id = c('John', 'Benny', 'Rose', 'Cath', 'Bill', 'Erin'))
data = tidyr::expand_grid(data, x = seq(1, 100, by=4))
data$id_numeric = as.numeric(as.factor(data$id))
data$y = 2.# or whatever signals 'numeric'. Will be replaced by simulation below.

empty = mcp(model, data, sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 40,
  cp_1_id = 7*(data$id_numeric - mean(data$id_numeric)),
  Intercept_1 = 15,
  x_1 = 3,
  x_2 = -2,
  sigma_1 = 25
)

# Run sampling
fit = mcp(model, data, sample = sample, warn = warn, seed = 42)

# Illustrative plot
if (plot) {
  set.seed(42)
  print(plot(fit, facet_by = 'id') + ggplot2::labs(title = 'plot(fit, facet_by = \"id\")'))
}"
  )

  # Run the code
  name = rlang::arg_match0(name, names(examples))
  eval(str2expression(examples[[name]]))

  fit$call = examples[[name]]
  class(fit$call) = c("mcptext", "character")
  fit
}


#' @aliases mcp_example_data
#' @export
#' @describeIn mcp_example Conveniently get simulated data only.
mcp_example_data = function(name) {
  mcp_example(name, sample = FALSE, plot = FALSE)$data
}
