#' Get example models and data
#'
#' @aliases mcp_example
#' @param name Name of the example. One of:
#'  * `"demo"`: Two change points between intercepts and joined/disjoined slopes.
#'  * `"ar"`: One change point in autoregressive residuals.
#'  * `"binomial"`: Binomial with two change points. Much like `"demo"` on a logit scale.
#'  * `"intercepts"`: An intercept-only change point.
#'  * `"multiple"`: Multiple regression with categorical predictors and interactions.
#'  * `"quadratic"`: A change point to a quadratic segment.
#'  * `"trigonometric"`: Trigonometric/seasonal data and model.
#'  * `"varying"`: Varying / hierarchical change points.
#'  * `"variance"`: A change in variance, including a variance slope.
#' @inheritParams mcp
#' @return An `mcpfit`, enriched with a `$call` field. It contains the code to
#'   reproduce the data and the fit.
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' \donttest{
#' fit = mcp_example("multiple")
#' plot(fit)
#' print(fit$call)  # See how the data was simulated
#'
#' # Without sampling
#' empty = mcp_example("binomial", sample = FALSE)
#' print(empty)
#' print(empty$call)
#' print(empty$data)
#' #'
#' #' # Now sample this model
#' fit2 = mcp(empty$model, empty$data, family = empty$family)
#' plot(fit2)
#'}
mcp_example = function(name, sample = "post") {
  assert_types(name, "character", len = 1)
  data = data.frame()  # To make R CMD Check happy.

  examples = list(
    ar = "# Define model
model = list(
  price ~ 1 + ar(2),
  ~ 0 + time + ar(1)
)

# Simulate data
set.seed(42)
data = data.frame(
  time = 1:200,
  price = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, sample = FALSE)
data$price = empty$simulate(empty, data,
  cp_1 = 120,
  Intercept_1 = 20,
  time_2 = 0.5,
  sigma_1 = 5,
  ar1_1 = 0.7,
  ar2_1 = 0.2,
  ar1_2 = -0.4
)

# Run sampling
fit = mcp(model, data, sample = sample)",



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
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, family = binomial(), sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 30,
  cp_2 = 70,
  Intercept_1 = 2,
  Intercept_3 = 0.4,
  x_2 = -0.2,
  x_3 = 0.05
)

# Run sampling
fit = mcp(model, data, family = binomial(), adapt = 5000, sample = sample)
",



demo = "# Define model
model = list(
  response ~ 1,
  ~ 0 + time,
  ~ 1 + time
)

# Simulate data
set.seed(40)
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
fit = mcp(model, data, sample = sample)",



intercepts = "# Define model
model = list(
  y ~ 1,
  ~ 1
)

# Simulate data
set.seed(40)
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
fit = mcp(model, data, par_x = 'x', sample = sample)",


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
fit = mcp(model, data, par_x = 'x', cores = 3, sample = sample)",


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
fit = mcp(model, data, sample = sample)",



trigonometric = "model = list(
  y ~ 1 + sin(x),
  ~ 1 + cos(x) + x
)

# Simulate data
set.seed(42)
data = data.frame(
  x = seq(0, 35, by = 0.2),
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp::mcp(model, data, sample = FALSE)
data$y = empty$simulate(empty, data,
  cp_1 = 17,
  Intercept_1 = 10,
  sinx_1 = 10,
  Intercept_2 = 10,
  x_2 = 3,
  cosx_2 = 8,
  sigma_1 = 3
)

# Run sampling
fit = mcp(model, data, sample = sample)",



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
    cp_1 = 25,
    cp_2 = 75,
    Intercept_1 = 20,
    x_3 = 2,
    sigma_1 = 7,
    sigma_2 = 25,
    sigma_x_2 = -0.45
  )

# Run sampling
fit = mcp(model, data, adapt = 3000, sample = sample)",



varying = "# Define model
model = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope
)

# Simulate data
set.seed(42)
data = tibble::tibble(id = c('John', 'Benny', 'Rose', 'Cath', 'Bill', 'Erin')) %>%
  tidyr::expand_grid(x = seq(1, 100, by=4)) %>%
  dplyr::mutate(
    id_numeric = as.numeric(as.factor(id)),
    y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
  )
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
fit = mcp(model, data, cores = 3, sample = sample)"
)

# Run the code in an environment
assert_value(name, allowed = names(examples))
example_env = new.env()
with(example_env, eval(str2expression(examples[[name]])))

example_env$fit$call = examples[[name]]
class(example_env$fit$call) = c("mcptext", "character")

example_env$fit
}


#' @aliases mcp_example_data
#' @export
#' @describeIn mcp_example Conveniently get simulated data only.
mcp_example_data = function(name) {
  mcp_example(name, sample = FALSE)$data
}
