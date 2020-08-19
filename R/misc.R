# Converts logical(0) to null. Returns x otherwise
logical0_to_null = function(x) {
  if (length(x) > 0)
    return(x)
  else return(NULL)
}


#' Extracts the order from ARMA parameter name(s)
#'
#' If several names are provided (vector), it returns the maximum. If `pars_arma`
#' is an empty string, it returns `0`.
#'
#' @aliases get_arma_order
#' @keywords internal
#' @param pars_arma Character vector
#' @return integer
get_arma_order = function(pars_arma) {
  if (length(pars_arma) > 0) {
    order_str = sub("(ma|ar)([0-9]+).*", "\\2", pars_arma)
    order_max = max(as.numeric(order_str))
    return(order_max)
  } else {
    return(0)
  }
}


# Ask reminder questions for CRAN export
release_questions = function() {
  c(
    "Have you run the test of fits? options(test_mcp_fits = TRUE)",
    "Have you built the README plots and checked them? source('vignettes/figures/make_README_plots.R')",
    "Have you re-built the site using pkgdown::build_site() AFTER deleting caches of articles?",
    "Have you checked all articles and plots after re-building the site?",
    "Have you run the script to insert the correct logo.png in the HTML meta?"
  )
}


#' Remove varying or population terms from a formula
#'
#' WARNING: removes response side from the formula
#'
#' @aliases remove_terms
#' @keywords internal
#' @param form A formula
#' @param remove Either "varying" or "population". These are removed.
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
remove_terms = function(form, remove) {
  assert_types(form, "formula")
  assert_value(remove, allowed = c("varying", "population"))

  # Find terms with "|"
  attrs = attributes(stats::terms(form))
  term.labels = attrs$term.labels
  varying_bool = stringr::str_detect(term.labels, "\\|")

  # Add parenthesis back to them
  term.labels[varying_bool] = paste0("(", term.labels[varying_bool], ")")

  # Remove non-matching types
  if (remove == "varying") {
    term.labels = term.labels[!varying_bool]
    term.labels = c(attrs$intercept, term.labels)  # Add intercept indicator
  } else if (remove == "population") {
    term.labels = term.labels[varying_bool]
  }

  # Build formula from terms and return
  if (length(term.labels) == 0) {
    return(NULL)
  } else {
    formula_terms = paste0(term.labels, collapse = " + ")
    formula_str = paste0("~", formula_terms)
    return(stats::as.formula(formula_str, env=globalenv()))
  }
}


#' Takes any formula-like input and returns a formula
#' @aliases to_formula
#' @keywords internal
#' @param form Formula or character (with or without initial tilde/"~")
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
to_formula = function(form) {
  assert_types(form, "character", "formula")
  if (is.character(form)) {
    # Add tilde
    if (!stringr::str_detect(form, "^(\\s|)~")) {
      form = paste0("~", form)
    }
    form = stats::as.formula(form)
  }

  return(form)
}

#' Expand samples with quantiles
#'
#' TO DO: implement using `fitted()` and `predict()` but avoid double-computing the samples? E.g.:
#' `get_quantiles2 = function(fit, quantiles, facet_by = NULL) {`
#'   `fitted(fit, probs = c(0.1, 0.5, 0.9), newdata = data.frame(x = c(11, 50, 100))) %>%`
#'   `tidyr::pivot_longer(tidyselect::starts_with("Q")) %>%`
#'   `dplyr::mutate(quantile = stringr::str_remove(name, "Q") %>% as.numeric() / 100)`
#' `}`
#'
#' @aliases get_quantiles
#' @keywords internal
#' @inheritParams plot.mcpfit
#' @param samples A tidybayes tibble
#' @param quantiles Vector of quantiles (0.0 to 1.0)
#' @param xvar An rlang::sym() with the name of the x-col in `samples`
#' @param yvar An rlang::sym() with the name of the response col in `samples`
#' @return A tidybayes long format tibble with the column "quantile"
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#'
get_quantiles = function(samples, quantiles, xvar, yvar, facet_by = NULL) {
  # Trick to declare no facet = common group for all
  if (is.null(facet_by))
    facet_by = xvar

  # Return data with added quantiles
  samples %>%
    tidyr::expand_grid(quantile = quantiles) %>%

    # Now compute the quantile for each parameter, quantile, and (optionally) facet:
    dplyr::group_by(!!xvar, .data$quantile) %>%
    dplyr::group_by(!!rlang::sym(facet_by), .add = TRUE) %>%
    dplyr::summarise(
      y = stats::quantile(!!yvar, probs = .data$quantile[1])
    )
}

# Hack to make R CMD pass for function geom_cp_density()
utils::globalVariables(c("value", "..scaled..", ".chain", "cp_name", "."))


#' Print mcplist
#'
#' Shows a list in a more condensed format using `str(list)`.
#' @aliases print.mcplist
#' @inheritParams print.mcpfit
print.mcplist = function(x, ...) {
  assert_types(x, "list")
  assert_ellipsis(...)

  # For all-formula list (typically fit$model)
  if (all(sapply(x, is.formula)) == TRUE) {
    cat(paste0("List of ", length(x), "\n"))
    for (i in x) {
      cat(" $ ")
      print(i)
    }
  } else {
    # Other lists
    utils::str(x, vec.len = Inf, give.head = FALSE, give.attr = FALSE)
  }
}


#' Nice printing texts
#'
#' Useful for `print(fit$jags_code)`, `print(mcp_demo$call)`, etc.
#'
#' @aliases print.mcptext
#' @param x Character, often with newlines.
#' @param ... Currently ignored.
#' @return NULL
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' mytext = "line1 = 2\n line2 = 'horse'"
#' class(mytext) = "mcptext"
#' print(mytext)
print.mcptext = function(x, ...) {
  assert_types(x, "character")
  assert_ellipsis(...)
  cat(x)
}


#' Get example models and data
#'
#' @aliases mcp_example
#' @param name Name of the example. One of:
#'  * `"demo"`: Two change points between intercepts and joined/disjoined slopes.
#'  * `"ar"`: One change point in autoregressive residuals.
#'  * `"binomial"`: Binomial with two change points. Much like `"demo"` on a logit scale.
#'  * `"intercepts"`: An intercept-only change point.
#'  * `rel_prior`: Relative parameterization and informative priors.
#'  * `"quadratic"`: A change point to a quadratic segment.
#'  * `"trigonometric"`: Trigonometric/seasonal data and model.
#'  * `"varying"`: Varying / hierarchical change points.
#'  * `"variance"`: A change in variance, including a variance slope.
#' @param sample TRUE (run `fit = mcp(model, data, ...)`) or FALSE.
#' @return List with
#'  * `model`: A list of formulas
#'  * `data`: The simulated data
#'  * `simulated`: The parameters used for simulating the data.
#'  * `fit`: an `mcpfit` if `sample = TRUE`,
#'  * `call`: the code to run the above.
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' \donttest{
#' ex = mcp_example("demo")
#' plot(ex$data)  # Plot data
#' print(ex$simulated)  # See true parameters used to simulate
#' print(ex$call)  # See how the data was simulated
#'
#' # Fit the model. Either...
#' fit = mcp(ex$model, ex$data)
#' plot(fit)
#'
#' ex_with_fit = mcp_example("demo", sample = TRUE)
#' plot(ex_with_fit$fit)
#'}
mcp_example = function(name, sample = FALSE) {
  assert_value(name, c("demo", "ar", "binomial", "intercepts", "rel_prior", "quadratic", "trigonometric", "varying", "variance"))
  data = data.frame()  # To make R CMD Check happy.
  if (name == "demo") {
    call = "# Define model
model = list(
  response ~ 1,
  ~ 0 + time,
  ~ 1 + time
)

# Simulate data
empty = mcp(model, sample = FALSE)
set.seed(40)
data = tibble::tibble(
  time = runif(100, 0, 100),
  response = empty$simulate(
    time,
    cp_1 = 30,
    cp_2 = 70,
    int_1 = 10,
    int_3 = 0,
    sigma = 4,
    time_2 = 0.5,
    time_3 = -0.2)
)

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data)"
  } else if (name == "ar") {
    call = "# Define model
model = list(
  price ~ 1 + ar(2),
  ~ 0 + time + ar(1)
)

# Simulate data
empty = mcp(model, sample = FALSE)
set.seed(42)
data = tibble::tibble(
  time = 1:200,
  price = empty$simulate(
    time,
    cp_1 = 120,
    int_1 = 20,
    time_2 = 0.5,
    sigma_1 = 5,
    ar1_1 = 0.7,
    ar2_1 = 0.2,
    ar1_2 = -0.4)
)

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data)"
  } else if (name == "binomial") {
    call = "# Define model
model = list(
  y | trials(N) ~ 1,  # constant rate
  ~ 0 + x,  # joined changing rate
  ~ 1 + x  # disjoined changing rate
)

# Simulate data
empty = mcp(model, family = binomial(), sample = FALSE)
set.seed(42)
data = tibble::tibble(
  x = 1:100,
  N = sample(10, length(x), replace=TRUE),
  y = empty$simulate(
    x = x,
    N,
    cp_1 = 30,
    cp_2 = 70,
    int_1 = 2,
    int_3 = 0.4,
    x_2 = -0.2,
    x_3 = 0.05)
)

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data, family = binomial(), adapt = 5000)
"
  } else if (name == "intercepts") {
    call = "# Define model
model = list(
  y ~ 1,
  ~ 1
)

# Simulate data
empty = mcp(model, sample = FALSE, par_x = \"x\")
set.seed(40)
data = tibble::tibble(
  x = runif(100, 0, 100),
  y = empty$simulate(
    x,
    cp_1 = 50,
    int_1 = 10,
    int_2 = 20,
    sigma = 8)
)

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data, par_x = \"x\")"
  } else if (name == "quadratic") {
    call = "# Define model
model = list(
  y ~ 1,
  ~ 0 + x + I(x^2)
)

# Simulate data
empty = mcp::mcp(model, sample = FALSE)
set.seed(42)
data = tibble::tibble(
  x = seq(0, 40, by = 0.5),
  y = empty$simulate(
    x,
    cp_1 = 15,
    int_1 = 10,
    sigma = 30,
    x_2 = -30,
    x_2_E2 = 1.5)
)

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data)"
  } else if (name == "rel_prior") {
    call = "# Define model
model = list(
  y ~ 1 + x,
  ~ rel(1) + rel(x),
  rel(1) ~ 0 + x
)

# Simulate data
empty = mcp::mcp(model, sample = FALSE)
set.seed(40)
data = tibble::tibble(
  x = 1:100,
  y = empty$simulate(
    x,
    cp_1 = 25,
    cp_2 = 40,
    int_1 = 25,
    int_2 = -10,
    sigma = 7,
    x_1 = 1,
    x_2 = -3,
    x_3 = 0.5)
)

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data)"
  } else if (name == "trigonometric") {
    call = "model = list(
  y ~ 1 + sin(x),
  ~ 1 + cos(x) + x
)

# Simulate data
empty = mcp::mcp(model, sample = FALSE)
set.seed(40)
data = tibble::tibble(
  x = seq(0, 35, by = 0.2),
  y = empty$simulate(
    x,
    cp_1 = 17,
    int_1 = 10,
    x_1_sin = 10,
    x_2_cos = 8,
    int_2 = 10,
    x_2 = 3,
    sigma = 3
  )
)

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data)"
  } else if (name == "variance") {
    call = "# Define model
model = list(
  y ~ 1,
  ~ 0 + sigma(1 + x),
  ~ 0 + x
)

# Simulate data
empty = mcp::mcp(model, sample = FALSE)
set.seed(30)
data = tibble::tibble(
  x = 1:100,
  y = empty$simulate(
    x,
    cp_1 = 25,
    cp_2 = 75,
    int_1 = 20,
    x_3 = 2,
    sigma_1 = 7,
    sigma_2 = 25,
    sigma_x_2 = -0.45)
  )

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data)"
  } else if (name == "varying") {
    call = "# Define model
model = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope
)

# Simulate data
empty = mcp::mcp(model, sample=FALSE)
data = tibble::tibble(id = c(\"John\", \"Benny\", \"Rose\", \"Cath\", \"Bill\", \"Erin\")) %>%
  tidyr::expand_grid(x = seq(1, 100, by=4)) %>%
  dplyr::mutate(
    id_numeric = as.numeric(as.factor(id)),
    y = empty$simulate(
      x,
      cp_1 = 40,
      cp_1_id = 7*(id_numeric - mean(id_numeric)),
      int_1 = 15,
      x_1 = 3,
      x_2 = -2,
      sigma = 25
    )
  )

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data)"
  } else {
    stop("Unknown `name`: ", name)
  }

  # Run the code
  eval(parse(text = call))

  # Get stuff ready for return
  for (i in seq_along(model))
    environment(model[[i]]) = parent.frame()

  class(call) = c("mcptext", "character")
  class(model) = c("mcplist", "list")

  last_col = dplyr::pull(data, -1)  # The response column is always the last column
  simulated = attr(last_col, "simulated")

  if (sample == FALSE)
    fit = NULL

  return(list(
    model = model,  # Bind them to global workspace for nicer display
    data = data,
    simulated = simulated,  # response is always the last column
    fit = fit,
    call = call
  ))
}
