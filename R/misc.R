#' RowSums of element-wise products.
#' Like an inner product, just vectorized for many y.
#' This ensures R <--> JAGS code compatibility
#'
#' @aliases inprod
#' @param x A matrix
#' @param y A matrix
#' @return A vector
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
inprod = function(x, y) {
  rowSums(x * y)
}


# Converts logical(0) to null. Returns x otherwise
logical0_to_null = function(x) {
  if (length(x) > 0)
    return(x)
  else
    return(NULL)
}


# if((a %in% b) == FALSE) --> if(a %notin% b)
`%notin%` = Negate(`%in%`)


# Is this a continuous vector?
is_continuous = function(x) {
  is.numeric(x) &
    length(unique(na.omit(x))) > 2
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


# Just returns cp_1, cp_2, etc.
get_cp_pars = function(pars) {
  pars$reg[stringr::str_detect(pars$reg, "^cp_[0-9]+$") & !stringr::str_detect(pars$reg, "^cp_[0-9]+_sd$")]  # cp_1 but not cp_1_sd
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
remove_terms = function(form, remove) {
  assert_types(form, "formula", len = c(2, 3))
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
  assert_types(form, "character", "formula", len = c(1, 3))
  if (is.character(form)) {
    # Add tilde
    if (!stringr::str_detect(form, "^(\\s|)~")) {
      form = paste0("~", form)
    }
    form = stats::as.formula(form)
  }

  form
}


#' Homogonize enumerating strings in mcp
#'
#' Nice for error messages.
#'
#' @aliases collapse_quote
#' @keywords internal
#' @param x A character vector
#' @return Character
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
and_collapse = function(x) {
  assert_types(x, "character")
  paste0(x, collapse = " and ")
}


#' Converts formula to string
#'
#' Note: this uses base R and circumvents the length-limitation of `deparse()`
#' and `format()`.
#'
#' @aliases formula_to_char
#' @keywords internal
#' @param form Any valid formula with any number of tildes.
#' @return A character.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
formula_to_char = function(form) {
  assert_types(form, "formula", len = c(1, 3))
  form_char = as.character(form)
  if (length(form_char) == 2 & form_char[1] == "~") {
    return(paste0(form_char, collapse = " "))
  } else if (length(form_char == 3) & form_char[1] == "~") {
    return(paste0(form_char[c(2, 3)], collapse = " ~ "))
  } else {
    stop_github("Could not decode formula ", deparse(form, width.cutoff = 500))
  }
}


#' Returns the right-hand-side of a formula
#'
#' @aliases get_rhs
#' @keywords internal
#' @param form Formula, e.g. `~x`, `y ~ x` or `y ~ z ~ x`
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs = function(form) {
  if (length(form) == 2) {
    return(form)
  } else if (length(form) == 3) {
    return(form[-2])
  }
}


#' Returns all vars in the RHS of an mcpmodel
#'
#' @aliases get_rhs_vars
#' @keywords internal
#' @inheritParams mcp
#' @return Character vector with unique term names
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs_vars = function(model) {
  assert_types(model, "mcpmodel")

  model %>%
    lapply(get_rhs) %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()
}

#' Returns all vars in the RHS of an mcpmodel
#'
#' @aliases get_model_vars
#' @keywords internal
#' @inheritParams mcp
#' @return Character vector with unique term names
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_model_vars = function(model) {
  assert_types(model, "mcpmodel")

  model %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()
}


#' Create model matrix from rhs_table
#'
#' cbinds rhs_table$matrix_data
#' @aliases get_rhs_matrix
#' @keywords internal
#' @param rhs_table The output of `get_rhs_table()`
#' @return matrix with one column for each row in `rhs_table`
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs_matrix = function(rhs_table) {
  suppressMessages(dplyr::bind_cols(rhs_table$matrix_data, .name_repair = "unique")) %>% # Suppress message about lacking column names
    as.matrix() %>%
    magrittr::set_colnames(rhs_table$code_name)
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
  assert_types(x, "character", len = 1)
  assert_ellipsis(...)
  cat(x)
}


# Set model environment to parent.frame() for prettier printing
# and because it was created in a different environment than inteded for use.
fix_model_environment = function(model) {
  assert_types(model, "mcpmodel")
  for (i in seq_along(model))
    environment(model[[i]]) = globalenv()
  model
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
data$price = NULL  # Simulate new
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
if (sample == TRUE)
  fit = mcp(model, data)",



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
  N = sample(10, 100, replace=TRUE),
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
if (sample == TRUE)
  fit = mcp(model, data, family = binomial(), adapt = 5000)
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
if (sample == TRUE)
  fit = mcp(model, data)",



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
if (sample == TRUE)
  fit = mcp(model, data, par_x = 'x')",



multiple = "# Define model
model = list(
  y ~ 1 + x:group + z,
  ~ 1 + x + group,
  ~ 0 + I(x^2)
)

# Simulate data
data = data.frame(
  x = 1:220,
  group = rep(c('A', 'B', 'C', 'D'), 55),
  z = rnorm(220, mean = 1:220, sd = 25),
  y = 2.  # or whatever signals 'numeric'. Will be replaced by simulation below.
)
empty = mcp(model, data, sample = FALSE, par_x = 'x')
data$y = empty$simulate(empty, data,
  cp_1 = 100,
  cp_2 = 180,

  Intercept_1 = 10,
  z_1 = 0.2,
  xgroupA_1 = -0.75,
  xgroupB_1 = -0.25,
  xgroupC_1 = 0.25,
  xgroupD_1 = 0.75,

  Intercept_2 = 40,
  x_2 = -0.8,
  groupB_2 = 15,
  groupC_2 = 30,
  groupD_2 = 45,

  xE2_3 = 0.1,

  sigma_1 = 5
)

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data, par_x = 'x', cores = 3)",



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
if (sample == TRUE)
  fit = mcp(model, data)",



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
if (sample == TRUE)
  fit = mcp(model, data)",



variance = "# Define model
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
    cp_1 = 25,
    cp_2 = 75,
    Intercept_1 = 20,
    x_3 = 2,
    sigma_1 = 7,
    sigma_2 = 25,
    sigma_x_2 = -0.45
  )

# Run sampling
if (sample == TRUE)
  fit = mcp(model, data)",



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
if (sample == TRUE)
  fit = mcp(model, data)"
  )

  # Run the code
  assert_value(name, allowed = names(examples))
  eval(parse(text = examples[[name]]))

  # Get stuff ready for return
  model = fix_model_environment(model)

  class(call) = c("mcptext", "character")
  class(model) = c("mcplist", "list")

  last_col = dplyr::pull(data, -1)  # The response column is always the last column
  simulated = attr(last_col, "simulated")

  if (sample == FALSE)
    fit = NULL

  # Return
  list(
    model = model,  # Bind them to global workspace for nicer display
    data = data,
    simulated = simulated,  # response is always the last column
    fit = fit,
    call = call
  )
}
