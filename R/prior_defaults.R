# Default prior specifications ---------------------------------------------

# Family defaults are stored as unresolved JAGS-scale rules. `.x` and `.y`
# are replaced by the actual change-point predictor and response names when a
# model is built. Keeping the rules symbolic makes family$default_prior useful
# to inspect and lets the same compiler handle defaults and user priors.

linked_response_rule = function(link) {
  switch(
    link,
    identity = ".y",
    log = "log(.y)",
    NULL
  )
}


normalize_family_default_priors = function(defaults) {
  assert_types(defaults, "tibble", "data.frame")
  assert_data_cols(defaults, c("dpar", "par_type", "prior"))
  defaults = tibble::as_tibble(defaults)
  if (!"description" %in% names(defaults))
    defaults$description = NA_character_
  if (!"condition" %in% names(defaults))
    defaults$condition = "always"
  defaults[, c("dpar", "par_type", "prior", "description", "condition")]
}


default_priors_gaussian = function(link) {
  link_y = linked_response_rule(link)
  if (is.null(link_y))
    return(NULL)

  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "mu", "Intercept", paste0("dt(median(", link_y, "), mad(", link_y, "), 3)"), "Current data-calibrated mean intercept prior", "always",
    "mu", "dummy", paste0("dt(0, mad(", link_y, "), 3)"), "Current data-calibrated categorical contrast prior", "always",
    "mu", "slope", paste0("dt(0, mad(", link_y, ") / segment_width(.x), 3)"), "Current slope prior scaled to the expected change-point segment width", "always",
    "sigma", "Intercept", paste0("dt(0, mad(", link_y, "), 3) T(0, )"), "Current residual-SD intercept prior", "always",
    "sigma", "dummy", paste0("dt(0, mad(", link_y, "), 3)"), "Current residual-SD contrast prior", "always",
    "sigma", "slope", paste0("dt(0, mad(", link_y, ") / segment_width(.x), 3)"), "Current residual-SD slope prior scaled to the expected change-point segment width", "always"
  )
}


default_priors_binomial = function(link) {
  if (link == "identity") {
    return(tibble::tribble(
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", "dbeta(1, 1)", "Uniform probability intercept", "always",
      "mu", "dummy", "dunif(-1, 1)", "Probability difference between levels", "always",
      "mu", "slope", "dt(0, 1 / segment_width(.x), 3)", "Current slope prior scaled to the expected change-point segment width", "always"
    ))
  }

  if (link %notin% c("logit", "probit"))
    return(NULL)

  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "mu", "Intercept", "dt(0, 2.5, 3)", "Current weak prior for the probability-scale intercept", "always",
    "mu", "dummy", "dt(0, 2.5, 3)", "Current weak prior for a categorical probability contrast", "always",
    "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Current slope prior scaled to the expected change-point segment width", "always"
  )
}


default_priors_poisson = function(link) {
  if (link == "identity") {
    return(tibble::tribble(
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", "dt(median(.y), median(.y), 3)", "Current identity-count intercept prior", "always",
      "mu", "dummy", "dt(0, median(.y), 3)", "Current data-calibrated categorical contrast prior", "always",
      "mu", "slope", "dt(0, median(.y) / segment_width(.x), 3)", "Current slope prior scaled to the expected change-point segment width", "always"
    ))
  }

  if (link != "log")
    return(NULL)

  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "mu", "Intercept", "dt(0, 10, 3)", "Current weak prior for the log-count intercept", "always",
    "mu", "dummy", "dt(0, 10, 3)", "Current weak prior for a categorical log-count contrast", "always",
    "mu", "slope", "dt(0, 10, 3)", "Current weak prior for a log-count slope", "always"
  )
}


default_priors_negbinomial = function(link) {
  if (link != "log")
    return(NULL)

  count_y = "log(pmax(.y, 0.1))"
  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "mu", "Intercept",
      paste0("dt(round(median(", count_y, "), 1), max(2.5, round(mad(", count_y, "), 1)), 3)"),
      "Current data-calibrated log-count intercept prior", "always",
    "mu", "dummy", "dt(0, 2.5, 3)", "Current link-scale categorical count contrast prior", "always",
    "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Current slope prior scaled to the expected change-point segment width", "always",
    "shape", "Intercept", "dloginvgamma(0.4, 0.3)", "brms-compatible negative-binomial shape prior", "shape_constant",
    "shape", "Intercept", "dt(0, 2.5, 3)", "Current link-scale shape intercept prior", "shape_modeled",
    "shape", "dummy", "dt(0, 2.5, 3)", "Current link-scale shape contrast prior", "always",
    "shape", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Current shape-slope prior scaled to the expected change-point segment width", "always"
  )
}


get_family_default_priors = function(family) {
  builder = switch(
    family$family,
    gaussian = default_priors_gaussian,
    binomial = default_priors_binomial,
    bernoulli = default_priors_binomial,
    poisson = default_priors_poisson,
    negbinomial = default_priors_negbinomial,
    NULL
  )
  defaults = if (is.null(builder)) NULL else builder(family$link)
  if (is.null(defaults)) {
    stop(
      "mcp has no default priors for ", family$family, "(link = \"", family$link,
      "\") so it's likely not supported. See `mcpfamily()` on how to create a custom family."
    )
  }
  defaults
}


default_arma_priors = function() {
  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "ar", "Intercept", "dunif(-1, 1)", "Bounded dependence coefficient", "always",
    "ar", "dummy", "dt(0, 1, 3)", "Current weak prior for a categorical dependence coefficient", "always",
    "ar", "slope", "dt(0, 1 / (max(.x) - min(.x)), 3)", "One coefficient-unit change across the observed change-point span", "always",
    "ma", "Intercept", "dunif(-1, 1)", "Bounded dependence coefficient", "always",
    "ma", "dummy", "dt(0, 1, 3)", "Current weak prior for a categorical dependence coefficient", "always",
    "ma", "slope", "dt(0, 1 / (max(.x) - min(.x)), 3)", "One coefficient-unit change across the observed change-point span", "always"
  )
}
