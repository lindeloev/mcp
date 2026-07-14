# Default prior specifications ---------------------------------------------

# Family defaults are stored as unresolved JAGS-scale rules. `.x` and `.y`
# are replaced by the actual change-point predictor and response names when a
# model is built. Keeping the rules symbolic makes family$default_prior useful
# to inspect and lets the same compiler handle defaults and user priors.

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
  response_location = "round(median(.y), 1)"
  response_scale = "max(2.5, round(mad(.y), 1))"

  if (link == "identity") {
    mu_location = response_location
    mu_scale = response_scale
  } else if (link == "log") {
    # The log link constrains the mean, not the observations. Calibrating from
    # log(y) would reject legitimate non-positive Gaussian responses.
    mu_location = "round(log(max(median(.y), 0.1)), 1)"
    mu_scale = "2.5"
  } else {
    return(NULL)
  }

  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "mu", "Intercept", paste0("dt(", mu_location, ", ", mu_scale, ", 3)"), "Robustly centered mean intercept with a brms-inspired minimum scale", "always",
    "mu", "dummy", paste0("dt(0, ", mu_scale, ", 3)"), "Regularizing mean contrast on the link scale", "always",
    "mu", "slope", paste0("dt(0, ", mu_scale, " / segment_width(.x), 3)"), "Regularizing mean slope scaled to the expected segment width", "always",
    "sigma", "Intercept", paste0("dt(0, ", response_scale, ", 3) T(0, )"), "Positive residual SD calibrated on the response scale", "always",
    "sigma", "dummy", paste0("dt(0, ", response_scale, ", 3)"), "Residual-SD contrast calibrated on the response scale", "always",
    "sigma", "slope", paste0("dt(0, ", response_scale, " / segment_width(.x), 3)"), "Residual-SD slope scaled to the expected segment width", "always"
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
    "mu", "Intercept", "dt(0, 2.5, 3)", "brms-inspired weak prior for the link-scale intercept", "always",
    "mu", "dummy", "dt(0, 2.5, 3)", "Regularizing categorical contrast on the link scale", "always",
    "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing link-scale slope scaled to the expected segment width", "always"
  )
}


default_priors_poisson = function(link) {
  response_location = "round(median(.y), 1)"
  response_scale = "max(2.5, round(mad(.y), 1))"

  if (link == "identity") {
    return(tibble::tribble(
      ~dpar, ~par_type, ~prior, ~description, ~condition,
      "mu", "Intercept", paste0("dt(", response_location, ", ", response_scale, ", 3) T(0, )"), "Positive count intercept calibrated on the response scale", "always",
      "mu", "dummy", paste0("dt(0, ", response_scale, ", 3)"), "Count contrast calibrated on the response scale", "always",
      "mu", "slope", paste0("dt(0, ", response_scale, " / segment_width(.x), 3)"), "Count slope scaled to the expected segment width", "always"
    ))
  }

  if (link != "log")
    return(NULL)

  count_y = "log(pmax(.y, 0.1))"
  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "mu", "Intercept", paste0("dt(round(median(", count_y, "), 1), max(2.5, round(mad(", count_y, "), 1)), 3)"), "Robustly centered log-count intercept with a brms-inspired minimum scale", "always",
    "mu", "dummy", "dt(0, 2.5, 3)", "Regularizing categorical contrast on the log scale", "always",
    "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing log-count slope scaled to the expected segment width", "always"
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
      "Robustly centered log-count intercept with a brms-inspired minimum scale", "always",
    "mu", "dummy", "dt(0, 2.5, 3)", "Regularizing categorical count contrast on the log scale", "always",
    "mu", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing log-count slope scaled to the expected segment width", "always",
    "shape", "Intercept", "dloginvgamma(0.4, 0.3)", "brms-compatible negative-binomial shape prior", "shape_constant",
    "shape", "Intercept", "dt(0, 2.5, 3)", "brms-inspired prior for a modeled log-shape intercept", "shape_modeled",
    "shape", "dummy", "dt(0, 2.5, 3)", "Regularizing shape contrast on the log scale", "always",
    "shape", "slope", "dt(0, 2.5 / segment_width(.x), 3)", "Regularizing log-shape slope scaled to the expected segment width", "always"
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
