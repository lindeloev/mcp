# ABOUT: Assemble the full JAGS model around the linear predictors from
# get_formula_jags(): priors, GARMA recursion, link functions, and likelihood
# -----------------

# JAGS names for the per-observation nodes used by the family's likelihood code,
# e.g., the response `response[i_]`, distributional parameters `sigma_[i_]`,
# and auxiliary response columns such as `trials(N)` as `N[i_]`.
get_jags_family_context = function(data_columns) {
  list(
    y = paste0(data_columns$response, "[i_]"),
    threshold = "garma_threshold_[i_]",
    dpar = function(name) paste0(name, "_[i_]"),
    local = function(name) paste0(name, "_[i_]"),
    aux = function(name, default = NULL) {
      if (is.null(data_columns[[name]])) {
        if (!is.null(default))
          return(default)
        stop_github("Missing response auxiliary ", name, "() in JAGS family context.")
      }
      paste0(data_columns[[name]], "[i_]")
    }
  )
}


#' Make JAGS code for Multiple Change Point model
#'
#' @param model_tables Model tables, as in `get_fit_model_tables(fit)`.
#' @param prior_table Resolved priors and metadata from `get_prior_table()`.
#' @param family An mcpfamily object.
#' @return String. A JAGS model.
#' @noRd
get_jags_code = function(model_tables, prior_table, family) {
  group_predictors = definition_rows(model_tables$group_predictors)
  par_x = model_tables$data_columns$par_x
  ar_order = get_arma_order(model_tables$predictors, "ar")
  ma_order = get_arma_order(model_tables$predictors, "ma")

  # Begin building JAGS model.
  # Change points come first and in order, from the fixed cp_0 = min(x) to cp_K = max(x).
  K = nrow(model_tables$cps) + 1
  jags_code = paste0("model {
  # Change points
  cp_0 = min(", par_x, ")\n")

  #######################
  # CHANGE POINT PRIORS #
  #######################
  # Each change point keeps its specified prior, and cp_order_ restricts them to ordered
  # change points with any minimum segment lengths (see add_nodes()). Group-level change
  # points come right after their population-level change point.
  cps = model_tables$cps
  prior_row = function(name) prior_table[match(name, prior_table$parameter), ]
  for (j in seq_len(nrow(cps))) {
    jags_code = paste0(jags_code, get_prior_str(prior_row(cps$name[j])))
    if (!is.na(cps$group_col[j]))
      jags_code = paste0(jags_code, get_prior_str(prior_row(cps$sd_name[j])), get_prior_str(prior_row(cps$group_name[j]), cps$group_col[j]))
  }
  jags_code = paste0(jags_code, "  cp_", K, " = max(", par_x, ")\n")

  # Observed nodes (data: cp_order_ = 1 and cp_weight_ = 0) that restrict the priors above
  cp_rows = prior_row(cps$name)
  if (nrow(cps) > 1) {
    min_lengths = cp_rows$min_length
    steps = paste0("step(cp_", 2:nrow(cps), " - cp_", 2:nrow(cps) - 1, ifelse(min_lengths[-1] == "0", "", paste0(" - ", min_lengths[-1])), ")")
    jags_code = paste0(jags_code, "  cp_order_ ~ dbern(", paste(steps, collapse = " * "), ")  # Ordered change points", if (any(min_lengths != "0")) " with minimum segment lengths", "\n")
  }
  alpha = if (nrow(cps) > 0 && cp_rows$parsed[[1]]$dist %in% "dirichlet") as.numeric(cp_rows$parsed[[1]]$args) else 1  # %in% is FALSE for dist = NA
  if (alpha != 1) {
    # Zeros trick: dpois(0 | lambda) adds -lambda = (alpha - 1) * sum(log(segment length / range of x)) to the log density
    log_lengths = paste0("log((cp_", seq_len(K), " - cp_", seq_len(K) - 1, ") / (cp_", K, " - cp_0))", collapse = " + ")
    jags_code = paste0(jags_code, "  cp_weight_ ~ dpois(-", format_prior_number(alpha - 1), " * (", log_lengths, "))  # dirichlet(", format_prior_number(alpha), ") weight on the segment lengths\n")
  }
  jags_code = paste0(jags_code, "
  # Priors for population-level effects\n")

  ################
  # OTHER PRIORS #
  ################
  # Split scalar/population priors from group-indexed coefficient priors.
  is_cp = prior_table$parameter %in% c(cps$name, cps$sd_name, cps$group_name)
  group_cols = stats::setNames(group_predictors$group_col, group_predictors$name)
  is_group = prior_table$parameter %in% names(group_cols)
  for (i in which(!is_cp & !is_group))
    jags_code = paste0(jags_code, get_prior_str(prior_table[i, ]))
  if (any(!is_cp & is_group)) {
    jags_code = paste0(jags_code, "\n  # Priors for group-level effects\n")
    for (i in which(!is_cp & is_group))
      jags_code = paste0(jags_code, get_prior_str(prior_table[i, ], group_cols[[prior_table$parameter[i]]]))
  }


  #########
  # GARMA #
  #########
  has_arma_terms = !is.na(ar_order) || !is.na(ma_order)
  if (has_arma_terms)
    jags_code = paste0(jags_code, get_garma_jags_code(ar_order, ma_order, par_x, series = !is.null(model_tables$data_columns$series)))


  ###########
  # FORMULA #
  ###########
  jags_code = paste0(jags_code, "
  # Model and likelihood
  for (i_ in 1:length(", par_x, ")) {")

  # Add JAGS code for fitted values and indent it
  jags_code = paste0(jags_code, gsub("\n", "\n    ", get_formula_jags(model_tables, family, stats::setNames(prior_table$node, prior_table$parameter))))


  ##############
  # LIKELIHOOD #
  ##############
  context = get_jags_family_context(model_tables$data_columns)
  jags_code = paste0(jags_code, "\n\n    # Likelihood for family = ", family$family, "()
    ")

  # Transform link-scale predictors to distribution-scale parameters. These are
  # deterministic JAGS nodes and are not monitored by default.
  for (dpar in family$dpar_specs$dpar) {
    dpar_spec = get_dpar_spec(family, dpar)
    link_code = paste0("link_", dpar, "_[i_]")
    if (dpar == "mu" && has_arma_terms)
      link_code = paste0(link_code, " + resid_garma_[i_]")

    linkinv_str = get_link_str(dpar_spec$link, inverse = TRUE)
    response_code = ifelse(
      linkinv_str == "",
      link_code,
      paste0(linkinv_str, "(", link_code, ")")
    )
    if (!is.na(dpar_spec$floor))  # Numerical floor, as in evaluate_model()
      response_code = paste0("max(", format(dpar_spec$floor, scientific = TRUE), ", ", response_code, ")")

    jags_code = paste0(jags_code, dpar, "_[i_] = ", response_code, "\n    ")
  }

  likelihood = family$backends$jags$likelihood(context)
  jags_code = paste0(jags_code, paste(likelihood, collapse = "\n    "))

  # Compute link-scale residuals for GARMA
  if (has_arma_terms) {
    garma_y = family$garma$observed_jags(context)
    linkfun_str = get_link_str(family$link)
    garma_link_y = if (linkfun_str == "") "garma_y_[i_]" else
      paste0(linkfun_str, "(garma_y_[i_])")
    jags_code = paste0(
      jags_code,
      "\n    garma_y_[i_] = ", garma_y,
      "\n    garma_link_y_[i_] = ", garma_link_y
    )
    jags_code = paste0(jags_code, "\n    resid_abs_[i_] = garma_link_y_[i_] - link_mu_[i_]")
    if (!is.na(ma_order))
      jags_code = paste0(jags_code, "\n    resid_ma_[i_] = garma_link_y_[i_] - link_mu_[i_] - resid_garma_[i_]")
  }

  # Close the likelihood loop and the model
  paste0(jags_code, "
  }
}")
}


# How each parameter is sampled in JAGS. Adds these columns to the prior table, which
# get_prior_str() and the initial values read:
# - node: the stochastic JAGS node. Usually the parameter itself; NA for constants and expressions.
# - node_prior: the parsed prior on `node` (see parse_prior()).
# - derive: JAGS expression for the parameter from `node`; NA where `node` is the parameter.
# - sampled_as: what `node` is, for the JAGS comment; NA where `node` is the parameter.
# - min_length: minimum length of the segment that ends at each change point, in cp_order_.
#   NA for other parameters.
# Here, `node` differs from the parameter for dinvgamma() priors and group-level change
# points. add_reparams() adds more. `group_names` are the group-level parameters.
add_nodes = function(prior_table, cps, group_names) {
  prior_table$node = ifelse(prior_table$kind == "distribution", prior_table$parameter, NA_character_)
  prior_table$node_prior = prior_table$parsed
  prior_table$derive = NA_character_
  prior_table$sampled_as = NA_character_
  prior_table$min_length = NA_character_

  # JAGS does not support dinvgamma. If inverse ~ Gamma(a, b), then 1 / inverse ~ InvGamma(a, b).
  for (i in which(vapply(prior_table$parsed, function(parsed) parsed$dist %in% "dinvgamma", logical(1)))) {  # %in% is FALSE for dist = NA
    parsed = prior_table$parsed[[i]]
    if (prior_table$parameter[i] %in% group_names)
      stop("dinvgamma() is currently only supported for population-level parameters.")
    if (length(parsed$args) != 2 || any(nzchar(parsed$trunc)))
      stop("Expected dinvgamma(shape, scale) without truncation. Got '", format_prior(parsed), "'.")
    parsed$dist = "dgamma"
    prior_table$node[i] = paste0(prior_table$parameter[i], "_inverse")
    prior_table$node_prior[[i]] = parsed
    prior_table$derive[i] = paste0("1 / ", prior_table$node[i])
    prior_table$sampled_as[i] = "its inverse"
  }

  # Change points are ordered in cp_order_, which also holds minimum segment lengths, e.g.,
  # T(cp_1 + 20, ) is T(cp_0, ) with a minimum length of 20. dirichlet(alpha) is the default
  # uniform prior with a weight on the segment lengths (see get_jags_code()).
  cp_rows = match(cps$name, prior_table$parameter)
  if (any(vapply(prior_table$parsed[cp_rows], function(parsed) parsed$dist %in% "dirichlet", logical(1)))) {  # %in% is FALSE for dist = NA
    # overlay_user_prior_specs() ensures that all or none are dirichlet()
    alpha = vapply(prior_table$parsed[cp_rows], function(parsed) if (length(parsed$args) == 1) suppressWarnings(as.numeric(parsed$args)) else NA_real_, numeric(1))
    if (any(!is.finite(alpha)) || any(alpha < 1))
      stop("`dirichlet(alpha)` requires one finite alpha >= 1.")
    if (any(alpha != alpha[1]))
      stop("All `dirichlet(alpha)` change point priors must use the same alpha.")
  }
  locations = character()  # Earlier group-level change points sampled as locations, e.g., c("cp_1 + cp_1_id[" = "cp_1_id_location[")
  for (j in seq_len(nrow(cps))) {
    # Population-level change points: dirichlet() is uniform here, and lower bounds at the
    # previous change point become cp_0 with a minimum length
    i = cp_rows[j]
    parsed = prior_table$parsed[[i]]
    prior_table$min_length[i] = "0"
    if (parsed$dist %in% "dirichlet") {  # %in% is FALSE for dist = NA
      prior_table$node_prior[[i]] = parse_prior(paste0("dunif(cp_0, cp_", nrow(cps) + 1, ")"))
    } else if (!is.na(parsed$dist)) {
      bounds = get_prior_bounds(parsed)
      min_length = get_min_segment_length(bounds, j, cps$name)
      if (!is.null(min_length)) {
        prior_table$min_length[i] = min_length
        prior_table$node_prior[[i]] = set_prior_bounds(parsed, c("cp_0", bounds[2]))
      }
    }

    # Group-level change points. Bounds at earlier group-specific change points, e.g.,
    # "(cp_1 + cp_1_id[id_])", use their sampled locations.
    group_row = match(cps$group_name[j], prior_table$parameter)
    parsed = prior_table$parsed[[group_row]]
    if (is.na(cps$group_col[j]) || is.na(parsed$dist))
      next
    for (deviation in names(locations))
      parsed$trunc = gsub(deviation, locations[[deviation]], parsed$trunc, fixed = TRUE)

    # dnorm() deviations are sampled as each group's absolute location, cp_1_id_location ~ dnorm(cp_1 + mean, sd)
    if (parsed$dist == "dnorm" && length(parsed$args) == 2) {
      population = cps$name[j]
      parsed$args[1] = if (parsed$args[1] == "0") population else paste0(population, " + ", parsed$args[1])
      # Deviation bounds such as "cp_0 - cp_1" are simply "cp_0" on the location
      population_suffix = paste0(" - ", population)
      parsed$trunc = vapply(parsed$trunc, function(bound) {
        if (!nzchar(bound))
          return("")
        if (endsWith(bound, population_suffix))
          return(substr(bound, 1, nchar(bound) - nchar(population_suffix)))
        paste0(population, " + (", bound, ")")
      }, character(1), USE.NAMES = FALSE)
      prior_table$node[group_row] = paste0(cps$group_name[j], "_location")
      prior_table$derive[group_row] = paste0(prior_table$node[group_row], "[", cps$group_col[j], "_] - ", population)
      prior_table$sampled_as[group_row] = "each group's absolute change point"
      locations[[paste0(population, " + ", cps$group_name[j], "[")]] = paste0(prior_table$node[group_row], "[")
    }
    prior_table$node_prior[[group_row]] = parsed
  }
  prior_table
}


# Get the JAGS code (one or more lines) for the prior on one parameter, given its row
# in the prior table. Writes the node with its prior and derives the parameter from
# it (see add_nodes()). Constants and expressions are written as they are.
# - group_col: NULL for population-level priors. Otherwise the grouping column, and
#   the prior is written in a loop over its levels.
get_prior_str = function(row, group_col = NULL) {
  name = row$parameter
  description = gsub("[\r\n]+", " ", row$description)
  index = if (!is.null(group_col)) paste0("[", group_col, "_]") else ""
  if (is.na(row$node)) {
    # Fixed constants/expressions for group-level effects must be vectors over group levels
    value = row$value
    if (!is.null(group_col))
      value = gsub(paste0("([A-Za-z0-9_]+_", group_col, ")\\b(?!\\[)"), paste0("\\1[", group_col, "_]"), value, perl = TRUE)
    lines = paste0(name, index, " = ", value, "  # ", description)
  } else {
    # Convert SD/scale to JAGS precision and write a stochastic node
    sampled_as = if (is.na(row$sampled_as)) "" else paste0("; sampled as ", row$sampled_as)
    lines = paste0(row$node, index, " ~ ", format_prior(prior_to_jags(row$node_prior[[1]])), "  # ", description, sampled_as)
    if (!is.na(row$derive))
      lines = c(lines, paste0(name, index, " = ", row$derive))
  }

  if (is.null(group_col))
    return(paste0("  ", lines, "\n", collapse = ""))
  paste0("  for (", group_col, "_ in 1:n_unique_", group_col, ") {\n", paste0("    ", lines, "\n", collapse = ""), "  }\n")
}


# Translate a scale-parameterized prior to JAGS ---------------------------

# mcp's legacy string syntax uses conventional distribution scales. JAGS
# uses inverse variance for normal-like distributions and inverse scale for
# double-exponential and logistic distributions. Takes and returns a parsed prior.
prior_to_jags = function(parsed) {
  if (identical(parsed$dist, "dcauchy"))
    stop("JAGS does not support `dcauchy()`; use `dt(location, scale, 1)` instead.")

  exponent = jags_scale_exponent(parsed$dist)
  if (is.null(exponent))
    return(parsed)
  if (length(parsed$args) < 2)
    stop("Expected a scale as the second argument of '", format_prior(parsed), "'.")

  denominator = paste0("(", gsub(" ", "", parsed$args[2]), ")")
  if (exponent == 2)
    denominator = paste0(denominator, "^2")
  parsed$args[2] = paste0("1/", denominator)
  parsed
}


# JAGS precision is scale^-exponent for these distributions.
jags_scale_exponent = function(dist) {
  switch(dist, dnorm = 2, dt = 2, dlnorm = 2, ddexp = 1, dlogis = 1, NULL)
}


#' Transform an mcp prior to the parameterization used by JAGS
#'
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated. JAGS uses precision rather than SD for some
#' distributions. For example, this function converts
#' `dnorm(4.2, 1.3)` into `dnorm(4.2, 1/1.3^2)`. It allows users to specify
#' priors using conventional scale parameters before they are translated to
#' JAGS code. Users normally do not need to call this function.
#'
#' @aliases sd_to_prec
#' @param prior_str String. A JAGS prior. Can be truncated, e.g.
#'   `dt(3, 2, 1) T(my_var, )`.
#' @return A string
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @encoding UTF-8
#' @export
sd_to_prec = function(prior_str) {
  lifecycle::deprecate_soft("0.4.0", "sd_to_prec()")
  format_prior(prior_to_jags(parse_prior(prior_str)))
}
