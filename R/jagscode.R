# ABOUT: These functions "pad" the linear predictors from get_formula_jags()
# resulting in a full JAGS model
# -----------------

# JAGS names for the per-observation nodes used by the family's likelihood code,
# e.g., the response `response[i_]`, distributional parameters `sigma_[i_]`,
# and auxiliary response columns such as `trials(N)` as `N[i_]`.
get_jags_family_context = function(segments) {
  list(
    y = paste0(segments$y[1], "[i_]"),
    threshold = "garma_threshold_[i_]",
    dpar = function(name) paste0(name, "_[i_]"),
    local = function(name) paste0(name, "_[i_]"),
    aux = function(name, default = NULL) {
      if (name %notin% names(segments) || all(is.na(segments[[name]]))) {
        if (!is.null(default))
          return(default)
        stop_github("Missing response auxiliary ", name, "() in JAGS family context.")
      }
      paste0(stats::na.omit(unique(segments[[name]]))[1], "[i_]")
    }
  )
}


#' Make JAGS code for Multiple Change Point model
#'
#' @param tables Model tables, as in `get_fit_model_tables(fit)`.
#' @param prior Resolved priors from `get_prior()`.
#' @param prior_table Resolved prior metadata from `get_prior()`.
#' @param family An mcpfamily object.
#' @return String. A JAGS model.
#' @noRd
get_jags_code = function(tables, prior, prior_table, family) {
  segments = tables$segments
  group_effects = definition_rows(tables$group_effects)
  par_x = tables$data_columns$par_x
  ar_order = get_arma_order(tables$predictors, "ar")
  ma_order = get_arma_order(tables$predictors, "ma")
  prior_description = stats::setNames(prior_table$description, prior_table$parameter)
  prior_kind_ = stats::setNames(prior_table$kind, prior_table$parameter)

  # Parameters sampled via an auxiliary node (see add_reparams())
  aux_table = prior_table[!is.na(prior_table$aux_suffix), , drop = FALSE]

  # Begin building JAGS model. `mm` is short for "mcp model".
  # Change points come first and in order, from the fixed cp_0 = min(x) to cp_K = max(x).
  K = max(segments$segment)
  mm = paste0("model {
  # Change points
  cp_0 = min(", par_x, ")\n")

  #######################
  # CHANGE POINT PRIORS #
  #######################
  # Each change point keeps its specified prior, and cp_order_ restricts them to ordered
  # change points with any minimum segment lengths. dirichlet(alpha) is the default uniform prior
  # expanded to the segment length.
  cps = prior[intersect(paste0("cp_", seq_len(K - 1L)), names(prior))]
  is_dirichlet = grepl("^dirichlet\\(", unlist(cps))
  if (any(is_dirichlet)) {
    # overlay_user_prior_specs() ensures that all or none are dirichlet()
    dirch_alpha = vapply(cps, function(x) {
      args = parse_prior_call(x)$args
      if (length(args) == 1) suppressWarnings(as.numeric(args)) else NA_real_
    }, numeric(1))
    if (any(!is.finite(dirch_alpha)) || any(dirch_alpha < 1))
      stop("`dirichlet(alpha)` requires one finite alpha >= 1.")
    if (any(dirch_alpha != dirch_alpha[1]))
      stop("All `dirichlet(alpha)` change point priors must use the same alpha.")
  }

  cp_prior_str = function(name, value = prior[[name]], ...) {
    get_prior_str(
      stats::setNames(list(value), name), 1,
      description = prior_description[[name]],
      kind = prior_kind_[[name]], ...
    )
  }
  cp_effects = group_effects[group_effects$part == "cp", , drop = FALSE]
  locations = cp_effects$name[vapply(prior[cp_effects$name], is_cp_location, logical(1))]
  min_lengths = rep("0", length(cps))
  for (i in seq_along(cps)) {
    name = names(cps)[i]
    value = cps[[i]]
    if (any(is_dirichlet)) {
      value = paste0("dunif(cp_0, cp_", K, ")")
    } else if (prior_kind_[[name]] == "distribution") {
      # A minimum segment length goes in cp_order_
      cp_prior = split_cp_bounds(value)
      min_length = get_min_segment_length(cp_prior$bounds, i, names(cps))
      if (!is.null(min_length)) {
        min_lengths[i] = min_length
        value = cp_prior$rebuild(c("cp_0", cp_prior$bounds[2]))
      }
    }
    mm = paste0(mm, cp_prior_str(name, value))

    # Group-level change points right after their population change point. Bounds at earlier
    # group-specific change points, e.g., "(cp_1 + cp_1_id[id_])", use their sampled locations.
    effect = cp_effects[cp_effects$population_name == name, , drop = FALSE]
    if (nrow(effect) == 1) {
      group_value = prior[[effect$name]]
      for (location in locations) {
        population_name = cp_effects$population_name[cp_effects$name == location]
        group_value = gsub(paste0(population_name, " + ", location, "["), paste0(location, "_location["), group_value, fixed = TRUE)
      }
      mm = paste0(mm, cp_prior_str(effect$sd_name), cp_prior_str(effect$name, group_value, group_col = effect$group_col, population_name = name))
    }
  }
  mm = paste0(mm, "  cp_", K, " = max(", par_x, ")\n")

  # Observed nodes (data: cp_order_ = 1 and cp_weight_ = 0) that restrict the priors above
  if (length(cps) > 1) {
    steps = paste0("step(cp_", 2:length(cps), " - cp_", 2:length(cps) - 1, ifelse(min_lengths[-1] == "0", "", paste0(" - ", min_lengths[-1])), ")")
    mm = paste0(mm, "  cp_order_ ~ dbern(", paste(steps, collapse = " * "), ")  # Ordered change points", if (any(min_lengths != "0")) " with minimum segment lengths", "\n")
  }
  if (any(is_dirichlet) && dirch_alpha[1] != 1) {
    # Zeros trick: dpois(0 | lambda) adds -lambda = (alpha - 1) * sum(log(segment length / range of x)) to the log density
    log_lengths = paste0("log((cp_", seq_len(K), " - cp_", seq_len(K) - 1, ") / (cp_", K, " - cp_0))", collapse = " + ")
    mm = paste0(mm, "  cp_weight_ ~ dpois(-", format_prior_number(dirch_alpha[1] - 1), " * (", log_lengths, "))  # dirichlet(", format_prior_number(dirch_alpha[1]), ") weight on the segment lengths\n")
  }
  mm = paste0(mm, "
  # Priors for population-level effects\n")

  ################
  # OTHER PRIORS #
  ################
  # Split scalar/population priors from group-indexed coefficient priors.
  is_cp = names(prior) %in% c(names(cps), cp_effects$name, cp_effects$sd_name)
  prior_pop = prior[!is_cp & !names(prior) %in% group_effects$name]
  prior_group = prior[!is_cp & names(prior) %in% group_effects$name]

  # Use get_prior_str() to add remaining population-level priors
  for (i in seq_along(prior_pop)) {
    name = names(prior_pop)[i]
    aux = aux_table[aux_table$parameter == name, , drop = FALSE]
    if (nrow(aux) == 1) {
      mm = paste0(mm, get_aux_prior_str(name, prior_pop[[i]], prior_description[[name]], aux))
    } else {
      mm = paste0(mm, get_prior_str(prior_pop, i, description = prior_description[[name]], kind = prior_kind_[[name]]))
    }
  }

  # Use get_prior_str() to add group-level priors.
  if (length(prior_group) > 0) {
    mm = paste0(mm, "\n  # Priors for group-level effects\n")
    for (i in seq_along(prior_group)) {
      effect = group_effects[group_effects$name == names(prior_group)[i], , drop = FALSE]
      if (nrow(effect) != 1)
        stop_github("Expected exactly one group-effect row for ", names(prior_group)[i], ".")
      mm = paste0(mm, get_prior_str(
        prior = prior_group,
        i = i,
        group_col = effect$group_col,
        description = prior_description[[names(prior_group)[i]]],
        kind = prior_kind_[[names(prior_group)[i]]]
      ))
    }
  }


  #########
  # GARMA #
  #########
  has_arma_terms = !is.na(ar_order) || !is.na(ma_order)
  if (has_arma_terms)
    mm = paste0(mm, get_garma_jagscode(ar_order, ma_order, par_x, series = !is.null(tables$data_columns$series)))


  ###########
  # FORMULA #
  ###########
  mm = paste0(mm, "
  # Model and likelihood
  for (i_ in 1:length(", par_x, ")) {")

  # Add JAGS code for fitted values and indent it
  mm = paste0(mm, gsub("\n", "\n    ", get_formula_jags(tables, family, locations)))


  ##############
  # LIKELIHOOD #
  ##############
  context = get_jags_family_context(segments)
  mm = paste0(mm, "\n\n    # Likelihood for family = ", family$family, "()
    ")

  # Transform link-scale predictors to distribution-scale parameters. These are
  # deterministic JAGS nodes and are not monitored by default.
  for (dpar in family$dpar_specs$dpar) {
    spec = get_dpar_spec(family, dpar)
    link_code = paste0("link_", dpar, "_[i_]")
    if (dpar == "mu" && has_arma_terms)
      link_code = paste0(link_code, " + resid_garma_[i_]")

    linkinv_str = get_link_str(spec$link, inverse = TRUE)
    response_code = ifelse(
      linkinv_str == "",
      link_code,
      paste0(linkinv_str, "(", link_code, ")")
    )
    if (!is.na(spec$floor))  # Numerical floor, as in evaluate_model()
      response_code = paste0("max(", format(spec$floor, scientific = TRUE), ", ", response_code, ")")

    mm = paste0(mm, dpar, "_[i_] = ", response_code, "\n    ")
  }

  likelihood = family$backends$jags$likelihood(context)
  mm = paste0(mm, paste(likelihood, collapse = "\n    "))

  # Compute link-scale residuals for GARMA
  if (has_arma_terms) {
    garma_y = family$garma$observed_jags(context)
    linkfun_str = get_link_str(family$link)
    garma_link_y = if (linkfun_str == "") "garma_y_[i_]" else
      paste0(linkfun_str, "(garma_y_[i_])")
    mm = paste0(
      mm,
      "\n    garma_y_[i_] = ", garma_y,
      "\n    garma_link_y_[i_] = ", garma_link_y
    )
    mm = paste0(mm, "\n    resid_abs_[i_] = garma_link_y_[i_] - link_mu_[i_]")
    if (!is.na(ma_order))
      mm = paste0(mm, "\n    resid_ma_[i_] = garma_link_y_[i_] - link_mu_[i_] - resid_garma_[i_]")
  }


  ###############
  # OTHER STUFF #
  ###############

  # Finish up
  paste0(mm, "
  }
}")
}


# Whether a group-level change point with this prior is sampled as its absolute location,
# e.g., cp_1_id_location ~ dnorm(cp_1 + mean, sd), rather than as a deviation (see get_prior_str()).
# The JAGS graph then uses the location directly; the deviation is a deterministic node.
is_cp_location = function(prior) {
  call = parse_prior_call(split_prior_truncation(prior)$distribution)
  !is.null(call) && call$name == "dnorm" && length(call$args) == 2
}


# Directly sampled parameters in `prior_table`, i.e., not constant, reparameterized (see add_reparams()),
# or sampled via a `_location` or `_inverse` node. Must agree with get_prior_str().
is_sampled_par = function(prior_table, parameters) {
  group_cps = parameters$parameter[parameters$part == "cp" & parameters$scope == "group"]
  is_location = prior_table$parameter %in% group_cps & vapply(prior_table$value, is_cp_location, logical(1))
  prior_table$kind == "distribution" & is.na(prior_table$aux_suffix) & !is_location & !grepl("^dinvgamma\\(", prior_table$value)
}


# Get the JAGS code (one or more lines) for the prior on `prior[[i]]`.
# - group_col: NULL for population-level priors. Otherwise the grouping column, and
#   the prior is written in a loop over its levels.
# - population_name: for group-level change points, the population change point, e.g., "cp_1".
# - description: comment to add in the JAGS code. kind: see prior_kind().
# Parameters sampled via another node must be excluded in is_sampled_par().
get_prior_str = function(prior, i, group_col = NULL, population_name = NULL,
                          description = "Prior", kind = NULL) {
  # Helpers
  value = prior[[i]]
  name = names(prior[i])
  description = gsub("[\r\n]+", " ", description)
  if (is.null(kind))
    kind = prior_kind(value, names(prior))

  # JAGS does not support dinvgamma. Write it manually.
  # If inverse_shape ~ Gamma(a, b), then shape = 1 / inverse_shape ~ InvGamma(a, b).
  if (stringr::str_detect(value, "^dinvgamma\\(")) {
    if (!is.null(group_col))
      stop("dinvgamma() is currently only supported for population-level parameters.")

    call = parse_prior_call(value)
    if (is.null(call) || call$name != "dinvgamma" || length(call$args) != 2)
      stop("Expected dinvgamma(shape, scale) without truncation. Got '", value, "'.")

    inverse_name = paste0(name, "_inverse")
    return(paste0(
      "  ", inverse_name, " ~ dgamma(", call$args[1], ", ", call$args[2], ")  # ", description, "\n",
      "  ", name, " = 1 / ", inverse_name, "  # Inverse-gamma parameter\n"
    ))
  }

  # Other distributions: convert SD/scale to JAGS precision and write a stochastic node
  if (kind == "distribution") {
    value = prior_to_jags(value)

    # Plain population-level parameter
    if (is.null(group_col)) {
      return(paste0("  ", name, " ~ ", value, "  # ", description, "\n"))

    # Ordinary group-level effect: one draw per group level
    } else if (is.null(population_name)) {
      return(paste0("  for (", group_col, "_ in 1:n_unique_", group_col, ") {
    ", name, "[", group_col, "_] ~ ", value, "  # ", description, "
  }\n"))

    # Group-level change point. For dnorm() priors, sample each group's absolute change point
    # location, cp_1_id_location ~ dnorm(cp_1 + mean, sd).
    } else {
      parts = split_prior_truncation(value)
      call = parse_prior_call(parts$distribution)
      if (is_cp_location(value)) {
        location_name = paste0(name, "_location")
        mean = if (call$args[1] == "0") population_name else paste0(population_name, " + ", call$args[1])
        distribution = paste0("dnorm(", mean, ", ", call$args[2], ")")
        truncation = ""
        if (!is.null(parts$truncation)) {
          trunc = parse_prior_call(parts$truncation)
          # Deviation bounds such as "cp_0 - cp_1" are simply "cp_0" on the location
          population_suffix = paste0(" - ", population_name)
          bounds = vapply(trunc$args, function(bound) {
            if (!nzchar(bound))
              return("")
            if (endsWith(bound, population_suffix))
              return(substr(bound, 1, nchar(bound) - nchar(population_suffix)))
            paste0(population_name, " + (", bound, ")")
          }, character(1))
          truncation = paste0(" T(", paste(bounds, collapse = ", "), ")")
        }
        return(paste0("  for (", group_col, "_ in 1:n_unique_", group_col, ") {
    ", location_name, "[", group_col, "_] ~ ", distribution, truncation, "  # ", description, "
    ", name, "[", group_col, "_] = ", location_name, "[", group_col, "_] - ", population_name, "  # deviation from population change point
  }\n"))
      }

      # Group-level for other distributions: sample the deviation directly
      return(paste0("  for (", group_col, "_ in 1:n_unique_", group_col, ") {
    ", name, "[", group_col, "_] ~ ", value, "  # ", description, "
  }\n"))
    }
  }

  # Fixed constants/expressions for group-level effects must be vectors over group
  # levels so they can be indexed as name[group_col_index_[i_]] in the likelihood.
  # References to other group-level effects are indexed too, e.g.,
  # "x_1_id * 2" becomes "x_1_id[id_] * 2".
  if (!is.null(group_col)) {
    value = gsub(paste0("([A-Za-z0-9_]+_", group_col, ")\\b(?!\\[)"), paste0("\\1[", group_col, "_]"), value, perl = TRUE)
    return(paste0("  for (", group_col, "_ in 1:n_unique_", group_col, ") {
    ", name, "[", group_col, "_] = ", value, "  # ", description, "
  }\n"))
  }

  # Merge it nicely and return
  paste0("  ", name, " = ", value, "  # ", description, "\n")
}


# Translate a scale-parameterized prior to JAGS ---------------------------

# mcp's legacy string syntax uses conventional distribution scales. JAGS
# uses inverse variance for normal-like distributions and inverse scale for
# double-exponential and logistic distributions.
prior_to_jags = function(prior_str) {
  parts = split_prior_truncation(prior_str)
  call = parse_prior_call(parts$distribution)
  if (is.null(call))
    return(prior_str)
  if (call$name == "dcauchy")
    stop("JAGS does not support `dcauchy()`; use `dt(location, scale, 1)` instead.")

  exponent = jags_scale_exponent(call$name)
  if (is.null(exponent))
    return(prior_str)
  if (length(call$args) < 2)
    stop("Expected a scale as the second argument of '", prior_str, "'.")

  denominator = paste0("(", gsub(" ", "", call$args[2]), ")")
  if (exponent == 2)
    denominator = paste0(denominator, "^2")
  call$args[2] = paste0("1/", denominator)
  converted = paste0(call$name, "(", paste(call$args, collapse = ", "), ") ")
  if (!is.null(parts$truncation))
    converted = paste0(converted, gsub("\\s+", "", parts$truncation))
  converted
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
  prior_to_jags(prior_str)
}
