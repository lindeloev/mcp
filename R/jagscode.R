# ABOUT: These functions "pad" the regression model from get_formula.R
# resulting in a full JAGS model
# -----------------

get_jags_family_context = function(segments) {
  list(
    y = paste0(segments$y[1], "[i_]"),
    boundary = "garma_boundary_[i_]",
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


#' Get (or create) a named JAGS data constant for a numeric value
#'
#' Deduplicates by `format_prior_number()`-formatted value so that the same
#' underlying quantity (e.g., `min(.x)`), however it is spelled out in the
#' model text, maps to a single constant. Critically, `value` is stored
#' as-is (full precision) rather than being re-parsed from formatted text,
#' since `format_prior_number()` is lossy (see `jagsify_constants()`).
#'
#' @keywords internal
#' @noRd
#' @param value A single raw numeric value.
#' @param registry An environment with `names` (formatted-value -> constant
#'   name), `values` (constant name -> raw value), and `counter`.
#' @return String. The constant's name in `registry`.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
register_jags_constant = function(value, registry) {
  key = format_prior_number(value)
  name = registry$names[[key]]
  if (is.null(name)) {
    registry$counter = registry$counter + 1L
    name = paste0("CONST", registry$counter, "_")
    registry$names[[key]] = name
    registry$values[[name]] = value
  }
  name
}


#' Replace numeric literals with named JAGS data constants
#'
#' JAGS's adaptive slice sampler for a change point prior can get stuck very
#' close to its initial value when the prior's `dunif()` bounds or `T()`
#' truncation bounds are literal numbers in the model code - even though the
#' value is numerically identical to a `data`-supplied constant. Routing
#' these numbers through named `data` constants instead avoids the problem.
#'
#' Numeric literals in `x` are themselves `format_prior_number()`-formatted
#' text, which is lossy (fewer significant digits than the underlying
#' double). This is harmless for display, but for a change point's `dunif()`
#' bounds - which double as indicator thresholds in the likelihood (e.g.
#' `x[i_] >= cp_0`) - losing precision can round a bound past the data's
#' true min/max, excluding that observation and destabilizing the sampler.
#' Callers should `register_jags_constant()` such values with their *raw*
#' (unformatted) value beforehand, so this function's regex match reuses
#' the already-registered exact value instead of re-parsing lossy text.
#'
#' @aliases jagsify_constants
#' @keywords internal
#' @noRd
#' @param x String with a (partly) resolved prior, e.g., `"dt(30.5, 2.4, 1) T(cp_1, 99.2)"`.
#' @param registry An environment with `names` (value -> constant name),
#'   `values` (constant name -> value), and `counter`.
#' @return `x` with numeric literals replaced by names of constants added to `registry`.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
jagsify_constants = function(x, registry) {
  numeric_pattern = "(?<![A-Za-z0-9_.])-?[0-9]+\\.?[0-9]*(?:[eE][-+]?[0-9]+)?"
  stringr::str_replace_all(x, numeric_pattern, function(literals) {
    vapply(literals, function(literal) {
      register_jags_constant(as.numeric(literal), registry)
    }, character(1), USE.NAMES = FALSE)
  })
}


#' Make JAGS code for Multiple Change Point model
#'
#' @aliases get_jags_code
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param formula_jags String. The formula string returned by `get_formula_jags()`.
#' @param segments Segment table returned by `get_segment_tables()`.
#' @param group_effects Group-effect table returned by `get_group_effects()`.
#' @param ar_order,ma_order NA or positive integer. The GARMA component orders.
#' @param prior_table Resolved prior metadata from `get_prior()`.
#' @param prior_context Data summaries used while resolving priors.
#' @return String. A JAGS model. Has attribute `jags_constants`: a named list
#'   of data constants referenced by change point priors (see `jagsify_constants()`).
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_jags_code = function(prior, segments, group_effects, formula_jags, ar_order, ma_order, family, par_x,
                          prior_table = NULL, prior_context = NULL, series = FALSE) {
  prior_description = if (is.null(prior_table)) {
    stats::setNames(rep("Prior", length(prior)), names(prior))
  } else {
    stats::setNames(prior_table$description, prior_table$parameter)
  }
  prior_kind_ = if (is.null(prior_table)) {
    NULL
  } else {
    stats::setNames(prior_table$kind, prior_table$parameter)
  }
  prior_context_ = prior_context

  # Change point priors are bounded/truncated, which is where a literal-number
  # JAGS sampler bug can strike (see jagsify_constants()). Route their
  # numbers through named data constants instead.
  jags_constants = new.env()
  jags_constants$names = list()
  jags_constants$values = list()
  jags_constants$counter = 0L
  if (!is.null(prior_context_)) {
    # Pre-register with the *raw* x_min/x_max so later jagsify_constants()
    # calls on the (lossy, formatted-text) cp priors reuse these exact
    # values instead of re-parsing a rounded literal (see docs above).
    register_jags_constant(prior_context_$x_min, jags_constants)
    register_jags_constant(prior_context_$x_max, jags_constants)
  }

  # Begin building JAGS model. `mm` is short for "mcp model".
  # Add fixed variables.
  mm = paste0("model {")

  ####################################
  # DIRICHLET PRIOR ON CHANGE POINTS #
  ####################################
  # Get change point priors and check if they are Dirichlet
  cps = prior[stringr::str_detect(names(prior), "^cp_[1-9]+$")]
  dirichlet_calls = lapply(cps, function(x) {
    call = parse_prior_call(x)
    if (!is.null(call) && call$name == "dirichlet") call else NULL
  })
  is_dirichlet = !vapply(dirichlet_calls, is.null, logical(1))
  if (any(is_dirichlet)) {
    if (!all(is_dirichlet))
      stop("All or none of the change point priors must be `dirichlet(alpha)`.")

    alpha = vapply(dirichlet_calls, function(call) {
      if (length(call$args) != 1)
        return(NA_real_)
      suppressWarnings(as.numeric(call$args))
    }, numeric(1))
    if (any(!is.finite(alpha)) || any(alpha <= 0))
      stop("`dirichlet(alpha)` requires one finite alpha > 0.")
    if (any(alpha != alpha[1]))
      stop("All `dirichlet(alpha)` change point priors must use the same alpha.")

    # Build JAGS code. cp_betas is a simplex. cp_i is scaled to the observed range of x.
    mm = paste0(mm, "
  # Scaled Dirichlet prior on change points
  cp_betas ~ ddirch(c(", paste(rep(format_prior_number(alpha[1]), length(cps) + 1L), collapse = ", "), "))  # Scaled Dirichlet prior on change points")
    for (i in seq_along(cps)) {
      mm = paste0(mm, "
  cp_", i, " = ", format_prior_number(prior_context_$x_min),
                  " + sum(cp_betas[1:", i, "]) * ", format_prior_number(prior_context_$x_span),
                  "  # Within the observed change-point span")
    }

    # Clean up. Remove any dirichlet priors from the list of priors
    prior[names(cps)] = NULL
  }

  ################
  # OTHER PRIORS #
  ################
  # ... also handles non-Dirichlet priors

  # Split scalar/population priors from group-indexed coefficient priors.
  prior_pop = prior[!names(prior) %in% group_effects$name]
  prior_group = prior[names(prior) %in% group_effects$name]

  # Use get_prior_str() to add population-level priors
  mm = paste0(mm, "
  # mcp helper values\n")

  # Helpers for change points:
  mm = paste0(mm, "  cp_0 = ", jagsify_constants(format_prior_number(prior_context_$x_min), jags_constants), "\n")
  mm = paste0(mm, "  cp_", max(segments$segment), " = ", jagsify_constants(format_prior_number(prior_context_$x_max), jags_constants), "

  # Priors for population-level effects\n")
  is_cp_name = function(name) grepl("^cp_[0-9]+$", name)
  for (i in seq_along(prior_pop)) {
    name = names(prior_pop)[i]
    if (is_cp_name(name))
      prior_pop[[i]] = jagsify_constants(prior_pop[[i]], jags_constants)
    mm = paste0(mm, get_prior_str(
      prior_pop, i,
      description = prior_description[[name]],
      kind = if (is.null(prior_kind_)) NULL else prior_kind_[[name]]
    ))
  }


  # Use get_prior_str() to add group-level priors.
  if (length(prior_group) > 0) {
    mm = paste0(mm, "\n  # Priors for group-level effects\n")
    for (i in 1:length(prior_group)) {
      effect = group_effects[group_effects$name == names(prior_group)[i], , drop = FALSE]
      if (nrow(effect) != 1)
        stop_github("Expected exactly one group-effect row for ", names(prior_group)[i], ".")
      prior_group[[i]] = jagsify_constants(prior_group[[i]], jags_constants)
      mm = paste0(mm, get_prior_str(
        prior = prior_group,
        i = i,
        group_col = effect$group_col,
        center = effect$part == "cp",
        description = prior_description[[names(prior_group)[i]]],
        kind = if (is.null(prior_kind_)) NULL else prior_kind_[[names(prior_group)[i]]]
      ))
    }
  }

  cp_effects = group_effects[group_effects$part == "cp", , drop = FALSE]
  if (nrow(cp_effects) > 0 && nrow(segments) > 2) {
    group_col = cp_effects$group_col[1]
    boundaries = segments$cp_code_form[-1]
    boundaries = gsub(
      "CP_[0-9]+_INDEX", paste0("[", group_col, "_]"), boundaries
    )
    mm = paste0(mm, "\n  # Order realized group-level change points\n",
      "  for (", group_col, "_ in 1:n_unique_", group_col, ") {\n")
    for (j in seq_len(length(boundaries) - 1L)) {
      difference = paste0(boundaries[j + 1L], " - ", boundaries[j])
      mm = paste0(
        mm, "    cp_order_[", group_col, "_, ", j, "] ~ dbern(step(",
        difference, ") * (1 - equals(", difference, ", 0)))\n"
      )
    }
    mm = paste0(mm, "  }\n")
  }


  #########
  # GARMA #
  #########
  has_arma_terms = !is.na(ar_order) || !is.na(ma_order)
  if (has_arma_terms)
    mm = paste0(mm, get_garma_jagscode(ar_order, ma_order, par_x, series))



  ###########
  # FORMULA #
  ###########
  # Transform formula_jags into JAGS format. Insert par_x and group indices.
  for (i in seq_len(max(segments$segment))) {
    formula_jags = gsub(paste0("CP_", i, "_INDEX"), paste0("[", segments$cp_group_col[i], "[i_]]"), formula_jags)
  }

  # Insert formula_jags
  mm = paste0(mm, "
  # Model and likelihood
  for (i_ in 1:length(", par_x, ")) {")

  # Add JAGS code for fitted values and indent it
  mm = paste0(mm, gsub("\n", "\n    ", formula_jags))


  ##############
  # LIKELIHOOD #
  ##############
  context = get_jags_family_context(segments)
  mm = paste0(mm, "\n\n    # Likelihood and log-density for family = ", family$family, "()
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
    if (!is.na(spec$lower))
      response_code = paste0("max(", format(spec$lower, scientific = TRUE), ", ", response_code, ")")

    mm = paste0(mm, dpar, "_[i_] = ", response_code, "\n    ")
  }

  likelihood = family$backends$jags$likelihood(context)
  mm = paste0(mm, paste(likelihood, collapse = "\n    "))

  # Compute link-scale residuals for GARMA
  if (has_arma_terms) {
    garma_y = family$garma$observed_jags(context)
    garma_link_y = if (family$linkfun_str == "") "garma_y_[i_]" else
      paste0(family$linkfun_str, "(garma_y_[i_])")
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
  mm = paste0(mm, "
  }
}")
  attr(mm, "jags_constants") = jags_constants$values
  mm
}


#' Get JAGS code for a prior
#'
#' @aliases get_prior_str
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param i The index in `prior` to get code for
#' @param group_col String or NULL. `NULL` indicates a population-level prior.
#'   A string indicates a group-level prior (one value for each group
#'   level).
#' @param center Logical. Exactly zero-center a group-indexed vector?
#' @param description Short comment to include in generated JAGS code.
#' @param kind One of distribution, alias, expression, or constant.
#' @return A string
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @encoding UTF-8
get_prior_str = function(prior, i, group_col = NULL, center = FALSE,
                          description = "Prior", kind = NULL) {
  # Helpers
  value = prior[[i]]
  name = names(prior[i])
  description = gsub("[\r\n]+", " ", description)
  if (is.null(kind))
    kind = prior_kind(value, names(prior))

  # JAGS does not support dinvgamma or dloginvgamma. Write it manually.
  # If inverse_shape ~ Gamma(a, b), then
  # shape = exp(-log(inverse_shape)) ~ InvGamma(a, b).
  if (stringr::str_detect(value, "^dloginvgamma\\(")) {
    if (!is.null(group_col))
      stop("dloginvgamma() is currently only supported for population-level parameters.")

    call = parse_prior_call(value)
    if (is.null(call) || call$name != "dloginvgamma" || length(call$args) != 2)
      stop("Expected dloginvgamma(shape, scale). Got '", value, "'.")

    inverse_name = paste0(name, "_inverse")
    return(paste0(
      "  ", inverse_name, " ~ dgamma(", call$args[1], ", ", call$args[2], ")  # ", description, "\n",
      "  ", name, " = -log(", inverse_name, ")  # Log of inverse-gamma parameter\n"
    ))
  }

  if (kind == "distribution") {
    # Convert to precision
    value = sd_to_prec(value)

    # ... and this is a population-level effect
    if (is.null(group_col)) {
      return(paste0("  ", name, " ~ ", value, "  # ", description, "\n"))
    } else if (!center) {
      return(paste0("  for (", group_col, "_ in 1:n_unique_", group_col, ") {
    ", name, "[", group_col, "_] ~ ", value, "  # ", description, "
  }\n"))
    } else {
      return(paste0("  for (", group_col, "_ in 1:n_unique_", group_col, ") {
    ", name, "_uncentered[", group_col, "_] ~ ", value, "  # ", description, "
  }
  ", name, " = ", name, "_uncentered - mean(", name, "_uncentered)  # vectorized zero-centering\n"))
    }
  }

  paste0("  ", name, " = ", value, "  # ", description, "\n")
}


#' Transform a JAGS Prior from SD to Precision.
#'
#' JAGS uses precision rather than SD. This function converts
#' `dnorm(4.2, 1.3)` into `dnorm(4.2, 1/1.3^2)`. It allows users to specify
#' priors using SD and then it's transformed for the JAGS code. It works for the
#' following distributions: dnorm|dt|dcauchy|ddexp|dlogis|dlnorm. In all of
#' these,
#' tau/sd is the second parameter.
#'
#' @aliases sd_to_prec
#' @param prior_str String. A JAGS prior. Can be truncated, e.g.
#'   `dt(3, 2, 1) T(my_var, )`.
#' @return A string
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @encoding UTF-8
#' @export
sd_to_prec = function(prior_str) {
  parts = split_prior_truncation(prior_str)
  call = parse_prior_call(parts$distribution)
  if (is.null(call) || call$name %notin% c("dnorm", "dt", "dcauchy", "ddexp", "dlogis", "dlnorm"))
    return(prior_str)
  if (length(call$args) < 2)
    stop("Expected a scale as the second argument of '", prior_str, "'.")

  call$args[2] = paste0("1/(", gsub(" ", "", call$args[2]), ")^2")
  converted = paste0(call$name, "(", paste(call$args, collapse = ", "), ") ")
  if (!is.null(parts$truncation))
    converted = paste0(converted, gsub("\\s+", "", parts$truncation))
  converted
}
