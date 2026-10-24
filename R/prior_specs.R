# Prior specifications -----------------------------------------------------

# Default builders return symbolic specifications only. Defaults and user
# priors are overlaid first and compiled together in one final pass.

empty_prior_specs = function() {
  tibble::tibble(
    parameter = character(),
    code = character(),
    description = character(),
    source = character()
  )
}


default_arma_specs = function() {
  tibble::tribble(
    ~dpar, ~par_type, ~prior, ~description, ~condition,
    "ar", "Intercept", "dnorm(0, 0.5) T(-1, 1)", "Zero-centered regularizing dependence coefficient", "always",
    "ar", "dummy", "dnorm(0, 0.25)", "Modest categorical change in a dependence coefficient", "always",
    "ar", "slope", "dnorm(0, 0.25 / predictor_scale())", "Modest dependence-coefficient change per predictor SD", "always",
    "ma", "Intercept", "dnorm(0, 0.5) T(-1, 1)", "Zero-centered regularizing dependence coefficient", "always",
    "ma", "dummy", "dnorm(0, 0.25)", "Modest categorical change in a dependence coefficient", "always",
    "ma", "slope", "dnorm(0, 0.25 / predictor_scale())", "Modest dependence-coefficient change per predictor SD", "always"
  )
}


# Default change point priors: uniform within the observed range, and for group-level
# change points, the SD and truncated deviations
default_cp_specs = function(cps, context) {
  n_cp = context$n_cp
  if (n_cp == 0)
    return(empty_prior_specs())

  prior_specs = list()
  for (j in seq_len(nrow(cps))) {
    name = cps$name[j]
    prior_specs[[name]] = tibble::tibble(
      parameter = name,
      code = "dunif(min(.x), max(.x))",
      description = paste0("Uniform within the observed range of ", context$x_display),
      source = "default"
    )
  }

  # Group-level change points
  for (j in seq_len(nrow(cps))) {
    if (is.na(cps$group_col[j]))
      next

    sd_name = cps$sd_name[j]
    group_name = cps$group_name[j]
    prior_specs[[sd_name]] = tibble::tibble(
      parameter = sd_name,
      code = "dnorm(0, (max(.x) - min(.x)) / (2 * (n_cp() + 1))) T(0, )",
      description = "Group-level change-point SD before truncation",
      source = "default"
    )

    prior_specs[[group_name]] = tibble::tibble(
      parameter = group_name,
      code = paste0("dnorm(0, ", sd_name, ") T(", paste(group_cp_truncation(cps, j, "min(.x)", "max(.x)"), collapse = ", "), ")"),
      description = "Group-level change-point deviations, truncated between the neighboring change points",
      source = "default"
    )
  }

  dplyr::bind_rows(prior_specs)
}


# Standard deviation of one model-matrix column for autoscaled priors (rstanarm).
# x-terms use the change-point variable: sd(x) and sd((x - min(x))^p) for powers.
default_predictor_scale = function(matrix_data, x_degree) {
  values = stats::na.omit(as.numeric(matrix_data))
  n_unique = length(unique(values))
  # For x-terms, a 0/1 indicator (as in x:state) only selects rows, so the column sd is close to sd(x).
  # Other two-valued columns, e.g., z in {100, 200}, carry units and are scaled.
  is_indicator = all(values %in% c(0, 1))
  data_scale = if (n_unique <= 1 || (is_indicator && x_degree > 0)) 1 else stats::sd(values)

  parts = character()
  if (x_degree > 0)
    parts = if (x_degree == 1) "sd(.x)" else paste0("sd((.x - min(.x))^", x_degree, ")")
  if (data_scale != 1)
    parts = c(parts, format_prior_number(data_scale))
  if (length(parts) == 0)
    parts = "1"
  paste0("(", paste(parts, collapse = " * "), ")")
}


# Replace `.offset` in default mu priors with the segment's offset symbol (see prior_context()),
# so priors are calibrated on the offset-adjusted response, e.g., mean(y - offset) or
# log(y) - offset. Without an offset, `.offset` terms are removed.
fill_offset = function(prior, offset) {
  no_offset = gsub("-min(.offset)", "0", prior, fixed = TRUE)  # e.g., T(-min(.offset), ) --> T(0, )
  no_offset = gsub(" - .offset", "", no_offset, fixed = TRUE)  # e.g., mean(.y - .offset) --> mean(.y)
  ifelse(is.na(offset), no_offset, stringr::str_replace_all(prior, stringr::fixed(".offset"), offset))
}


# Keep default-prior rows whose `condition` applies: "always", or "modeled"/"constant"
# depending on whether the model includes a formula for that dpar, e.g., sigma(1 + x).
filter_default_conditions = function(defaults, family) {
  for (dpar in unique(defaults$dpar[defaults$condition == "modeled"])) {
    condition = if (get_dpar_spec(family, dpar)$modeled) "modeled" else "constant"
    defaults = defaults[defaults$dpar != dpar | defaults$condition %in% c("always", condition), ]
  }
  defaults
}


# Replace predictor_scale() in default prior strings with each coefficient's scale.
fill_predictor_scale = function(prior, matrix_data, x_degree) {
  vapply(seq_along(prior), function(i) {
    gsub("predictor_scale()", default_predictor_scale(matrix_data[[i]], x_degree[i]), prior[i], fixed = TRUE)
  }, character(1))
}


# The par_type ("Intercept", "dummy", or "slope") whose default prior each row of `coefs`
# gets. Usually its own, but a factor in a segment that starts at zero, i.e., without an earlier
# intercept or x-term to join, is coded as all intercepts, e.g., `y ~ 0 + state` or segment 2 in
# `list(y ~ 0, ~ 0 + state)`. In other joined segments, `0 + state` gives jumps, which keep the contrast prior.
get_prior_type = function(coefs) {
  key = paste(coefs$dpar, coefs$order)
  joins = coefs$par_type == "Intercept" | coefs$x_degree > 0
  starts_at_zero = vapply(seq_len(nrow(coefs)), function(i) {
    !any(joins & key == key[i] & coefs$segment < coefs$segment[i])
  }, logical(1))

  term = paste(key, coefs$segment, coefs$term_key)
  sums_to_one = vapply(split(coefs$matrix_data, term), function(columns) {
    all(abs(Reduce(`+`, columns) - 1) < 1e-8, na.rm = TRUE)
  }, logical(1))

  is_factor_intercept = coefs$par_type == "dummy" & coefs$x_degree == 0 & starts_at_zero & sums_to_one[term]
  replace(coefs$par_type, is_factor_intercept, "Intercept")
}


# Default priors for population-level coefficients, looked up by dpar and par_type
default_coef_specs = function(coefs, family, context) {
  defaults = filter_default_conditions(dplyr::bind_rows(family$default_prior, default_arma_specs()), family)

  keys = paste(defaults$dpar, defaults$par_type)
  if (anyDuplicated(keys))
    stop_github("Default prior specifications are not unique by dpar and par_type.")

  joined = coefs |>
    dplyr::left_join(defaults, by = c("dpar", prior_type = "par_type"))
  is_factor_intercept = joined$prior_type == "Intercept" & joined$par_type != "Intercept"
  joined$description[is_factor_intercept] = paste0(joined$description[is_factor_intercept], "; one per factor level")
  if (any(is.na(joined$prior))) {
    stop_github(
      "mcp could not find a default prior for ",
      and_collapse(joined$name[is.na(joined$prior)])
    )
  }

  # Calibrate mu priors on the offset-adjusted response. Only mu offsets are used.
  offset = ifelse(joined$dpar == "mu", context$offsets[as.character(joined$segment)], NA_character_)
  is_adjusted = !is.na(offset) & grepl(".offset", joined$prior, fixed = TRUE)
  joined$description[is_adjusted] = paste0(joined$description[is_adjusted], "; adjusted for the offset")
  joined$prior = fill_offset(joined$prior, offset)

  joined$prior = fill_predictor_scale(joined$prior, joined$matrix_data, joined$x_degree)

  tibble::tibble(
    parameter = joined$name,
    code = joined$prior,
    description = joined$description,
    source = "default"
  )
}


# Default priors for group-level SDs and deviations
default_group_specs = function(group_coefs, family) {
  if (nrow(group_coefs) == 0)
    return(empty_prior_specs())

  defaults = filter_default_conditions(family$default_prior, family)
  joined = group_coefs |>
    dplyr::left_join(
      dplyr::select(defaults, "dpar", "par_type", "group_sd_prior"),
      by = c("dpar", "par_type")
    )
  if (any(is.na(joined$group_sd_prior))) {
    stop_github(
      "mcp could not find a default group-level SD prior for ",
      and_collapse(joined$name[is.na(joined$group_sd_prior)])
    )
  }

  joined$group_sd_prior = fill_predictor_scale(joined$group_sd_prior, joined$matrix_data, joined$x_degree)
  coefficient_description = ifelse(
    joined$par_type == "Intercept",
    "intercept",
    paste0("coefficient `", joined$matrix_name, "`")
  )

  dplyr::bind_rows(
    tibble::tibble(
      parameter = joined$sd_name,
      code = joined$group_sd_prior,
      description = paste0(
        "SD of group-level ", joined$dpar, " ",
        coefficient_description, " deviations"
      ),
      source = "default"
    ),
    tibble::tibble(
      parameter = joined$name,
      code = paste0("dnorm(0, ", joined$sd_name, ")"),
      description = paste0(
        "Zero-mean group-level ", joined$dpar, " ",
        coefficient_description, " deviations"
      ),
      source = "default"
    )
  )
}


# Truncation bounds keeping group-level change point j between the previous realized change point
# (or `x_min`) and the next population-level change point (or `x_max`), as deviations from cp_j.
group_cp_truncation = function(cps, j, x_min, x_max) {
  previous_cp = if (j == 1) {
    x_min
  } else if (is.na(cps$group_name[j - 1])) {
    cps$name[j - 1]
  } else {
    paste0("(", cps$name[j - 1], " + ", cps$group_name[j - 1], "[", cps$group_col[j], "_])")
  }
  later_population = which(seq_len(nrow(cps)) > j & is.na(cps$group_col))
  next_cp = if (length(later_population) == 0) x_max else cps$name[later_population[1]]
  paste0(c(previous_cp, next_cp), " - ", cps$name[j])
}


# Minimum length of the segment that ends at change point j, written as the lower
# bound `cp_1 + 20` in dunif(cp_1 + 20, max(x)) or T(cp_1 + 20, ). Returns the length as a
# string. "0" for just `cp_1`, or NULL. Change points are ordered anyway, so other references
# to change points in the bounds are errors.
get_min_segment_length = function(bounds, j, cp_names) {
  bounds = lapply(bounds, function(bound) if (nzchar(bound)) str2lang(bound))
  previous = if (j > 1) as.name(cp_names[j - 1])
  lower = bounds[[1]]
  min_length = NULL
  if (!is.null(previous) && identical(lower, previous)) {
    min_length = "0"
    bounds[[1]] = NULL
  } else if (!is.null(previous) && is.call(lower) && identical(lower[[1]], as.name("+")) && identical(lower[[2]], previous)) {
    min_length = deparse_prior(lower[[3]])
    bounds[[1]] = lower[[3]]
  }

  references = intersect(unlist(lapply(bounds, all.vars)), cp_names)
  if (length(references) > 0) {
    stop(
      "The lower or upper bound of the prior on `", cp_names[j], "` (in dunif() or T()) refers to ", and_collapse(paste0("`", references, "`")),
      ". Change points are always ordered, so the bounds can only refer to the previous change point to set a minimum segment length, ",
      "e.g., `cp_2 = \"dunif(cp_1 + 20, max(x))\"` or `cp_2 = \"dnorm(50, 10) T(cp_1 + 20, )\"`.",
      call. = FALSE
    )
  }
  min_length
}


# User priors on change points are used unchanged, restricted to ordered change points (see
# get_jags_code()). Priors without dunif() or T() are truncated to the observed range of x.
# An empty upper bound next to a minimum segment length, e.g., T(cp_1 + 20, ), is also max(x).
truncate_cp_prior = function(cps, j, parsed, context) {
  if (is.na(parsed$dist) || parsed$dist == "dirichlet")
    return(parsed)

  x_range = paste0(c("min(", "max("), context$x_display, ")")
  bounds = get_prior_bounds(parsed)
  if (!any(nzchar(bounds)))
    return(set_prior_bounds(parsed, x_range))

  min_length = get_min_segment_length(bounds, j, cps$name)
  if (identical(min_length, "0")) {
    warning(
      "`", cps$name[j], "` = \"", format_prior(parsed), "\": `", cps$name[j - 1], "` as the lower bound is already embodied in the mcp model and hence not needed. ",
      "Each change point keeps its specified prior, restricted to ordered change points. ",
      "To set a minimum segment length, use e.g. `", cps$name[j - 1], " + 20` as the lower bound.",
      call. = FALSE
    )
  }
  if (!is.null(min_length) && !nzchar(bounds[2]))
    return(set_prior_bounds(parsed, c(bounds[1], x_range[2])))
  parsed
}


# Bound a parsed prior below at `lower` unless it is already bounded.
truncate_lower_prior = function(parsed, lower, name = NULL) {
  # Only stochastic distributions without existing truncation need bounding.
  # dinvgamma() is positive and written as a transformed dgamma(), which cannot take T().
  if (is.na(parsed$dist) || any(nzchar(parsed$trunc)) || parsed$dist == "dinvgamma")
    return(parsed)

  # dunif cannot take T() in JAGS; bound through its lower parameter instead
  if (parsed$dist == "dunif" && length(parsed$args) == 2) {
    lower_val = suppressWarnings(as.numeric(parsed$args[1]))
    upper_val = suppressWarnings(as.numeric(parsed$args[2]))
    if (!is.na(upper_val) && upper_val <= lower)
      stop("Prior", if (!is.null(name)) paste0(" for '", name, "'"), " must allow values above ", lower, ".")
    if (!is.na(lower_val) && lower_val < lower)
      parsed$args[1] = format_prior_number(lower)
    return(parsed)
  }

  # Truncate unbounded distributions at the lower bound
  parsed$trunc = c(format_prior_number(lower), "")
  parsed
}


# Check that priors are set on the original parameter, not on names that same() replaced.
# - coefs, group_coefs: occurrence tables, i.e., including occurrences from same().
assert_prior_on_definitions = function(prior, coefs, group_coefs) {
  occurrences = dplyr::bind_rows(
    dplyr::transmute(
      dplyr::filter(coefs, .data$segment != .data$definition_segment),
      replaced = stringr::str_replace(.data$name, "_[0-9]+$", paste0("_", .data$segment)),
      original = .data$name
    ),
    dplyr::transmute(
      dplyr::filter(group_coefs, .data$segment != .data$definition_segment),
      replaced = stringr::str_replace(.data$sd_name, paste0("_", .data$definition_segment, "_"), paste0("_", .data$segment, "_")),
      original = .data$sd_name
    )
  ) |> dplyr::distinct()
  bad_prior = occurrences[occurrences$replaced %in% names(prior), , drop = FALSE]
  if (nrow(bad_prior) > 0) {
    guidance = if (nrow(bad_prior) == 1) {
      paste0("`", bad_prior$original, "` instead of `", bad_prior$replaced, "`")
    } else {
      paste0("`", bad_prior$original, "` (not `", bad_prior$replaced, "`)")
    }
    stop("`same()` reuses earlier parameters instead of creating new ones. Set prior on ",
      and_collapse(guidance), ".", call. = FALSE)
  }
  invisible(TRUE)
}


# Replace default prior specs with the user's parsed `prior`, adding the bounds that the model
# requires, e.g., the observed range for change points and sigma > 0. Updates descriptions accordingly.
overlay_user_prior_specs = function(prior_specs, prior, cps, context, coefs, family) {
  name_matches = names(prior) %in% prior_specs$parameter
  if (any(!name_matches)) {
    stop(
      "Prior(s) were specified for parameter name(s) that are not part of the model: ",
      and_collapse(names(prior)[!name_matches])
    )
  }

  # dirichlet(alpha) applies to all change points, so propagate it to unassigned change points
  user_cp_names = intersect(names(prior), cps$name)
  user_cp_is_dirichlet = vapply(prior[user_cp_names], function(parsed) parsed$dist %in% "dirichlet", logical(1))  # %in% is FALSE for dist = NA
  if (any(user_cp_is_dirichlet)) {
    if (!all(user_cp_is_dirichlet))
      stop("All or none of the change point priors must be `dirichlet(alpha)`.")
    source_desc = and_collapse(paste0("`", user_cp_names, "`"))
    for (u_name in setdiff(cps$name, user_cp_names)) {
      i = match(u_name, prior_specs$parameter)
      prior_specs$parsed[[i]] = prior[[user_cp_names[1]]]
      prior_specs$source[i] = "user"
      prior_specs$description[i] = paste0("Inherited from user-specified ", source_desc)
    }
  }

  # Bound change points within the observed range
  auto_truncated = character()
  for (j in seq_len(nrow(cps))) {
    name = cps$name[j]
    if (name %in% names(prior)) {
      original = prior[[name]]
      prior[[name]] = truncate_cp_prior(cps, j, original, context)
      if (!identical(prior[[name]], original))
        auto_truncated = c(auto_truncated, name)
    }

    group_name = cps$group_name[j]
    if (!is.na(cps$group_col[j]) && group_name %in% names(prior)) {
      parsed = prior[[group_name]]
      if (!is.na(parsed$dist) && !any(nzchar(get_prior_bounds(parsed)))) {
        x_display = context$x_display
        prior[[group_name]]$trunc = group_cp_truncation(cps, j, paste0("min(", x_display, ")"), paste0("max(", x_display, ")"))
        auto_truncated = c(auto_truncated, group_name)
      }
    }
  }

  # Enforce the support of constant dpars, e.g., sigma > 0
  bounds = get_constant_dpar_bounds(family)
  for (name in intersect(names(prior), coefs$name[coefs$dpar %in% names(bounds)])) {
    dpar = coefs$dpar[match(name, coefs$name)]
    original = prior[[name]]
    prior[[name]] = truncate_lower_prior(original, bounds[[dpar]], name)
    if (!identical(prior[[name]], original))
      auto_truncated = c(auto_truncated, name)
  }

  # Overlay
  for (name in names(prior)) {
    i = match(name, prior_specs$parameter)
    prior_specs$parsed[[i]] = prior[[name]]
    prior_specs$source[i] = "user"
    prior_specs$description[i] = if (name %in% auto_truncated && name %in% cps$name) {
      paste0("User-specified prior, truncated to the observed range of ", context$x_display)
    } else if (name %in% auto_truncated) {
      "User-specified prior with required bounds added by mcp"
    } else {
      NA_character_
    }
  }

  # Change points with a prior distribution are ordered (see get_jags_code())
  if (nrow(cps) > 1) {
    is_ordered = prior_specs$parameter %in% cps$name & vapply(prior_specs$parsed, function(parsed) !is.na(parsed$dist), logical(1))
    description = dplyr::coalesce(prior_specs$description[is_ordered], "User-specified prior")
    prior_specs$description[is_ordered] = paste0(description, ", restricted to ordered change points")
    if (length(user_cp_names) == 0)
      prior_specs$description[is_ordered] = "Flat Dirichlet on the segment lengths"
  }
  prior_specs
}


# Segment-1 x-terms that change the coefficient `name` along x. These have the same
# model-matrix column apart from x, e.g., `stateb.x` for `stateb` in `1 + state * x`.
# Intercepts (incl. each level in `0 + state`) also have x-terms without other variables, e.g., `x`.
get_x_partners = function(coefs, name) {
  i = match(name, coefs$name)
  own = coefs$matrix_data[[i]]
  candidates = coefs[coefs$segment == 1 & coefs$dpar == coefs$dpar[i] & coefs$x_degree > 0, ]
  is_partner = vapply(candidates$matrix_data, function(values) {
    (coefs$prior_type[i] == "Intercept" && all(values == 1, na.rm = TRUE)) || isTRUE(all.equal(values, own))
  }, logical(1))
  candidates[is_partner, ]
}


# Place default segment-1 priors on the value at the segment start, min(x), like
# brms and rstanarm place intercept priors on centered predictors. The reported
# coefficients remain at x = 0 (as in lm()), so their priors are conditional on the
# coefficients' x-terms, e.g., "dnorm(mean(y) - x_1 * min(x), 2.5 * sd(y))" or
# "dnorm(0 - stateb.x_1 * min(x), ...)". Other predictors are at 0.
add_intercept_references = function(prior_specs, coefs, context) {
  if (context$x_min == 0)
    return(prior_specs)

  targets = coefs[coefs$segment == 1 & coefs$x_degree == 0 & coefs$dpar %notin% c("ar", "ma"), ]
  for (i in seq_len(nrow(targets))) {
    row = match(targets$name[i], prior_specs$parameter)
    if (is.na(row) || prior_specs$source[row] != "default")
      next
    parsed = prior_specs$parsed[[row]]
    if (parsed$dist %notin% c("dnorm", "dt", "dunif"))
      next

    terms = get_x_partners(coefs, targets$name[i])
    if (nrow(terms) == 0)
      next

    # Subtract the x-terms at min(x) from the location and the truncation bounds
    reference = paste0(" - ", terms$name, " * min(", context$x_display, ")", x_power(terms$x_degree), collapse = "")
    prior_specs$parsed[[row]] = shift_prior(parsed, function(arg) if (arg == "0") reference else paste0(arg, reference))
    what = if (targets$prior_type[i] == "Intercept") "value" else "difference"
    prior_specs$description[row] = paste0(prior_specs$description[row], "; on the ", what, " at min(", context$x_display, ")")
  }
  prior_specs
}
