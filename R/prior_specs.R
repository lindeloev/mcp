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


default_cp_specs = function(cps, context) {
  n_cp = context$n_cp
  if (n_cp == 0)
    return(empty_prior_specs())

  specs = list()
  for (j in seq_len(nrow(cps))) {
    name = cps$name[j]
    specs[[name]] = tibble::tibble(
      parameter = name,
      code = "dunif(min(.x), max(.x))",
      description = paste0("Uniform within the observed range of ", context$x_display),
      source = "default"
    )
  }

  for (j in seq_len(nrow(cps))) {
    if (is.na(cps$group_col[j]))
      next

    sd_name = cps$sd_name[j]
    group_name = cps$group_name[j]
    specs[[sd_name]] = tibble::tibble(
      parameter = sd_name,
      code = "dnorm(0, (max(.x) - min(.x)) / (2 * (n_cp() + 1))) T(0, )",
      description = "Group-level change-point SD before truncation",
      source = "default"
    )

    specs[[group_name]] = tibble::tibble(
      parameter = group_name,
      code = paste0("dnorm(0, ", sd_name, ") ", group_cp_truncation(cps, j, "min(.x)", "max(.x)")),
      description = "Group-level change-point deviations, truncated between the neighboring change points",
      source = "default"
    )
  }

  dplyr::bind_rows(specs)
}


# Standard deviation of one model-matrix column for autoscaled priors (rstanarm).
# x-terms use the change-point variable: sd(x) and sd((x - min(x))^p) for powers.
default_predictor_scale = function(matrix_data, x_factor) {
  values = stats::na.omit(as.numeric(matrix_data))
  n_unique = length(unique(values))
  # For x-terms, a 0/1 indicator (as in x:state) only selects rows, so the column sd is close to sd(x).
  # Other two-valued columns, e.g., z in {100, 200}, carry units and are scaled.
  is_indicator = all(values %in% c(0, 1))
  data_scale = if (n_unique <= 1 || (is_indicator && x_factor != "1")) 1 else stats::sd(values)

  parts = character()
  if (x_factor != "1") {
    degree = if (x_factor == "x") 1 else as.numeric(sub("x^", "", x_factor, fixed = TRUE))
    parts = if (degree == 1) "sd(.x)" else paste0("sd((.x - min(.x))^", degree, ")")
  }
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
fill_predictor_scale = function(prior, matrix_data, x_factor) {
  vapply(seq_along(prior), function(i) {
    gsub("predictor_scale()", default_predictor_scale(matrix_data[[i]], x_factor[i]), prior[i], fixed = TRUE)
  }, character(1))
}


# Which default prior each coefficient gets. Default priors are looked up by par_type
# ("Intercept", "dummy", or "slope"), and this returns the par_type to use for each row of
# `predictors`. It is the coefficient's own par_type, with one exception: a factor without an
# intercept (y ~ 0 + state) is coded as all intercepts. This mostly concerns segment 1,
# but applies wherever a segment starts at zero, i.e., has no earlier intercept or x-term to
# join, e.g., segment 2 in `list(y ~ 0, ~ 0 + state)`. In other joined segments, `0 + state`
# gives jumps from where the earlier segment left off, which keep the contrast prior.
get_prior_type = function(predictors) {
  key = paste(predictors$dpar, predictors$order)
  joins = predictors$par_type == "Intercept" | predictors$x_factor != "1"
  starts_at_zero = vapply(seq_len(nrow(predictors)), function(i) {
    !any(joins & key == key[i] & predictors$segment < predictors$segment[i])
  }, logical(1))

  term = paste(key, predictors$segment, predictors$term_key)
  sums_to_one = vapply(split(predictors$matrix_data, term), function(columns) {
    all(abs(Reduce(`+`, columns) - 1) < 1e-8, na.rm = TRUE)
  }, logical(1))

  is_factor_intercept = predictors$par_type == "dummy" & predictors$x_factor == "1" & starts_at_zero & sums_to_one[term]
  replace(predictors$par_type, is_factor_intercept, "Intercept")
}


default_predictor_specs = function(predictors, family, context) {
  defaults = filter_default_conditions(dplyr::bind_rows(family$default_prior, default_arma_specs()), family)

  keys = paste(defaults$dpar, defaults$par_type)
  if (anyDuplicated(keys))
    stop_github("Default prior specifications are not unique by dpar and par_type.")

  joined = predictors %>%
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

  joined$prior = fill_predictor_scale(joined$prior, joined$matrix_data, joined$x_factor)

  tibble::tibble(
    parameter = joined$name,
    code = joined$prior,
    description = joined$description,
    source = "default"
  )
}


default_group_specs = function(group_effects, family) {
  effects = group_effects %>%
    dplyr::filter(.data$part == "predictor")
  if (nrow(effects) == 0)
    return(empty_prior_specs())

  defaults = filter_default_conditions(family$default_prior, family)
  joined = effects %>%
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

  joined$group_sd_prior = fill_predictor_scale(joined$group_sd_prior, joined$matrix_data, joined$x_factor)
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


# Truncation keeping group-level change point j between the previous realized change point
# (or `x_min`) and the next population-level change point (or `x_max`), as deviations from cp_j.
group_cp_truncation = function(cps, j, x_min, x_max) {
  previous_cp = if (j == 1) x_min else sub("CP_[0-9]+_INDEX", paste0("[", cps$group_col[j], "_]"), cps$code[j - 1])
  later_population = which(seq_len(nrow(cps)) > j & is.na(cps$group_col))
  next_cp = if (length(later_population) == 0) x_max else cps$name[later_population[1]]
  paste0("T(", previous_cp, " - ", cps$name[j], ", ", next_cp, " - ", cps$name[j], ")")
}


# Only stochastic distributions can be truncated, not fixed values like "cp_1 + 5" or "cp_1_id"
is_distribution = function(prior_value) {
  call = if (is.character(prior_value)) parse_prior_call(split_prior_truncation(prior_value)$distribution)
  !is.null(call) && grepl("^d[A-Za-z]", call$name)
}


# The lower and upper bounds of a change point prior: the arguments of dunif() or T(), or NULL
# without them. `rebuild(bounds)` writes the prior with other T(low, high) bounds.
split_cp_bounds = function(prior_value) {
  parts = split_prior_truncation(prior_value)
  call = parse_prior_call(parts$distribution)
  if (call$name == "dunif") {
    if (!is.null(parts$truncation))
      stop("JAGS cannot truncate dunif(). Set its lower and upper bounds as its arguments instead, e.g., \"dunif(20, 50)\".", call. = FALSE)
    return(list(bounds = call$args, rebuild = function(bounds) paste0("dunif(", bounds[1], ", ", bounds[2], ")")))
  }
  bounds = if (!is.null(parts$truncation)) parse_prior_call(parts$truncation)$args
  list(bounds = bounds, rebuild = function(bounds) paste0(parts$distribution, " T(", bounds[1], ", ", bounds[2], ")"))
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
truncate_cp_prior = function(cps, j, prior_value, context) {
  if (!is_distribution(prior_value) || grepl("^\\s*dirichlet\\s*\\(", prior_value))
    return(prior_value)

  x_range = paste0(c("min(", "max("), context$x_display, ")")
  cp_prior = split_cp_bounds(prior_value)
  if (is.null(cp_prior$bounds))
    return(cp_prior$rebuild(x_range))

  min_length = get_min_segment_length(cp_prior$bounds, j, cps$name)
  if (identical(min_length, "0")) {
    warning(
      "`", cps$name[j], "` = \"", prior_value, "\": `", cps$name[j - 1], "` as the lower bound is already embodied in the mcp model and hence not needed. ",
      "Each change point keeps its specified prior, restricted to ordered change points. ",
      "To set a minimum segment length, use e.g. `", cps$name[j - 1], " + 20` as the lower bound.",
      call. = FALSE
    )
  }
  if (!is.null(min_length) && !nzchar(cp_prior$bounds[2]))
    return(cp_prior$rebuild(c(cp_prior$bounds[1], x_range[2])))
  prior_value
}


truncate_lower_prior = function(prior_value, lower, name = NULL) {
  # Leave non-string priors (e.g. numeric constants) unchanged
  if (!is.character(prior_value))
    return(prior_value)

  # Parse distribution call and isolate any user-supplied T() clause
  parts = split_prior_truncation(prior_value)
  call = parse_prior_call(parts$distribution)

  # Only stochastic distributions without existing truncation need bounding.
  # dinvgamma() is positive and written as a transformed dgamma(), which cannot take T().
  if (is.null(call) || !grepl("^d[A-Za-z]", call$name) || !is.null(parts$truncation) || call$name == "dinvgamma")
    return(prior_value)

  # dunif cannot take T() in JAGS; bound through its lower parameter instead
  if (call$name == "dunif" && length(call$args) == 2) {
    lower_val = suppressWarnings(as.numeric(call$args[1]))
    upper_val = suppressWarnings(as.numeric(call$args[2]))
    if (!is.na(upper_val) && upper_val <= lower)
      stop("Prior", if (!is.null(name)) paste0(" for '", name, "'"), " must allow values above ", lower, ".")
    if (!is.na(lower_val) && lower_val < lower)
      return(paste0("dunif(", lower, ", ", call$args[2], ")"))
    return(prior_value)
  }

  # Truncate unbounded distributions at the lower bound
  paste0(prior_value, " T(", lower, ", )")
}


# Check that priors are set on the original parameter, not on names that same() replaced.
# - predictors, group_effects: occurrence tables, i.e., including rows reused by same().
assert_prior_not_shared = function(prior, predictors, group_effects) {
  shared = dplyr::bind_rows(
    dplyr::transmute(
      dplyr::filter(predictors, .data$segment != .data$definition_segment),
      replaced = stringr::str_replace(.data$name, "_[0-9]+$", paste0("_", .data$segment)),
      original = .data$name
    ),
    dplyr::transmute(
      dplyr::filter(group_effects, .data$segment != .data$definition_segment),
      replaced = stringr::str_replace(.data$sd_name, paste0("_", .data$definition_segment, "_"), paste0("_", .data$segment, "_")),
      original = .data$sd_name
    )
  ) %>% dplyr::distinct()
  bad_prior = shared[shared$replaced %in% names(prior), , drop = FALSE]
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


overlay_user_prior_specs = function(specs, prior, cps, context, predictors, family) {
  name_matches = names(prior) %in% specs$parameter
  if (any(!name_matches)) {
    stop(
      "Prior(s) were specified for parameter name(s) that are not part of the model: ",
      and_collapse(names(prior)[!name_matches])
    )
  }

  # dirichlet(alpha) applies to all change points, so propagate it to unassigned change points
  user_cp_names = intersect(names(prior), cps$name)
  user_cp_is_dirichlet = grepl("^\\s*dirichlet\\s*\\(", as.character(prior[user_cp_names]))
  if (any(user_cp_is_dirichlet)) {
    if (!all(user_cp_is_dirichlet))
      stop("All or none of the change point priors must be `dirichlet(alpha)`.")
    source_desc = and_collapse(paste0("`", user_cp_names, "`"))
    for (u_name in setdiff(cps$name, user_cp_names)) {
      i = match(u_name, specs$parameter)
      specs$code[i] = prior[[user_cp_names[1]]]
      specs$source[i] = "user"
      specs$description[i] = paste0("Inherited from user-specified ", source_desc)
    }
  }

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
      original = prior[[group_name]]
      is_bounded = is.character(original) && stringr::str_detect(
        original, "^\\s*dunif\\s*\\(|T\\s*\\("
      )
      if (is_distribution(original) && !is_bounded) {
        x_display = context$x_display
        prior[[group_name]] = paste0(
          original, " ", group_cp_truncation(cps, j, paste0("min(", x_display, ")"), paste0("max(", x_display, ")"))
        )
        auto_truncated = c(auto_truncated, group_name)
      }
    }
  }

  # Enforce the support of constant dpars, e.g., sigma > 0
  bounds = get_constant_dpar_bounds(family)
  for (name in intersect(names(prior), predictors$name[predictors$dpar %in% names(bounds)])) {
    dpar = predictors$dpar[match(name, predictors$name)]
    original = prior[[name]]
    prior[[name]] = truncate_lower_prior(original, bounds[[dpar]], name)
    if (!identical(prior[[name]], original))
      auto_truncated = c(auto_truncated, name)
  }

  for (name in names(prior)) {
    i = match(name, specs$parameter)
    specs$code[i] = prior[[name]]
    specs$source[i] = "user"
    specs$description[i] = if (name %in% auto_truncated && name %in% cps$name) {
      paste0("User-specified prior, truncated to the observed range of ", context$x_display)
    } else if (name %in% auto_truncated) {
      "User-specified prior with required bounds added by mcp"
    } else {
      NA_character_
    }
  }

  # Change points with a prior distribution are ordered (see get_jags_code())
  if (nrow(cps) > 1) {
    is_ordered = specs$parameter %in% cps$name & vapply(specs$code, is_distribution, logical(1))
    description = dplyr::coalesce(specs$description[is_ordered], "User-specified prior")
    specs$description[is_ordered] = paste0(description, ", restricted to ordered change points")
    if (length(user_cp_names) == 0)
      specs$description[is_ordered] = "Flat Dirichlet on the segment lengths"
  }
  specs
}


# Segment-1 x-terms that change the coefficient `name` along x. These have the same
# model-matrix column apart from x, e.g., `stateb.x` for `stateb` in `1 + state * x`.
# Intercepts (incl. each level in `0 + state`) also have x-terms without other variables, e.g., `x`.
get_x_partners = function(predictors, name) {
  i = match(name, predictors$name)
  own = predictors$matrix_data[[i]]
  is_partner = vapply(predictors$matrix_data, function(values) {
    isTRUE(all.equal(values, own)) || (predictors$prior_type[i] == "Intercept" && all(values == 1, na.rm = TRUE))
  }, logical(1))
  predictors[predictors$segment == 1 & predictors$dpar == predictors$dpar[i] & predictors$x_factor != "1" & is_partner, ]
}


# Place default segment-1 priors on the value at the segment start, min(x), like
# brms and rstanarm place intercept priors on centered predictors. The reported
# coefficients remain at x = 0 (as in lm()), so their priors are conditional on the
# coefficients' x-terms, e.g., "dnorm(mean(.y) - x_1 * min(.x), 2.5 * sd(.y))" or
# "dnorm(0 - stateb.x_1 * min(.x), ...)". Other predictors are at 0.
add_intercept_references = function(specs, predictors, context) {
  if (context$x_min == 0)
    return(specs)

  coefs = predictors[predictors$segment == 1 & predictors$x_factor == "1" & predictors$dpar %notin% c("ar", "ma"), ]
  for (i in seq_len(nrow(coefs))) {
    row = match(coefs$name[i], specs$parameter)
    if (is.na(row) || specs$source[row] != "default")
      next
    parts = split_prior_truncation(specs$code[row])
    call = parse_prior_call(parts$distribution)
    if (is.null(call) || call$name %notin% c("dnorm", "dt", "dunif"))
      next

    terms = get_x_partners(predictors, coefs$name[i])
    if (nrow(terms) == 0)
      next

    # Subtract the x-terms at min(x) from the location and the truncation bounds
    power = ifelse(terms$x_factor == "x", "", sub("x", "", terms$x_factor, fixed = TRUE))
    reference = paste0(" - ", terms$name, " * min(.x)", power, collapse = "")
    shift = function(arg) if (arg == "") "" else if (arg == "0") reference else paste0(arg, reference)
    location = if (call$name == "dunif") 1:2 else 1
    call$args[location] = vapply(call$args[location], shift, character(1))
    specs$code[row] = paste0(call$name, "(", paste(call$args, collapse = ", "), ")")
    if (!is.null(parts$truncation)) {
      bounds = vapply(parse_prior_call(parts$truncation)$args, shift, character(1))
      specs$code[row] = paste0(specs$code[row], " T(", paste(bounds, collapse = ", "), ")")
    }
    what = if (coefs$prior_type[i] == "Intercept") "value" else "difference"
    specs$description[row] = paste0(specs$description[row], "; on the ", what, " at min(", context$x_display, ")")
  }
  specs
}
