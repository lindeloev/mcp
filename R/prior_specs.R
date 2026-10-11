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
      code = "dirichlet(1)",
      description = "Uniform order statistics (flat Dirichlet) within the observed change-point span",
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
      description = "Ordered group-level change-point deviations from the population location",
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
  # For x-terms, a factor dummy (as in x:state) only selects rows, so the column sd is close to sd(x).
  data_scale = if (n_unique <= 1 || (n_unique == 2 && x_factor != "1")) 1 else stats::sd(values)

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


# Retrieve the offset spec included for a specific dpar and segment, if any.
get_segment_offset = function(design_specs, dpar, segment) {
  specs = Filter(function(s) isTRUE(s$has_offset) && identical(s$dpar, dpar) && s$segment == segment, design_specs)
  if (length(specs) == 0) return(NULL)
  spec = specs[[1]]
  if (is.null(spec$offset_data) || all(spec$offset_data == 0)) NULL else spec
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


default_predictor_specs = function(predictors, family, design_specs = list()) {
  defaults = filter_default_conditions(dplyr::bind_rows(family$default_prior, default_arma_specs()), family)

  keys = paste(defaults$dpar, defaults$par_type)
  if (anyDuplicated(keys))
    stop_github("Default prior specifications are not unique by dpar and par_type.")

  joined = predictors %>%
    dplyr::left_join(defaults, by = c("dpar", "par_type"))
  if (any(is.na(joined$prior))) {
    stop_github(
      "mcp could not find a default prior for ",
      and_collapse(joined$code_name[is.na(joined$prior)])
    )
  }

  # Offset adjustments for log-link mu intercepts
  if (identical(family$link, "log")) {
    offset_specs = Filter(function(s) isTRUE(s$has_offset) && identical(s$dpar, "mu") && any(s$offset_data != 0), design_specs)
    for (i in which(joined$dpar == "mu")) {
      spec = get_segment_offset(design_specs, "mu", joined$segment[i])
      if (!is.null(spec)) {
        sym = if (length(offset_specs) <= 1) "offset" else paste0("offset_", spec$segment)
        if (joined$par_type[i] == "Intercept")
          joined$prior[i] = gsub("log(pmax(.y, 0.1))", paste0("log(pmax(.y, 0.1)) - ", sym), joined$prior[i], fixed = TRUE)
        joined$description[i] = gsub("log-count|log-mean", "log-rate", joined$description[i])
      }
    }
  }

  joined$prior = fill_predictor_scale(joined$prior, joined$matrix_data, joined$x_factor)

  tibble::tibble(
    parameter = joined$code_name,
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
    paste0("coefficient `", joined$display_name, "`")
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


truncate_cp_prior = function(cps, j, prior_value, context) {
  if (is.numeric(prior_value))
    return(prior_value)
  is_bounded = stringr::str_detect(prior_value, "^\\s*(dunif|dirichlet)\\s*\\(")
  is_truncated = stringr::str_detect(prior_value, "T\\s*\\(")
  if (is_bounded || is_truncated)
    return(prior_value)

  lower = if (j == 1) {
    paste0("min(", context$x_display, ")")
  } else {
    cps$name[j - 1]
  }
  paste0(prior_value, " T(", lower, ", max(", context$x_display, "))")
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
      replaced = stringr::str_replace(.data$code_name, "_[0-9]+$", paste0("_", .data$segment)),
      original = .data$code_name
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

  user_cp_names = intersect(names(prior), cps$name)
  if (length(user_cp_names) > 0) {
    user_cp_is_dirichlet = grepl("^\\s*dirichlet\\s*\\(", as.character(prior[user_cp_names]))

    if (any(!user_cp_is_dirichlet)) {
      if (any(user_cp_is_dirichlet)) {
        stop("All or none of the change point priors must be `dirichlet(alpha)`.")
      }
      # When user specifies non-dirichlet cp prior(s), any unassigned default
      # change points fall back to direct priors.
      unassigned_cps = setdiff(cps$name, user_cp_names)
      for (u_name in unassigned_cps) {
        j = match(u_name, cps$name)
        lower = if (j == 1) "min(.x)" else cps$name[j - 1]
        i = match(u_name, specs$parameter)
        specs$code[i] = paste0("dunif(", lower, ", max(.x))")
        specs$description[i] = if (j == 1) {
          "Within the observed change-point span"
        } else {
          paste0("Ordered after ", cps$name[j - 1], " within the observed change-point span")
        }
      }
    } else if (all(user_cp_is_dirichlet) && length(user_cp_names) < nrow(cps)) {
      # Propagate user-specified alpha to unassigned change points
      first_dirichlet = prior[[user_cp_names[1]]]
      source_desc = and_collapse(paste0("`", user_cp_names, "`"))
      unassigned_cps = setdiff(cps$name, user_cp_names)
      for (u_name in unassigned_cps) {
        i = match(u_name, specs$parameter)
        specs$code[i] = first_dirichlet
        specs$source[i] = "user"
        specs$description[i] = paste0("Inherited from user-specified ", source_desc)
      }
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
      if (is.character(original) && !is_bounded) {
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
  for (name in intersect(names(prior), predictors$code_name[predictors$dpar %in% names(bounds)])) {
    dpar = predictors$dpar[match(name, predictors$code_name)]
    original = prior[[name]]
    prior[[name]] = truncate_lower_prior(original, bounds[[dpar]], name)
    if (!identical(prior[[name]], original))
      auto_truncated = c(auto_truncated, name)
  }

  for (name in names(prior)) {
    i = match(name, specs$parameter)
    specs$code[i] = prior[[name]]
    specs$source[i] = "user"
    specs$description[i] = if (name %in% auto_truncated) {
      "User-specified prior with required bounds added by mcp"
    } else {
      NA_character_
    }
  }
  specs
}


# Place default segment-1 intercept priors on the level at the segment start,
# min(x), like brms and rstanarm place them on centered predictors. The reported
# intercept remains at x = 0 (as in lm()). `reference` is the value of the
# segment's x-terms at min(x), e.g., "x_1 * 101 + xE2_1 * 10201". Other
# predictors are at 0. Returns `table` with a `reference` column.
add_intercept_references = function(table, predictors, context) {
  table$reference = NA_character_
  if (context$x_min == 0)
    return(table)

  intercepts = predictors[predictors$segment == 1 & predictors$par_type == "Intercept" & predictors$dpar %notin% c("ar", "ma"), ]
  for (i in seq_len(nrow(intercepts))) {
    row = match(intercepts$code_name[i], table$parameter)
    if (is.na(row) || table$source[row] != "default" || table$kind[row] != "distribution")
      next

    # x-terms of the same dpar without other variables (a model-matrix column of ones)
    terms = predictors[predictors$segment == 1 & predictors$dpar == intercepts$dpar[i] & predictors$x_factor != "1", ]
    is_pure = vapply(terms$matrix_data, function(values) all(values == 1, na.rm = TRUE), logical(1))
    terms = terms[is_pure, ]
    if (nrow(terms) == 0)
      next

    degree = as.numeric(ifelse(terms$x_factor == "x", "1", sub("x^", "", terms$x_factor, fixed = TRUE)))
    table$reference[row] = paste0(terms$code_name, " * ", sprintf("%.15g", context$x_min^degree), collapse = " + ")
    table$description[row] = paste0(table$description[row], "; on the level at min(", context$x_display, ")")
  }
  table
}
