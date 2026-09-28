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
    if (!cps$varying[j])
      next

    sd_name = cps$sd_name[j]
    group_name = cps$group_name[j]
    specs[[sd_name]] = tibble::tibble(
      parameter = sd_name,
      code = "dnorm(0, 2 * (max(.x) - min(.x)) / n_cp()) T(0, )",
      description = "Group-level change-point variation",
      source = "default"
    )

    previous_cp = if (j == 1) "min(.x)" else cps$code[j - 1]
    previous_cp = stringr::str_replace(
      previous_cp, "CP_[0-9]+_INDEX", paste0("[", cps$group_col[j], "_]")
    )
    later_population = which(seq_len(nrow(cps)) > j & !cps$varying)
    next_cp = if (length(later_population) == 0) "max(.x)" else
      cps$name[later_population[1]]
    lower = paste0(previous_cp, " - ", cps$name[j])
    upper = paste0(next_cp, " - ", cps$name[j])
    specs[[group_name]] = tibble::tibble(
      parameter = group_name,
      code = paste0("dnorm(0, ", sd_name, ") T(", lower, ", ", upper, ")"),
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


default_predictor_specs = function(predictors, family, design_specs = list()) {
  defaults = dplyr::bind_rows(family$default_prior, default_arma_specs())

  modeled_dpars = unique(defaults$dpar[defaults$condition == "modeled"])
  for (dpar in modeled_dpars) {
    is_modeled = get_dpar_spec(family, dpar)$modeled
    condition = if (is_modeled) "modeled" else "constant"
    defaults = defaults %>%
      dplyr::filter(
        .data$dpar != .env$dpar |
          .data$condition == "always" |
          .data$condition == .env$condition
      )
  }

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

  scaled = grepl("predictor_scale()", joined$prior, fixed = TRUE)
  joined$prior[scaled] = vapply(which(scaled), function(i) {
    gsub(
      "predictor_scale()",
      default_predictor_scale(joined$matrix_data[[i]], joined$x_factor[i]),
      joined$prior[i],
      fixed = TRUE
    )
  }, character(1))

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

  defaults = family$default_prior
  modeled_dpars = unique(defaults$dpar[defaults$condition == "modeled"])
  for (dpar in modeled_dpars) {
    is_modeled = get_dpar_spec(family, dpar)$modeled
    condition = if (is_modeled) "modeled" else "constant"
    defaults = defaults %>%
      dplyr::filter(
        .data$dpar != .env$dpar |
          .data$condition == "always" |
          .data$condition == .env$condition
      )
  }

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

  scaled = grepl("predictor_scale()", joined$group_sd_prior, fixed = TRUE)
  joined$group_sd_prior[scaled] = vapply(which(scaled), function(i) {
    gsub(
      "predictor_scale()",
      default_predictor_scale(joined$matrix_data[[i]], joined$x_factor[i]),
      joined$group_sd_prior[i],
      fixed = TRUE
    )
  }, character(1))
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


truncate_sigma_prior = function(prior_value, lower_floor, name = NULL) {
  # Leave non-string priors (e.g. numeric constants) unchanged
  if (!is.character(prior_value))
    return(prior_value)

  # Parse distribution call and isolate any user-supplied T() clause
  parts = split_prior_truncation(prior_value)
  call = parse_prior_call(parts$distribution)

  # Only stochastic distributions without existing truncation need bounding
  if (is.null(call) || !grepl("^d[A-Za-z]", call$name) || !is.null(parts$truncation))
    return(prior_value)

  # dunif cannot take T() in JAGS; bound through its lower parameter instead
  if (call$name == "dunif" && length(call$args) == 2) {
    lower_val = suppressWarnings(as.numeric(call$args[1]))
    upper_val = suppressWarnings(as.numeric(call$args[2]))
    if (!is.na(upper_val) && upper_val <= lower_floor)
      stop("Prior", if (!is.null(name)) paste0(" for '", name, "'"), " must allow values of at least ", format_prior_number(lower_floor), ".")
    if (!is.na(lower_val) && lower_val < lower_floor)
      return(paste0("dunif(", format_prior_number(lower_floor), ", ", call$args[2], ")"))
    return(prior_value)
  }

  # Truncate unbounded distributions at the family floor
  paste0(prior_value, " T(", format_prior_number(lower_floor), ", )")
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
    if (cps$varying[j] && group_name %in% names(prior)) {
      original = prior[[group_name]]
      is_bounded = is.character(original) && stringr::str_detect(
        original, "^\\s*dunif\\s*\\(|T\\s*\\("
      )
      if (is.character(original) && !is_bounded) {
        previous_cp = if (j == 1) paste0("min(", context$x_display, ")") else cps$code[j - 1]
        previous_cp = stringr::str_replace(
          previous_cp, "CP_[0-9]+_INDEX", paste0("[", cps$group_col[j], "_]"))
        later_population = which(seq_len(nrow(cps)) > j & !cps$varying)
        next_cp = if (length(later_population) == 0) {
          paste0("max(", context$x_display, ")")
        } else {
          cps$name[later_population[1]]
        }
        prior[[group_name]] = paste0(
          original, " T(", previous_cp, " - ", cps$name[j], ", ",
          next_cp, " - ", cps$name[j], ")"
        )
        auto_truncated = c(auto_truncated, group_name)
      }
    }
  }

  # Enforce positive floor on unmodeled residual standard deviation priors
  sigma_pars = predictors$code_name[predictors$dpar == "sigma"]
  if (length(sigma_pars) > 0 && !get_dpar_spec(family, "sigma")$modeled) {
    floor = get_dpar_spec(family, "sigma")$lower
    for (name in intersect(names(prior), sigma_pars)) {
      original = prior[[name]]
      prior[[name]] = truncate_sigma_prior(original, floor, name)
      if (!identical(prior[[name]], original))
        auto_truncated = c(auto_truncated, name)
    }
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


# Improve efficiency of Gibbs-like sampling by keeping each segment's line relatively
# independent/constant when a change point moves. This is a reparameterization - not
# a change in model. Otherwise, the sampler must change slopes and intercepts jointly 
# with change points, which mixes slowly. Two strategies are used:
# 
# - Joined segments (no intercept): sample each default slope as its rise over the
#   segment, rise ~ normal(0, s * width^p), and derive slope = rise / width^p.
# 
# - Disjoined segments k >= 2: sample the default intercept as the level at the segment
#   end, end ~ normal(m + change, s), and derive Intercept_k = end - change, where
#   change is the segment's x-terms over its width.
# 
# Both imply exactly the default standard priors; only the sampling geometry changes.
# width = cp_k - cp_(k-1). Skipped for segments bounded by group-varying change points,
# where the data use group-specific widths. Returns `table` with columns `rise_width`
# and `end_change` (NA where not applicable).
add_segment_anchors = function(table, predictors, cps) {
  table$rise_width = NA_character_
  table$end_change = NA_character_
  is_pure = vapply(predictors$matrix_data, function(values) all(values == 1, na.rm = TRUE), logical(1))
  slopes = predictors[predictors$dpar == "mu" & predictors$x_factor != "1" & is_pure & predictors$segment == predictors$definition_segment, ]
  is_default_normal = function(row) {
    call = parse_prior_call(table$value[row])
    !is.na(row) && table$source[row] == "default" && !is.null(call) && call$name == "dnorm" && !grepl("T(", table$value[row], fixed = TRUE)
  }
  width = function(k, x_factor) paste0("(cp_", k, " - cp_", k - 1, ")", ifelse(x_factor == "x", "", sub("x", "", x_factor, fixed = TRUE)))

  # Apply to each segment
  for (k in unique(slopes$segment)) {
    # Don't change anything for segments bounded by group-varying change points
    if (any(cps$varying[cps$name %in% paste0("cp_", c(k - 1, k))]))
      next

    # Joined segments (no intercept): sample each default slope as its rise over the segment
    terms = slopes[slopes$segment == k, ]
    intercept = predictors$code_name[predictors$segment == k & predictors$dpar == "mu" & predictors$par_type == "Intercept"]
    if (length(intercept) == 0) {
      for (i in seq_len(nrow(terms))) {
        row = match(terms$code_name[i], table$parameter)
        if (is_default_normal(row))
          table$rise_width[row] = width(k, terms$x_factor[i])
      }
    
    # Disjoined segments: sample the default intercept as the level at the segment end
    } else if (k > 1) {
      row = match(intercept, table$parameter)
      if (is_default_normal(row))
        table$end_change[row] = paste0(terms$code_name, " * ", width(k, terms$x_factor), collapse = " + ")
    }
  }
  table
}
