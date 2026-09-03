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
    "ar", "slope", "dnorm(0, 0.25 / predictor_scale())", "Modest dependence-coefficient change over a representative predictor change", "always",
    "ma", "Intercept", "dnorm(0, 0.5) T(-1, 1)", "Zero-centered regularizing dependence coefficient", "always",
    "ma", "dummy", "dnorm(0, 0.25)", "Modest categorical change in a dependence coefficient", "always",
    "ma", "slope", "dnorm(0, 0.25 / predictor_scale())", "Modest dependence-coefficient change over a representative predictor change", "always"
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
      description = "Zero-mean ordered group-level change-point offsets",
      source = "default"
    )
  }

  dplyr::bind_rows(specs)
}


# Express a reference change in one model-matrix column for prior scaling.
default_predictor_scale = function(matrix_data, x_factor) {
  values = stats::na.omit(as.numeric(matrix_data))
  unique_values = unique(values)
  data_scale = if (length(unique_values) <= 1) {
    1
  } else if (length(unique_values) == 2) {
    diff(range(unique_values))
  } else {
    2 * stats::sd(values)
  }
  if (!is.finite(data_scale) || data_scale <= 0)
    stop_github("Could not derive a positive scale from a model-matrix column.")

  parts = character()
  if (x_factor != "1") {
    x_expr = if (grepl("\\^", x_factor)) "(max(.x) - min(.x))" else "max(.x) - min(.x)"
    parts = gsub("x", x_expr, x_factor, fixed = TRUE)
    if (data_scale != 1 && !grepl("^\\(", parts))
      parts = paste0("(", parts, ")")
  }
  if (data_scale != 1)
    parts = c(parts, format_prior_number(data_scale))
  if (length(parts) == 0)
    parts = "1"
  paste0("(", paste(parts, collapse = " * "), ")")
}


default_predictor_specs = function(predictors, family) {
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


overlay_user_prior_specs = function(specs, prior, cps, context) {
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
      unassigned_cps = setdiff(cps$name, user_cp_names)
      for (u_name in unassigned_cps) {
        i = match(u_name, specs$parameter)
        specs$code[i] = first_dirichlet
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
  }

  for (name in names(prior)) {
    i = match(name, specs$parameter)
    specs$code[i] = prior[[name]]
    specs$source[i] = "user"
    specs$description[i] = if (name %in% auto_truncated) {
      "User-specified prior with ordered change-point bounds added by mcp"
    } else {
      NA_character_
    }
  }
  specs
}
