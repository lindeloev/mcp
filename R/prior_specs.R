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
    "ar", "Intercept", "dunif(-1, 1)", "Bounded dependence coefficient", "always",
    "ar", "dummy", "dt(0, 1, 3)", "Current weak prior for a categorical dependence coefficient", "always",
    "ar", "slope", "dt(0, 1 / (max(.x) - min(.x)), 3)", "One coefficient-unit change across the observed change-point span", "always",
    "ma", "Intercept", "dunif(-1, 1)", "Bounded dependence coefficient", "always",
    "ma", "dummy", "dt(0, 1, 3)", "Current weak prior for a categorical dependence coefficient", "always",
    "ma", "slope", "dt(0, 1 / (max(.x) - min(.x)), 3)", "One coefficient-unit change across the observed change-point span", "always"
  )
}


default_cp_specs = function(CP, context) {
  n_cp = context$n_cp
  if (n_cp == 0)
    return(empty_prior_specs())

  specs = list()
  for (j in seq_len(nrow(CP))) {
    name = CP$name[j]
    if (n_cp == 1) {
      code = "dunif(min(.x), max(.x))"
    } else {
      lower = if (j == 1) "min(.x)" else CP$name[j - 1]
      code = paste0(
        "dt(min(.x), (max(.x) - min(.x)) / n_cp(), n_cp() - 1) T(",
        lower, ", max(.x))"
      )
    }
    specs[[name]] = tibble::tibble(
      parameter = name,
      code = code,
      description = if (j == 1) {
        "Within the observed change-point span"
      } else {
        paste0("Ordered after ", CP$name[j - 1], " within the observed change-point span")
      },
      source = "default"
    )
  }

  for (j in seq_len(nrow(CP))) {
    if (!CP$varying[j])
      next

    sd_name = CP$sd_name[j]
    group_name = CP$group_name[j]
    specs[[sd_name]] = tibble::tibble(
      parameter = sd_name,
      code = "dnorm(0, 2 * (max(.x) - min(.x)) / n_cp()) T(0, )",
      description = "Group-level change-point variation",
      source = "default"
    )

    lower = if (j == 1) {
      paste0("min(.x) - ", CP$name[j])
    } else {
      paste0(CP$name[j - 1], " - ", CP$name[j])
    }
    upper = if (j == nrow(CP)) {
      paste0("max(.x) - ", CP$name[j])
    } else {
      paste0(CP$name[j + 1], " - ", CP$name[j])
    }
    specs[[group_name]] = tibble::tibble(
      parameter = group_name,
      code = paste0("dnorm(0, ", sd_name, ") T(", lower, ", ", upper, ")"),
      description = "Zero-centered group-level change-point offsets",
      source = "default"
    )
  }

  dplyr::bind_rows(specs)
}


default_rhs_specs = function(rhs_table, family) {
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

  joined = rhs_table %>%
    dplyr::left_join(defaults, by = c("dpar", "par_type"))
  if (any(is.na(joined$prior))) {
    stop_github(
      "mcp could not find a default prior for ",
      and_collapse(joined$code_name[is.na(joined$prior)])
    )
  }

  tibble::tibble(
    parameter = joined$code_name,
    code = joined$prior,
    description = joined$description,
    source = "default"
  )
}


truncate_cp_prior = function(CP, j, prior_value, context) {
  if (is.numeric(prior_value))
    return(prior_value)
  is_bounded = stringr::str_detect(prior_value, "^\\s*(dunif|dirichlet)\\s*\\(")
  is_truncated = stringr::str_detect(prior_value, "T\\s*\\(")
  if (is_bounded || is_truncated)
    return(prior_value)

  lower = if (j == 1) {
    paste0("min(", context$x_display, ")")
  } else {
    CP$name[j - 1]
  }
  paste0(prior_value, " T(", lower, ", max(", context$x_display, "))")
}


overlay_user_prior_specs = function(specs, prior, CP, context) {
  name_matches = names(prior) %in% specs$parameter
  if (any(!name_matches)) {
    stop(
      "Prior(s) were specified for parameter name(s) that are not part of the model: ",
      and_collapse(names(prior)[!name_matches])
    )
  }

  auto_truncated = character()
  for (j in seq_len(nrow(CP))) {
    name = CP$name[j]
    if (name %in% names(prior)) {
      original = prior[[name]]
      prior[[name]] = truncate_cp_prior(CP, j, original, context)
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
