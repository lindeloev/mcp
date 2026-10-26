# ABOUT: Building tables describing model segments, change points, and parameters
# ------------------------------------------------------------



# Build the model tables for `mcp()`
get_model_tables = function(model, data, family, par_x, series = NULL) {
  # Pass the same canonical segment formulas to all model-building consumers
  parsed = parse_segments(model, data, family, par_x, series)
  coef_tables = get_coef_tables(parsed$model, data, family, par_x)
  coef_definitions = definition_rows(coef_tables$coefs)
  group_definitions = definition_rows(coef_tables$group_coefs)
  family = resolve_dpar_specs(family, coef_tables$components)
  parameters = get_pars_table(coef_definitions, parsed$cps, group_definitions, family)
  if (nrow(parameters) == 0)
    stop("The model does not contain any parameters to estimate.", call. = FALSE)

  # Fitted-model methods evaluate the designs from .GlobalEnv, so formula functions must be reachable from there
  formula_functions = unique(unlist(lapply(coef_tables$design_specs, function(design_spec) called_functions(attr(design_spec$terms, "predvars")))))
  unreachable = formula_functions[!vapply(formula_functions, exists, logical(1), envir = globalenv(), mode = "function")]
  if (length(unreachable) > 0)
    stop(
      "Function(s) used in the formulas are not available from the global environment: ", and_collapse(paste0("`", unreachable, "`")),
      ". Define them in the global environment, or add the transformed column(s) to `data`.",
      call. = FALSE
    )

  if (any(coef_definitions$dpar %in% c("ar", "ma"))) {
    check_arma_model(family, coef_definitions, parsed$data_columns, data)
    warn_arma_not_included(coef_tables$components)
  }

  list(
    model = parsed$model,
    family = family,
    model_tables = list(
      data_columns = parsed$data_columns,
      cps = parsed$cps,
      coefs = coef_tables$coefs,
      group_coefs = coef_tables$group_coefs,
      parameters = parameters,
      design_specs = coef_tables$design_specs
    )
  )
}


#' Parse a list of segment formulas
#'
#' Used internally for most mcp functions.
#'
#' @inheritParams mcp
#' @param series The `series` column from `ar()`/`ma()` or `NULL`.
#' @return A list with `model` (the canonical segment formulas), `cps` (one row per
#'   estimated change point, naturally indexed by `cp = 1, 2, ...`), and
#'   `data_columns` (the data-column roles returned by `mcp_columns()`).
#' @noRd
parse_segments = function(model, data, family, par_x, series = NULL) {
  checkmate::assert_string(par_x)

  ############################################
  # CANONICALIZE ISOLATED SEGMENT FORMULAS   #
  ############################################
  default_response = NULL
  cp_group_cols = character()  # One per change point. NA for population-level change points.
  for (i in seq_along(model)) {
    env = environment(model[[i]])
    canonical = canonicalize_segment(model[[i]], i, default_response)
    if (i == 1) {
      default_response = canonical$response
      response_columns = unpack_y(canonical$response, i, family, env)
    } else {
      cp_group_cols[i - 1] = unpack_cp(canonical$cp, i, env)
    }
    model[[i]] = canonical$form
  }

  # Every segment shares one response and one change-point dimension, including multiple-regression models.
  data_columns = c(list(par_x = par_x, response = response_columns$response, series = series), response_columns[-1])



  ###########################
  # CHECK SEGMENTS AND DATA #
  ###########################

  # Group-level effects
  if (length(unique(stats::na.omit(cp_group_cols))) > 1) {
    stop(
      "All group-level change points must use the same grouping factor. Found ",
      and_collapse(paste0("'", unique(stats::na.omit(cp_group_cols)), "'")), "."
    )
  }

  # Check y
  y_col = data_columns$response
  if (any(is.na(data[[y_col]]))) {
    if (!is.null(data_columns$weights) && isTRUE(family$response$is_discrete))
      stop("Weights with missing responses (NA in '", y_col, "') is not supported for family = ", family$family, "().")
    message("NA values detected in '", y_col, "'. JAGS will treat them as latent responses and impute them during sampling.")
  }

  assert_response_data(family, data_columns, data)


  ###############################################
  # CHANGE POINTS                               #
  ###############################################
  # `cps`: one row per *estimated* change point (one fewer than segments), naturally
  # indexed by `cp` instead of borrowing the segment index.
  # `segment` is an explicit join key: the segment that this change point starts.
  cps = tibble::tibble(
    cp = seq_along(cp_group_cols),
    segment = .data$cp + 1L,
    name = paste0("cp_", .data$cp),
    group_col = cp_group_cols,
    sd_name = dplyr::if_else(!is.na(.data$group_col), paste0(.data$name, "_", .data$group_col, "_sd"), NA_character_),
    group_name = dplyr::if_else(!is.na(.data$group_col), paste0(.data$name, "_", .data$group_col), NA_character_)
  )

  list(model = model, cps = cps, data_columns = data_columns)
}


#' Build the canonical table of model parameter definitions
#'
#' Provides the canonical display order used by `summary()`, `fixef()`,
#' `ranef()`, `prior_summary()`, and `plot_pars(pars = "population")`: change
#' points first (including their SD/group-level hyperparameters), then `mu`,
#' then the other distributional parameters in the order defined by the
#' family, then `ar`/`ma` components, each ascending by segment.
#'
#' @param coefs The `coefs` table from `get_coef_tables()`.
#' @param cps A table of change points from `parse_segments()`.
#' @param group_coefs The `group_coefs` table from `get_coef_tables()`.
#' @param family An `mcpfamily` object.
#' @return A tibble with one row per model parameter (population and group),
#'   with columns `parameter`, `part`, `scope`, `role`, `segment`, `dpar`,
#'   `group_col`, and `population_parameter`.
#' @noRd
get_pars_table = function(coefs, cps, group_coefs, family) {
  group_cps = cps[!is.na(cps$group_col), , drop = FALSE]

  # One tibble per kind of parameter. Zero-row tables recycle to zero rows. `.tie` orders each
  # change point and coefficient (by model-matrix column) before its group-level SD and deviations.
  pars = dplyr::bind_rows(
    tibble::tibble(parameter = cps$name, part = "cp", scope = "population", role = "change_point", segment = cps$segment,
      dpar = "cp", group_col = NA_character_, population_parameter = NA_character_, .tie = 0),
    tibble::tibble(parameter = group_cps$sd_name, part = "cp", scope = "population", role = "group_sd", segment = group_cps$segment,
      dpar = "cp", group_col = group_cps$group_col, population_parameter = group_cps$name, .tie = 1),
    tibble::tibble(parameter = group_cps$group_name, part = "cp", scope = "group", role = "change_point", segment = group_cps$segment,
      dpar = "cp", group_col = group_cps$group_col, population_parameter = group_cps$name, .tie = 2),
    tibble::tibble(parameter = coefs$name, part = "predictor", scope = "population", role = "coefficient", segment = coefs$segment,
      dpar = coefs$dpar, order = coefs$order, .tie = coefs$matrix_col),
    tibble::tibble(parameter = group_coefs$sd_name, part = "predictor", scope = "population", role = "group_sd", segment = group_coefs$segment,
      dpar = group_coefs$dpar, order = group_coefs$order, group_col = group_coefs$group_col,
      population_parameter = group_coefs$population_parameter, .tie = group_coefs$matrix_col + 0.1),
    tibble::tibble(parameter = group_coefs$name, part = "predictor", scope = "group", role = "coefficient", segment = group_coefs$segment,
      dpar = group_coefs$dpar, order = group_coefs$order, group_col = group_coefs$group_col,
      population_parameter = group_coefs$population_parameter, .tie = group_coefs$matrix_col + 0.2)
  )

  # Canonical group order: cp, mu, other family dpars (family order), then
  # ar/ma components (ar before ma) and their lag order.
  dpar_levels = c("cp", "mu", setdiff(family$dpar_specs$dpar, "mu"), "ar", "ma")

  # Public `dpar` values are those accepted by `dpar` arguments, e.g., "ar1". Change points have none.
  pars |>
    dplyr::arrange(factor(.data$dpar, levels = dpar_levels), .data$order, .data$segment, .data$.tie) |>
    dplyr::mutate(dpar = dplyr::case_when(
      .data$part == "cp" ~ NA_character_,
      .data$dpar %in% c("ar", "ma") ~ paste0(.data$dpar, .data$order),
      .default = .data$dpar
    )) |>
    dplyr::select("parameter", "part", "scope", "role", "segment", "dpar", "group_col", "population_parameter")
}


#' Model parameters
#'
#' `r lifecycle::badge("experimental")`
#'
#' Return the canonical parameter definitions for an `mcpfit` object. This
#' works before sampling, so you can discover the model parameters before setting priors or inits.
#'
#' @param fit An `mcpfit` object.
#' @param scope Optional parameter scope(s): `"population"` or `"group"`.
#' @param role Optional parameter role(s): `"change_point"`, `"coefficient"`
#'   (regression coefficients of all distributional parameters, e.g., `mu`, `sigma`, and `ar1`),
#'   or `"group_sd"`. Combine with `scope` to select population-level parameters or their
#'   group-level deviations, e.g., `scope = "group", role = "change_point"` for group-level change points.
#' @return A tibble with one row per parameter definition. `parameter` is the
#'   model parameter name; `dpar` is its distributional parameter (e.g., `"mu"`,
#'   `"sigma"`, or `"ar1"`; see [mcp-formula]), as accepted by the `dpar` argument of
#'   [fitted.mcpfit()], [fixef.mcpfit()], etc., and `NA` for change points;
#'   `group_col` and `population_parameter` describe group-level effects.
#' @export
#' @examples
#' # Show every parameter in the model
#' print(mcp_pars(demo_fit), width = Inf)
#'
#' # Select population-level coefficients
#' mcp_pars(demo_fit, scope = "population", role = "coefficient")$parameter
mcp_pars = function(fit, scope = NULL, role = NULL) {
  checkmate::assert_class(fit, "mcpfit")
  if (!is.null(scope))
    scope = rlang::arg_match(scope, c("population", "group"), multiple = TRUE)
  if (!is.null(role))
    checkmate::assert_character(role, any.missing = FALSE)

  parameters = get_fit_model_tables(fit)$parameters
  keep = rep(TRUE, nrow(parameters))
  if (!is.null(scope))
    keep = keep & parameters$scope %in% scope
  if (!is.null(role))
    keep = keep & parameters$role %in% role
  parameters[keep, , drop = FALSE]
}


#' Model data columns
#'
#' `r lifecycle::badge("experimental")`
#'
#' Return the resolved data-column roles for an `mcpfit` object. In particular,
#' `par_x` is the change-point variable chosen automatically or supplied to
#' [mcp()]. Family-defined response auxiliaries, such as `trials` and
#' `weights`, are listed when relevant.
#'
#' @param fit An `mcpfit` object.
#' @return A named list with `par_x`, `response`, and `series`, plus any
#'   family-defined response-auxiliary column roles.
#' @export
#' @examples
#' # Show the predictor, response, and auxiliary-column roles
#' mcp_columns(demo_fit)
mcp_columns = function(fit) {
  checkmate::assert_class(fit, "mcpfit")
  get_fit_model_tables(fit)$data_columns
}


# Get column names that serve as predictors in the model
get_predictor_cols = function(fit) {
  model_tables = get_fit_model_tables(fit)
  design_vars = unlist(lapply(model_tables$design_specs, function(design_spec) {
    if (!is.null(design_spec$terms)) all.vars(design_spec$terms) else character()
  }))
  unique(c(mcp_columns(fit)$par_x, setdiff(design_vars, ".mcp_local_x")))
}


# Rows that define parameters (and priors). Coefficients in same() also
# occur in later segments, where `segment != definition_segment`.
definition_rows = function(table) {
  table[which(table$segment == table$definition_segment), , drop = FALSE]
}


# Names of population-level change points starting a disjoined segment, i.e., a segment that
# includes an intercept for every dpar with x-terms in it, and which no x-terms from the previous
# segment join, e.g., `sigma(1 + x)` in segment 1 applies as `min(x, cp_1)` in a segment 2 without sigma().
# Excludes change points with group-level change points.
get_disjoined_cps = function(fit) {
  model_tables = get_fit_model_tables(fit)
  terms = dplyr::bind_rows(model_tables$coefs, model_tables$group_coefs)
  cps = model_tables$cps[is.na(model_tables$cps$group_col), ]
  is_disjoined = vapply(cps$segment, function(segment) {
    segment_terms = terms[terms$segment == segment, ]
    joined_x = terms$segment == segment - 1 & terms$x_degree > 0 & dplyr::coalesce(terms$next_segment > segment, TRUE)
    !any(joined_x) && all(segment_terms$dpar[segment_terms$x_degree > 0] %in% segment_terms$dpar[segment_terms$par_type == "Intercept"])
  }, logical(1))
  cps$name[is_disjoined]
}


# Get model metadata tables from a fitted model
get_fit_model_tables = function(fit) {
  check_mcpfit_version(fit)
  fit$.internal$model_tables
}
