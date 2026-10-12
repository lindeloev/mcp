# ABOUT: Building tables describing model segments, change points, and parameters
# ------------------------------------------------------------



#' Build tables describing a list of segments
#'
#' Used internally for most mcp functions.
#'
#' @inheritParams mcp
#' @return A list with `segments` (one row per segment, including change-point
#'   code for the change point starting that segment) and `cps` (one row per
#'   estimated change point, naturally indexed by `cp = 1, 2, ...`).
#' @noRd
get_segment_tables = function(model, data, family, par_x) {
  checkmate::assert_string(par_x)

  ############################################
  # BUILD SEGMENTS FROM ISOLATED FORMULAS    #
  ############################################
  segments = tibble::tibble()
  default_response = NULL
  for (i in seq_along(model)) {
    env = environment(model[[i]])
    canonical = canonicalize_segment(model[[i]], i, default_response)
    if (i == 1)
      default_response = canonical$response

    row = tibble::tibble(
      segment = i,
      form = formula_to_char(canonical$form),
      form_env = list(env)
    )
    row = dplyr::bind_cols(row, unpack_y(canonical$response, i, family, env))
    row = dplyr::bind_cols(row, unpack_cp(canonical$cp, i, env))

    segments = dplyr::bind_rows(segments, row)
  }

  # Every segment shares one change-point dimension, including multiple-regression models.
  segments$x = par_x



  ###########################
  # CHECK SEGMENTS AND DATA #
  ###########################

  # Response variables
  derived_y = unique(stats::na.omit(segments$y))
  if (length(derived_y) != 1)
    stop("There should be exactly one response variable. Found ", and_collapse(derived_y), " across segments.")

  # Group-level effects
  cp_group_cols = unique(stats::na.omit(segments$cp_group_col))

  # Check y
  y_col = segments$y[1]
  if (any(is.na(data[[y_col]]))) {
    has_weights = any(!is.na(segments$weights))
    if (has_weights && isTRUE(family$response$is_discrete))
      stop("Weights with missing responses (NA in '", y_col, "') is not supported for family = ", family$family, "().")
    message("NA values detected in '", y_col, "'. JAGS will treat them as latent responses and impute them during sampling.")
  }

  # Check grouping factors
  for (group_col in cp_group_cols) {
    group_data = data[, group_col]
    if (is.numeric(group_data) && any(group_data != floor(group_data)))
      stop("Grouping factor '", group_col, "' must be integer, character, or factor.")
  }

  assert_response_data(family, segments, data)


  ###############################################
  # SPLIT INTO `segments` AND `cps`              #
  ###############################################
  # `segments`: one row per segment, shared metadata only.
  segment_metadata = segments %>%
    dplyr::select("segment", "form", "form_env", "y", "trials", "weights", "x")

  # `cps`: one row per *estimated* change point (nrow(segments) - 1 of
  # them), naturally indexed by `cp` instead of borrowing the segment index.
  # `segment` is an explicit join key: the segment that this change point starts.
  cps = segments %>%
    dplyr::filter(.data$segment > 1) %>%
    dplyr::transmute(
      cp = .data$segment - 1,
      segment = .data$segment,
      name = paste0("cp_", .data$cp),
      group_col = .data$cp_group_col,
      sd_name = ifelse(!is.na(.data$group_col), paste0(.data$name, "_", .data$group_col, "_sd"), NA),
      group_name = ifelse(!is.na(.data$group_col), paste0(.data$name, "_", .data$group_col), NA)
    ) %>%
    dplyr::mutate(
      code = ifelse(!is.na(.data$group_name), yes = paste0("(", .data$name, " + ", .data$group_name, "CP_", .data$segment, "_INDEX)"), no = .data$name)
    )

  if (length(cp_group_cols) > 1) {
    stop(
      "All group-level change points must use the same grouping factor. Found ",
      and_collapse(paste0("'", cp_group_cols, "'")), "."
    )
  }

  # `segments`: shared segment metadata with the code of the change points that start
  # and end each segment. Segment 1 starts at the fixed lower data limit `cp_0` and
  # the last segment ends at the upper data limit `cp_K`.
  segments = segment_metadata %>%
    dplyr::left_join(dplyr::select(cps, "segment", cp_group_col = "group_col", cp_start = "code"), by = "segment") %>%
    dplyr::mutate(
      cp_start = dplyr::coalesce(.data$cp_start, "cp_0"),
      cp_end = c(.data$cp_start[-1], paste0("cp_", dplyr::n()))
    )

  list(segments = segments, cps = cps)
}


#' Build the canonical table of model parameter definitions
#'
#' Provides the canonical display order used by `summary()`, `fixef()`,
#' `ranef()`, `prior_summary()`, and `plot_pars(pars = "population")`: change
#' points first (including their SD/group-level hyperparameters), then `mu`,
#' then the other distributional parameters in the order defined by the
#' family, then `ar`/`ma` components, each ascending by segment.
#'
#' @param predictors A table from `get_predictors()`.
#' @param cps A table of change points from `get_segment_tables()`.
#' @param group_effects A table from `get_group_effects()`.
#' @param family An `mcpfamily` object.
#' @return A tibble with one row per model parameter (population and group),
#'   with columns `parameter`, `part`, `scope`, `role`, `segment`, `dpar`,
#'   `group_col`, and `population_parameter`.
#' @noRd
get_pars_table = function(predictors, cps, group_effects, family) {
  cp_pars = tibble::tibble(
    parameter = character(), part = character(), scope = character(), role = character(),
    segment = integer(), dpar = character(), order = integer(),
    group_col = character(), population_parameter = character(), .tie = integer()
  )
  if (nrow(cps) > 0) {
    cp_pars = tibble::tibble(
      parameter = cps$name, part = "cp", scope = "population", role = "change_point",
      segment = cps$segment, dpar = "cp", order = NA_integer_,
      group_col = NA_character_, population_parameter = NA_character_, .tie = 0L
    )
    group_cps = cps[!is.na(cps$group_col), , drop = FALSE]
    if (nrow(group_cps) > 0) {
      cp_pars = dplyr::bind_rows(
        cp_pars,
        tibble::tibble(
          parameter = group_cps$sd_name, part = "cp", scope = "population", role = "group_sd",
          segment = group_cps$segment, dpar = "cp", order = NA_integer_,
          group_col = group_cps$group_col, population_parameter = group_cps$name, .tie = 1L
        ),
        tibble::tibble(
          parameter = group_cps$group_name, part = "cp", scope = "group", role = "change_point",
          segment = group_cps$segment, dpar = "cp", order = NA_integer_,
          group_col = group_cps$group_col, population_parameter = group_cps$name, .tie = 2L
        )
      )
    }
  }

  predictor_pars = predictors %>%
    dplyr::transmute(
      parameter = .data$name,
      part = "predictor",
      scope = "population",
      role = "coefficient",
      segment = .data$segment,
      dpar = .data$dpar,
      order = .data$order,
      group_col = NA_character_,
      population_parameter = NA_character_,
      .tie = .data$matrix_col
    )

  predictor_group_pars = group_effects %>%
    dplyr::filter(.data$part == "predictor") %>%
    dplyr::transmute(
      parameter = .data$sd_name,
      part = "predictor",
      scope = "population",
      role = "group_sd",
      segment = .data$segment,
      dpar = .data$dpar,
      order = .data$order,
      group_col = .data$group_col,
      population_parameter = .data$population_name,
      .tie = .data$matrix_col + 0.1
    ) %>%
    dplyr::bind_rows(
      group_effects %>%
        dplyr::filter(.data$part == "predictor") %>%
        dplyr::transmute(
          parameter = .data$name,
          part = "predictor",
          scope = "group",
          role = "coefficient",
          segment = .data$segment,
          dpar = .data$dpar,
          order = .data$order,
          group_col = .data$group_col,
          population_parameter = .data$population_name,
          .tie = .data$matrix_col + 0.2
        )
    )

  # Canonical group order: cp, mu, other family dpars (family order), then
  # ar/ma components (ar before ma) and their lag order.
  dpar_levels = c("cp", "mu", setdiff(family$dpar_specs$dpar, "mu"), "ar", "ma")

  # Public `dpar` values are those accepted by `dpar` arguments, e.g., "ar1". Change points have none.
  dplyr::bind_rows(cp_pars, predictor_pars, predictor_group_pars) %>%
    dplyr::arrange(factor(.data$dpar, levels = dpar_levels), .data$order, .data$segment, .data$.tie) %>%
    dplyr::mutate(dpar = dplyr::case_when(
      .data$part == "cp" ~ NA_character_,
      .data$dpar %in% c("ar", "ma") ~ paste0(.data$dpar, .data$order),
      .default = .data$dpar
    )) %>%
    dplyr::select("parameter", "part", "scope", "role", "segment", "dpar", "group_col", "population_parameter")
}


#' Model parameters
#'
#' `r lifecycle::badge("experimental")`
#'
#' Return the canonical parameter definitions for an `mcpfit` object. This
#' works before sampling and is the stable way to discover model parameters.
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
#' mcp_pars(demo_fit)
#'
#' # Select population-level coefficients
#' mcp_pars(demo_fit, scope = "population", role = "coefficient")
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
#' `par_x` is the change-point predictor chosen automatically or supplied to
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
  design_vars = unlist(lapply(model_tables$design_specs, function(spec) {
    if (!is.null(spec$terms)) all.vars(spec$terms) else character()
  }))
  unique(c(mcp_columns(fit)$par_x, setdiff(design_vars, ".mcp_local_x")))
}


#' Build the table of group-level effects
#'
#' @param cps A table from `get_segment_tables()`.
#' @param predictor_group_effects Predictor-side rows returned by
#'   `get_predictor_tables()`.
#' @return A tibble with one row per group-level effect.
#' @noRd
get_group_effects = function(cps, predictor_group_effects = NULL) {
  if (is.null(predictor_group_effects) || nrow(predictor_group_effects) == 0) {
    predictor_group_effects = tibble::tibble(
      population_name = character(),
      name = character(),
      part = character(),
      group_col = character(),
      segment = integer(),
      dpar = character(),
      sd_name = character(),
      par_type = character(),
      matrix_name = character(),
      order = integer(),
      x_factor = character(),
      design_id = character(),
      design_col = integer(),
      matrix_col = integer(),
      matrix_data = list(),
      next_segment = integer(),
      definition_segment = integer()
    )
  }

  cp_group_effects = cps %>%
    dplyr::filter(!is.na(.data$group_col)) %>%
    dplyr::transmute(
      population_name = as.character(.data$name),
      name = as.character(.data$group_name),
      part = "cp",
      group_col = as.character(.data$group_col),
      segment = .data$segment,
      dpar = "cp",
      sd_name = as.character(.data$sd_name),
      definition_segment = .data$segment
    )

  dplyr::bind_rows(cp_group_effects, predictor_group_effects)
}


# Rows that define parameters (and priors). Coefficients shared via same() also
# occur in later segments, where `segment != definition_segment`.
definition_rows = function(table) {
  table[which(table$segment == table$definition_segment), , drop = FALSE]
}


# Get model metadata tables from a fitted model
get_fit_model_tables = function(fit) {
  check_mcpfit_version(fit)
  fit$.internal$model_tables
}
