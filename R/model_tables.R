# ABOUT: Building tables describing model segments, change points, and parameters
# ------------------------------------------------------------



#' Build tables describing a list of segments
#'
#' Used internally for most mcp functions.
#'
#' @aliases get_segment_tables
#' @keywords internal
#' @inheritParams mcp
#' @return A list with `segments` (one row per segment, including change-point
#'   code for the change point starting that segment) and `cps` (one row per
#'   estimated change point, naturally indexed by `cp = 1, 2, ...`).
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @noRd
get_segment_tables = function(model, data = NULL, family = gaussian(), par_x) {
  checkmate::assert_string(par_x)

  ############################################
  # BUILD SEGMENTS FROM ISOLATED FORMULAS    #
  ############################################
  segments = tibble::tibble()
  for (i in seq_along(model)) {
    row = tibble::tibble(segment = i)
    row = dplyr::bind_cols(row, unpack_tildes(model[[i]], i))
    row = dplyr::bind_cols(row, unpack_y(row$form_y, i, family))
    row = dplyr::bind_cols(row, unpack_cp(row$form_cp, i))

    segments = dplyr::bind_rows(segments, row)
  }

  # Fill y and trials, where not explicit.
  # Build "full" formula (with explicit intercepts) and insert instead of the old
  segments = segments %>%
    tidyr::fill("y", "form_y", "trials", "weights", .direction = "downup") %>%  # Usually only provided in segment 1
    dplyr::mutate(form = ifelse(.data$segment == 1, .data$form, paste0(.data$form_y, .data$form_cp, " ~ ", .data$form_rhs))) %>%  # build full formula
    dplyr::select(-"form_y", -"form_cp", -"form_rhs")  # Not needed anymore

  # Every segment shares one change-point dimension, including multiple-regression models.
  segments$x = par_x



  ###########################
  # CHECK SEGMENTS AND DATA #
  ###########################

  # Check segment 1: change point not possible here
  if (any(segments[1, c("cp_int", "cp_varying", "cp_group_col")] != FALSE, na.rm = TRUE))
    stop("Change point defined in first segment. This should not be possible. Submit bug report in the GitHub repo.")

  # Response variables
  derived_y = unique(stats::na.omit(segments$y))
  if (length(derived_y) != 1)
    stop("There should be exactly one response variable. Found ", and_collapse(derived_y), " across segments.")

  aux_columns = get_family_aux_columns(family, segments)

  # Group-level effects
  derived_varying = unique(stats::na.omit(segments$cp_group_col))

  # Check data types
  if (!is.null(data)) {
    # Check y
    if (!is.numeric(data[, segments$y[1]]))
      stop("Data column '", segments$y[1], "' has to be numeric.")
    if (any(is.na(data[, segments$y[1]])))
      message("NA values detected in '", segments$y[1], "'. JAGS will treat them as latent responses and impute them during sampling.")

    # Check varying
    if (length(derived_varying) > 0) {
      for (varying_col in derived_varying) {
        data_varying = data[, varying_col]
        if (!is.character(data_varying) && !is.factor(data_varying))
          if (!all(data_varying == floor(data_varying)))
            stop("Varying group '", varying_col, "' has to be integer, character, or factor.")
      }
    }

    response_data = get_family_response_data(family, segments, data)
    response_columns = c(list(y = segments$y[1]), as.list(aux_columns))
    family$response$validate(data[[segments$y[1]]], response_data, response_columns)
  }


  ###############################################
  # SPLIT INTO `segments` AND `cps`              #
  ###############################################
  # `segments`: one row per segment, shared metadata only.
  segment_metadata = segments %>%
    dplyr::select("segment", "form", "y", "trials", "weights", "x")

  # `cps`: one row per *estimated* change point (nrow(segments) - 1 of
  # them), naturally indexed by `cp` instead of borrowing the segment index.
  # `segment` is an explicit join key: the segment that this change point starts.
  cps = segments %>%
    dplyr::filter(.data$segment > 1) %>%
    dplyr::transmute(
      cp = .data$segment - 1,
      segment = .data$segment,
      name = paste0("cp_", .data$cp),
      varying = .data$cp_varying,
      group_col = .data$cp_group_col,
      sd_name = ifelse(.data$varying, paste0(.data$name, "_sd"), NA),
      group_name = ifelse(.data$varying, paste0(.data$name, "_", .data$group_col), NA)
    ) %>%
    dplyr::mutate(
      code = ifelse(!is.na(.data$group_name), yes = paste0(.data$name, " + ", .data$group_name, "CP_", .data$segment, "_INDEX"), no = .data$name),
      code = format_code(.data$code, na_col = .data$name)
    )

  cp_group_cols = unique(stats::na.omit(cps$group_col))
  if (length(cp_group_cols) > 1) {
    stop(
      "All group-level change points must use the same grouping factor. Found ",
      and_collapse(paste0("'", cp_group_cols, "'")), "."
    )
  }

  # `segments`: shared segment metadata with change-point info joined back in
  # by `segment`. Segment 1 has no row in `cps` (there is no change point before
  # it), so it is given the fixed lower boundary `cp_0` explicitly. Kept in
  # this per-segment shape because most consumers look up the change point
  # *starting* a given segment, not the change point by its own number.
  segments = segment_metadata %>%
    dplyr::left_join(
      dplyr::select(
        cps, "segment",
        cp_name = "name", cp_group_col = "group_col",
        cp_sd = "sd_name", cp_group = "group_name", cp_code_form = "code"
      ),
      by = "segment"
    ) %>%
    dplyr::mutate(
      cp_name = dplyr::coalesce(.data$cp_name, "cp_0"),
      cp_code_form = dplyr::coalesce(.data$cp_code_form, "cp_0")
    )

  list(segments = segments, cps = cps)
}


#' Build a table of parameter names with their segment and dpar
#'
#' Provides the canonical display order used by `summary()`, `fixef()`,
#' `ranef()`, `prior_summary()`, and `plot_pars(pars = "population")`: change
#' points first (including their SD/group-level hyperparameters), then `mu`,
#' then the other distributional parameters in the order declared by the
#' family, then `ar`/`ma` components (combined with their lag, e.g. `"ar1"`),
#' each ascending by segment.
#'
#' @aliases get_pars_table
#' @keywords internal
#' @noRd
#' @param predictors A table from `get_predictors()`.
#' @param cps A table of change points from `get_segment_tables()`.
#' @param group_effects A table from `get_group_effects()`.
#' @param family An `mcpfamily` object.
#' @return A tibble with one row per model parameter (population and
#'   group), with columns `name`, `part`, `scope`, `role`, `segment`, and
#'   `dpar`.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_pars_table = function(predictors, cps, group_effects, family) {
  cp_pars = tibble::tibble(
    name = character(), part = character(), scope = character(), role = character(),
    segment = integer(), dpar = character(), .tie = integer()
  )
  if (nrow(cps) > 0) {
    cp_pars = tibble::tibble(
      name = cps$name, part = "cp", scope = "population", role = "change_point",
      segment = cps$segment, dpar = "cp", .tie = 0L
    )
    varying_cp = cps[cps$varying, , drop = FALSE]
    if (nrow(varying_cp) > 0) {
      cp_pars = dplyr::bind_rows(
        cp_pars,
        tibble::tibble(
          name = varying_cp$sd_name, part = "cp", scope = "population", role = "group_sd",
          segment = varying_cp$segment, dpar = "cp", .tie = 1L
        ),
        tibble::tibble(
          name = varying_cp$group_name, part = "cp", scope = "group", role = "group_deviation",
          segment = varying_cp$segment, dpar = "cp", .tie = 2L
        )
      )
    }
  }

  predictor_pars = predictors %>%
    dplyr::transmute(
      name = .data$code_name,
      part = "predictor",
      scope = "population",
      role = as.character(ifelse(.data$dpar %in% family$dpar_specs$dpar, "fixed_effect", "arma")),
      segment = .data$segment,
      dpar = as.character(ifelse(.data$dpar %in% c("ar", "ma"), paste0(.data$dpar, .data$order), .data$dpar)),
      .tie = .data$matrix_col
    )

  predictor_group_pars = group_effects %>%
    dplyr::filter(.data$part == "predictor") %>%
    dplyr::transmute(
      name = .data$sd_name,
      part = "predictor",
      scope = "population",
      role = "group_sd",
      segment = .data$segment,
      dpar = .data$dpar,
      .tie = .data$matrix_col + 0.1
    ) %>%
    dplyr::bind_rows(
      group_effects %>%
        dplyr::filter(.data$part == "predictor") %>%
        dplyr::transmute(
          name = .data$name,
          part = "predictor",
          scope = "group",
          role = "group_deviation",
          segment = .data$segment,
          dpar = .data$dpar,
          .tie = .data$matrix_col + 0.2
        )
    )

  # Canonical group order: cp, mu, other family dpars (declared order), then
  # ar/ma labels sorted by component ("ar" before "ma") and then lag order.
  arma_labels = unique(predictor_pars$dpar[predictor_pars$dpar %notin% c("mu", family$dpar_specs$dpar)])
  arma_labels = arma_labels[order(
    match(sub("[0-9]+$", "", arma_labels), c("ar", "ma")),
    as.integer(sub("^[a-z]+", "", arma_labels))
  )]
  dpar_levels = c("cp", "mu", setdiff(family$dpar_specs$dpar, "mu"), arma_labels)

  dplyr::bind_rows(cp_pars, predictor_pars, predictor_group_pars) %>%
    dplyr::mutate(dpar = factor(.data$dpar, levels = dpar_levels)) %>%
    dplyr::arrange(.data$dpar, .data$segment, .data$.tie) %>%
    dplyr::mutate(dpar = as.character(.data$dpar)) %>%
    dplyr::select("name", "part", "scope", "role", "segment", "dpar")
}


#' Build the table of group-level effects
#'
#' @keywords internal
#' @noRd
#' @param cps A table from `get_segment_tables()`.
#' @param predictor_group_effects Predictor-side rows returned by
#'   `get_predictor_tables()`.
#' @return A tibble with one row per group-level effect.
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
      display_name = character(),
      order = integer(),
      x_factor = character(),
      design_id = character(),
      design_col = integer(),
      matrix_col = integer(),
      matrix_data = list(),
      next_segment = integer(),
      correlated = logical()
    )
  }

  cp_group_effects = cps %>%
    dplyr::filter(.data$varying) %>%
    dplyr::transmute(
      population_name = as.character(.data$name),
      name = as.character(.data$group_name),
      part = "cp",
      group_col = as.character(.data$group_col),
      segment = .data$segment,
      dpar = "cp",
      sd_name = as.character(.data$sd_name)
    )

  dplyr::bind_rows(cp_group_effects, predictor_group_effects)
}


#' Get model metadata tables from a fitted model
#'
#' Fits store consistently named tables in `.internal$model_tables`.
#'
#' @keywords internal
#' @noRd
#' @param fit An `mcpfit` object.
#' @return A list with `segments`, `cps`, `predictors`, `group_effects`, `pars`,
#'   and fitted `design_specs`.
get_fit_model_tables = function(fit) {
  check_mcpfit_version(fit)
  fit$.internal$model_tables
}



#' Format code with one or multiple terms
#'
#' Take a value like "a + b" and
#' (1) replace it with NA if na_col == NA.
#' (2) Change to "(a + b)" if there is a "+"
#' (3) Return itself otherwise, e.g., "a" --> "a".
#'
#' @aliases format_code
#' @keywords internal
#' @noRd
#' @param col A column
#' @param na_col If this column is NA, return NA
#' @return string
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
format_code = function(col, na_col) {
  # Replace with NA
  col = ifelse(is.na(na_col), NA, col)

  # Add parenthesis if this is an equation
  col = ifelse(stringr::str_detect(col, "\\+"), paste0("(", col, ")"), col)
  col
}
