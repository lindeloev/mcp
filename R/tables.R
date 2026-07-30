#' Takes a formula and returns a string representation of y, cp, and rhs
#' @aliases unpack_tildes
#' @keywords internal
#' @noRd
#' @param form A formula
#' @param i The segment number
#' @return A one-row tibble with columns:
#'   * `form`: String. The full formula for this segment.
#'   * `form_y`: String. The expression for y (without tilde)
#'   * `form_cp`: String. The formula for the change point.
#'   * `form_rhs`: String. The predictor formula. Only used to build the
#'     formula representation in `summary.mcpfit()`.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_tildes = function(form, i) {
  has_LHS = attributes(stats::terms(form))$response == 1
  form_str = formula_to_char(form)
  if (has_LHS == FALSE && i == 1) {
    stop("No response variable in segment 1.")
  } else if (has_LHS == FALSE && i > 1) {
    # If no LHS, add a change point "intercept"
    form_str = paste("1", form_str)
  }

  # List of strings for each section
  chunks = stringr::str_trim(strsplit(form_str, "~")[[1]])

  if (length(chunks) == 2) {
    # Only one tilde. This is the first segment or y is implicit from earlier segment(s)
    return(tibble::tibble(
      form = form_str,
      form_y = ifelse(i == 1, chunks[1], NA),
      form_cp = ifelse(i == 1, NA, paste0(" ~ ", chunks[1])),
      form_rhs = chunks[2]
    ))
  } else if (length(chunks) == 3) {
    if (i == 1)
      stop("The first segment must have exactly one tilde. Got two.")

    return(tibble::tibble(
      form = form_str,
      form_y = chunks[1],
      form_cp = paste0(" ~ ", chunks[2]),
      form_rhs = chunks[3]
    ))
  } else {
    stop("Error in segment ", i, ": Got none or more than two ~ in a segment formula.")
  }
}


#' Unpacks y variable name
#'
#' @aliases unpack_y
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param form_y Character representation of formula
#' @param i Segment number
#' @return A one-row tibble with the response and auxiliary-data columns.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_y = function(form_y, i, family) {
  declared = names(family$response$auxiliary)
  aux_names = unique(c("trials", "weights", declared))
  response = stats::setNames(as.list(rep(NA_character_, length(aux_names) + 1)), c("y", aux_names))

  # If NA and not segment 1, just return empty
  if (is.na(form_y)) {
    if (i == 1)
      stop("A response must be defined in segment 1, e.g., 'y ~ 1'")

    return(tibble::as_tibble(response))
  }


  # Split by |
  y_split = strsplit(form_y, "\\|")[[1]]
  if (length(y_split) > 2)
    stop("There can only be zero or one pipe in response. Got '", form_y, "' in segment ", i)

  # RESPONSE
  lhs = y_split[1]
  y_col = attr(stats::terms(to_formula(lhs)), "term.labels")
  if (length(y_col) != 1)
    stop("There should be exactly one response variable. Got ", length(y_col), " in segment ", i)
  response$y = y_col

  term_labels = character()
  if (length(y_split) == 2) {
    rhs = y_split[2]
    term_labels = attr(stats::terms(to_formula(rhs)), "term.labels")
    ok_terms = if (length(declared) == 0) rep(FALSE, length(term_labels)) else
      vapply(term_labels, function(term) any(stringr::str_detect(term, paste0("^", declared, "\\("))), logical(1))
    if (!all(ok_terms))
      stop(
        "Only ", if (length(declared) == 0) "no terms are" else and_collapse(paste0("`", declared, "()`")),
        " allowed after the pipe for family = ", family$family, "(). Got '", rhs, "'."
      )
  }

  for (name in declared) {
    term_index = stringr::str_detect(term_labels, paste0("^", name, "\\("))
    got_term = any(term_index)
    if (family$response$auxiliary[[name]]$required && !got_term)
      stop("Error in response of segment ", i, ": need a valid ", name, "() specification.")
    if (!got_term)
      next
    if (sum(term_index) > 1)
      stop("Only one ", name, "() term is allowed in segment ", i, ".")

    term = term_labels[term_index]
    content = get_term_content(term)
    column = attr(stats::terms(content), "term.labels")
    if (length(column) != 1)
      stop("There must be exactly one term inside ", name, "(). Got ", term, " in segment ", i)
    response[[name]] = column
  }

  tibble::as_tibble(response)
}


#' Takes a cp formula (as a string) and returns its properties
#'
#' @aliases unpack_cp
#' @keywords internal
#' @noRd
#' @param form_cp Segment formula as string.
#' @param i segment number
#' @return A one-row tibble with columns:
#'   * `cp_int`: bool. Whether there is an intercept change in the change point.
#'   * `cp_varying`: bool or NA. Is there a group-level intercept on the change point?
#'   * `cp_group_col`: char or NA. Which data column defines the grouping factor?
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
unpack_cp = function(form_cp, i) {
  if (is.na(form_cp)) {
    return(tibble::tibble(
      cp_int = FALSE,
      cp_varying = FALSE,
      cp_group_col = NA
    ))
  }
  form_cp = stats::as.formula(form_cp)

  # Group-level effects
  form_varying = remove_terms(form_cp, "population")

  if (!is.null(form_varying)) {
    varying_terms = attr(stats::terms(form_varying), "term.labels")
    if (length(varying_terms) > 1)
      stop("Error in segment", i, " (change point): only one varying effect allowed. Found ", form_cp)

    varying_parts = strsplit(gsub(" ", "", varying_terms), "\\|")[[1]]
    if (!varying_parts[1] == "1")
      stop("Error in segment ", i, " (change point): Only plain intercepts are allowed in varying effects, e.g., (1|id).", i)

    if (!grepl("^[A-Za-z._0-9]+$", varying_parts[2]))
      stop("Error in segment ", i, " (change point): invalid format of grouping variable in varying effects. Got: ", varying_parts[2])
  }

  # Population-level effects
  attrs = attributes(stats::terms(remove_terms(form_cp, "varying")))
  if (length(attrs$term.labels) > 0)
    stop("Error in segment ", i, " (change point): Only intercepts (1) are allowed in population-level effects.")

  if (is.null(form_varying) && attrs$intercept == 0)
    stop("Error in segment ", i, " (change point): no intercept or varying effect. You can do e.g., ~ 1 or ~ (1 |id).")

  # Return as list.
  if (!is.null(form_varying)) {
    # If there is a varying effect
    return(tibble::tibble(
      cp_int = attrs$intercept == 1,
      cp_varying = ifelse(varying_parts[1] == "1", TRUE, NA),  # placeholder for later
      cp_group_col = varying_parts[2]
    ))
  } else {
    # If there is no varying effect
    return(tibble::tibble(
      cp_int = attrs$intercept == 1,
      cp_varying = FALSE,
      cp_group_col = NA
    ))
  }
}



# #' Unpack varying effects
# #'
# #' @aliases unpack_varying_term
# #' @keywords internal
# #' @param term
# #' @return A "one-row" list describing a varying intercept.
# #' @encoding UTF-8
# #' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
# unpack_varying_term = function(term, i) {
#   parts = stringr::str_trim(strsplit(term, "\\|")[[1]])  # as c("lhs", "group")
#
#   # Check that there is just one grouping term
#   if (!grepl("^[A-Za-z._0-9]+$", parts[2]))
#     stop("Error in segment ", i, " (right-hand side): Grouping variable in varying effects for change points.")
#
#   # LHS: Split intercepts and variable
#   preds = strsplit(gsub(" ", "", parts[1]), "\\+")[[1]]
#   slope = preds[(preds %in% c("0", "1")) == FALSE]
#   if (length(slope) > 1)
#     stop("Error in segment ", i, " (right-hand side): More than one variable in LHS of varying effect.")
#   else if (length(slope) == 0)
#     # If not slope
#     slope = NA
#
#   # Return
#   list(
#     int = "0" %notin% preds,  # bool. Is intercept present?
#     slope = slope,
#     group = parts[2]
#   )
# }



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
#'   group), with columns `name`, `part`, `scope`, `segment`, and `dpar`.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_pars_table = function(predictors, cps, group_effects, family) {
  cp_pars = tibble::tibble(
    name = character(), part = character(), scope = character(),
    segment = integer(), dpar = character(), .tie = integer()
  )
  if (nrow(cps) > 0) {
    cp_pars = tibble::tibble(
      name = cps$name, part = "cp", scope = "population",
      segment = cps$segment, dpar = "cp", .tie = 0L
    )
    varying_cp = cps[cps$varying, , drop = FALSE]
    if (nrow(varying_cp) > 0) {
      cp_pars = dplyr::bind_rows(
        cp_pars,
        tibble::tibble(
          name = varying_cp$sd_name, part = "cp", scope = "population",
          segment = varying_cp$segment, dpar = "cp", .tie = 1L
        ),
        tibble::tibble(
          name = varying_cp$group_name, part = "cp", scope = "group",
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
    dplyr::select("name", "part", "scope", "segment", "dpar")
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
#' New fits store consistently named tables in `.internal$model_tables`. The
#' fallback keeps methods working for legacy saved fits.
#'
#' @keywords internal
#' @noRd
#' @param fit An `mcpfit` object.
#' @return A list with `segments`, `cps`, `predictors`, `group_effects`, and
#'   `pars`.
get_fit_model_tables = function(fit) {
  if (!is.null(fit$.internal$model_tables))
    return(fit$.internal$model_tables)

  segments = fit$.internal$ST
  if (is.null(segments) && !is.null(fit$.other))
    segments = fit$.other$ST
  cps = fit$.internal$CP
  if (is.null(cps) && !is.null(segments)) {
    cps = segments %>%
      dplyr::filter(.data$segment > 1) %>%
      dplyr::transmute(
        cp = .data$segment - 1,
        segment = .data$segment,
        name = .data$cp_name,
        varying = !is.na(.data$cp_group),
        group_col = .data$cp_group_col,
        sd_name = .data$cp_sd,
        group_name = .data$cp_group,
        code = .data$cp_code_form
      )
  }

  pars = fit$.internal$pars_table
  if (!is.null(pars) && "part" %notin% names(pars)) {
    group_names = stats::na.omit(cps$group_name)
    pars$part = ifelse(pars$dpar == "cp", "cp", "predictor")
    pars$scope = ifelse(pars$name %in% group_names, "group", "population")
    pars = dplyr::select(pars, "name", "part", "scope", "segment", "dpar")
  }

  list(
    segments = segments,
    cps = cps,
    predictors = fit$.internal$rhs_table,
    group_effects = get_group_effects(cps),
    pars = pars
  )
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
