# ABOUT: These functions build the change point regression model for R and JAGS
# -----------------


#' Call `get_formula_jags_dpar` for multiple dpars and paste strings
#'
#' @aliases get_formula_jags
#' @keywords internal
#' @noRd
#' @inheritParams get_formula_jags_dpar
#' @inheritParams mcp
#' @param dpars A character vector of dpars to including in model building
#' @return A string with JAGS code.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_formula_jags = function(segments, predictors, group_effects, par_x, family) {
  # Explicit segment -> boundary-code lookup instead of relying on row
  # order matching segment numbers.
  boundary_code = stats::setNames(segments$cp_code_form, segments$segment)

  # Add X-helpers which code the X relative to the start of each segment.
  local_x_str = "\n# par_x local to each segment"
  for (i in seq_len(nrow(segments))) {
    segment_start = ifelse(i > 1, yes = paste0(" - ", boundary_code[[as.character(i)]]), no = "")  #
    segment_end = ifelse(i < nrow(segments), yes = boundary_code[[as.character(i + 1)]], no = paste0("cp_", i))  # infinite if last segment.

    local_x_str = paste0(local_x_str, "\nx_local_", i, "_[i_] = min(", par_x, "[i_], ", segment_end, ")", segment_start)
  }

  # Build formula for each dpar (note plural "_dpars")
  this_cp_lookup = dplyr::select(segments, "segment", "form", this_cp = "cp_code_form")
  next_cp_lookup = dplyr::select(segments, next_intercept = "segment", next_cp = "cp_code_form")

  predictor_group_effects = group_effects %>%
    dplyr::filter(.data$part == "predictor") %>%
    dplyr::transmute(
      dpar = .data$dpar,
      segment = .data$segment,
      matrix_col = .data$matrix_col,
      code_name = paste0(.data$name, "[", .data$group_col, "[i_]]"),
      order = .data$order,
      x_factor = .data$x_factor,
      next_intercept = .data$next_segment
    )
  formula_predictors = predictors %>%
    dplyr::select(
      "dpar", "segment", "matrix_col", "code_name", "order", "x_factor",
      "next_intercept"
    ) %>%
    dplyr::bind_rows(predictor_group_effects)

  formula_jags_dpars = formula_predictors %>%
    dplyr::left_join(this_cp_lookup, by = "segment") %>%
    dplyr::left_join(next_cp_lookup, by = "next_intercept") %>%
    dplyr::mutate(
      # One dpar per ar order: (ar, 1) --> ar1
      dpar = paste0(.data$dpar, tidyr::replace_na(as.character(.data$order), ""))
    ) %>%

    # Build formula for each dpar
    dplyr::group_by(.data$dpar) %>%
    tidyr::nest() %>%
    dplyr::rowwise() %>%
    dplyr::summarise(
      formula_jags_dpar = get_formula_jags_dpar(.data$data, .data$dpar, par_x, family)
    ) %>%
    dplyr::pull(.data$formula_jags_dpar) %>%
    paste0(collapse = "\n\n")

  garma_boundary_str = get_garma_boundary_jagscode(segments, predictors, par_x)

  # Concatenate and return
  formula_jags = paste0(local_x_str, garma_boundary_str, "\n\n", formula_jags_dpars)

  # Special case when no terms are present for a given dpar (all ~0): insert "dpar = 0".
  for (dpar in family$dpar_specs$dpar) {
    if (dpar %notin% formula_predictors$dpar)
      formula_jags = paste0(formula_jags, "\n\n# All segments are ~ 0 for this par:\nlink_", dpar, "_[i_] = 0")
  }

  # Return with nicer printing
  class(formula_jags) = c("mcptext", "character")
  formula_jags
}


#' Build an R formula (as string) for a dpar
#'
#' @aliases get_formula_jags_dpar
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param dpar_table Rows of `predictors` for one `(dpar, order)` pair.
#' @param family An mcpfamily object with distributional-parameter metadata.
#' @return A string with JAGS code.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_formula_jags_dpar = function(dpar_table, dpar, par_x, family) {
  # Build this! Initiate
  formula_str = paste0("# Formula for ", dpar)
  is_distributional = dpar %in% family$dpar_specs$dpar
  node_name = ifelse(is_distributional, paste0("link_", dpar), dpar)
  formula_str = paste0(formula_str, "\n", node_name, "_[i_] =\n")


  df_code_strs = dpar_table %>%
    # Code dpar_data_segment column indices
    dplyr::arrange(.data$segment) %>%  # just to make sure
    dplyr::group_by(.data$segment) %>%

    # Summarise
    dplyr::group_by(.data$segment, .data$x_factor, .data$next_cp) %>%
    dplyr::summarise(
      # The parts
      indicator_this = paste0("  (", par_x, "[i_] >= ", dplyr::first(.data$this_cp), ")"),
      indicator_next = dplyr::if_else(is.na(dplyr::first(.data$next_cp)) == TRUE, "", paste0(" * (", par_x, "[i_] < ", dplyr::first(.data$next_cp), ")")),
      inprod = paste0(" * inprod(rhs_matrix_[i_, c(", paste0(.data$matrix_col, collapse = ", "), ")], c(", paste0(.data$code_name, collapse = ", "), "))"),
      x_factor = gsub("x(?!p\\()", paste0("x_local_", dplyr::first(.data$segment), "_[i_]"), dplyr::first(.data$x_factor), perl = TRUE),  # "x" but not "exp("

      # All together
      segment_code = paste0(.data$indicator_this, .data$indicator_next, .data$inprod, " * ", .data$x_factor),
      form = dplyr::first(.data$form)
    ) %>%

    # Add title-comment
    dplyr::group_by(.data$segment) %>%
    dplyr::mutate(
      title = dplyr::if_else(dplyr::row_number() == 1, paste0("\n  # Segment ", dplyr::first(.data$segment), ": ", dplyr::first(.data$form), "\n"), ""),
      segment_code = paste0(.data$title, .data$segment_code)
    )

  # Return
  all_predictors = paste0(df_code_strs$segment_code, collapse = " + \n")
  formula_str = paste0(formula_str, all_predictors)

  formula_str
}


#' Convert for-looped JAGS code to vectorized R code
#'
#' @aliases get_formula_r
#' @keywords internal
#' @noRd
#' @param formula_jags Character, often residing in `fit$.internal$formula_jags`.
#' @param predictors Output of `get_predictors()`.
#' @param group_effects Output of `get_group_effects()`.
#' @param cps Output of `get_segment_tables()`.
#' @param par_x Name of the change-point predictor.
#' @return Character
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_formula_r = function(formula_jags, predictors, group_effects, cps, par_x) {
  sim_pars = get_sim_pars(cps, predictors, group_effects)
  predictor_pars = predictors$code_name
  predictor_group_effects = group_effects[group_effects$part == "predictor", , drop = FALSE]
  group_pars = predictor_group_effects$name
  cp_pars = setdiff(sim_pars, c(predictor_pars, group_pars))

  # Replacements that turns rowwise JAGS code into vectorized R code
  replace_args = c(
    # Group-level predictor coefficients
    if (length(group_pars) > 0) stats::setNames(
        paste0("args$", group_pars),
        paste0(group_pars, "[", predictor_group_effects$group_col, "]")
      ),

    # Predictor part
    stats::setNames(paste0(", args$", sim_pars), paste0(", ", sim_pars)),
    if (length(predictor_pars) > 0) stats::setNames(
        paste0("args$", predictor_pars, ", "),
        paste0(predictor_pars, ", ")
      ),
    stats::setNames(paste0("cbind(args$"), "cbind("),
    stats::setNames("args$", "args$args$"),  # Fix double-inserting args$ above

    # Change points
    stats::setNames(paste0("args$", par_x, " >="), paste0(par_x, " >=")),
    stats::setNames(paste0("args$", par_x, " <"), paste0(par_x, " <")),

    # General
    stats::setNames("pmin(args$", paste0("pmin("))
  )
  if (length(cp_pars) > 0) {
    replace_args = c(
      replace_args,
      stats::setNames(paste0(" (args$", cp_pars, " + "), paste0(" (", cp_pars, " + ")),  # group-level change point
      stats::setNames(paste0(" + args$", cp_pars, ")"), paste0(" + ", cp_pars, ")")),  # group-level change point
      stats::setNames(paste0(">= args$", cp_pars), paste0(">= ", cp_pars)),
      stats::setNames(paste0("< args$", cp_pars), paste0("< ", cp_pars)),
      stats::setNames(paste0(") - args$", cp_pars), paste0(") - ", cp_pars))
    )
  }

  formula_r = formula_jags %>%
    stringr::str_remove_all("\\[i_\\]") %>%  # Vectorized
    stringi::stri_replace_all_fixed("[i_,", "[,") %>%  # Vectorized
    stringi::stri_replace_all_fixed("min(", "pmin(") %>%  # Vectorized
    stringi::stri_replace_all_fixed("max(", "pmax(") %>%  # Vectorized
    stringi::stri_replace_all_fixed(")],", "), drop = FALSE],") %>%  # Prevent reducing matrix to vector for one-column indexing
    stringi::stri_replace_all_fixed("], c(", "], cbind(") %>%  # Vectorized
    stringr::str_remove_all("CP_[0-9]+_INDEX") %>%  # Only used for JAGS code; not in R.
    stringi::stri_replace_all_fixed(names(replace_args), replace_args, vectorize_all = FALSE)  # obs: fixed to not interpret $ as regex

  class(formula_r) = c("mcptext", "character")  # Nicer printing
  formula_r
}
