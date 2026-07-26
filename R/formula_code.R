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
get_formula_jags = function(ST, rhs_table, par_x, family) {
  # Explicit segment -> boundary-code lookup instead of relying on ST's row
  # order matching segment numbers.
  boundary_code = stats::setNames(ST$cp_code_form, ST$segment)

  # Add X-helpers which code the X relative to the start of each segment.
  local_x_str = "\n# par_x local to each segment"
  for (i in seq_len(nrow(ST))) {
    segment_start = ifelse(i > 1, yes = paste0(" - ", boundary_code[[as.character(i)]]), no = "")  #
    segment_end = ifelse(i < nrow(ST), yes = boundary_code[[as.character(i + 1)]], no = paste0("cp_", i))  # infinite if last segment.

    local_x_str = paste0(local_x_str, "\nx_local_", i, "_[i_] = min(", par_x, "[i_], ", segment_end, ")", segment_start)
  }

  # Build formula for each dpar (note plural "_dpars")
  this_cp_lookup = dplyr::select(ST, "segment", "form", this_cp = "cp_code_form")
  next_cp_lookup = dplyr::select(ST, next_intercept = "segment", next_cp = "cp_code_form")

  formula_jags_dpars = rhs_table %>%
    dplyr::select(-"matrix_data") %>%  # Throw less data around
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

  garma_boundary_str = get_garma_boundary_jagscode(ST, rhs_table, par_x)

  # Concatenate and return
  formula_jags = paste0(local_x_str, garma_boundary_str, "\n\n", formula_jags_dpars)

  # Special case when no terms are present for a given dpar (all ~0): insert "dpar = 0".
  for (dpar in family$dpar_specs$dpar) {
    if (dpar %notin% rhs_table$dpar)
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
#' @param dpar_table A rhs_table with only one (dpar, order) combo
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
    dplyr::group_by(.data$segment, .data$x_factor) %>%
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


#' Build the observation-boundary formula used by GARMA terms
#'
#' A boundary supplied with an AR or MA term remains active until the next such
#' term. The first supplied boundary also applies to earlier observations so
#' they can safely be used as lags.
#'
#' @keywords internal
#' @noRd
get_garma_boundary_jagscode = function(ST, rhs_table, par_x) {
  boundary_table = rhs_table %>%
    dplyr::filter(.data$dpar %in% c("ar", "ma"), !is.na(.data$boundary)) %>%
    dplyr::distinct(.data$segment, .data$boundary) %>%
    dplyr::arrange(.data$segment)

  if (nrow(boundary_table) == 0)
    return("")
  if (anyDuplicated(boundary_table$segment))
    stop_github("Found multiple GARMA boundaries in one segment.")

  boundary_code = stats::setNames(ST$cp_code_form, ST$segment)
  boundary_parts = character(nrow(boundary_table))
  for (i in seq_len(nrow(boundary_table))) {
    lower = if (i == 1) "" else paste0("(", par_x, "[i_] >= ", boundary_code[[as.character(boundary_table$segment[i])]], ") * ")
    upper = if (i == nrow(boundary_table)) "" else paste0("(", par_x, "[i_] < ", boundary_code[[as.character(boundary_table$segment[i + 1])]], ") * ")
    boundary_value = sprintf("%.15g", boundary_table$boundary[i])
    boundary_parts[i] = paste0("  ", lower, upper, boundary_value)
  }

  paste0(
    "\n\n# GARMA observation boundary\n",
    "garma_boundary_[i_] =\n",
    paste0(boundary_parts, collapse = " +\n")
  )
}


#' Convert for-looped JAGS code to vectorized R code
#'
#' @aliases get_formula_r
#' @keywords internal
#' @noRd
#' @param formula_jags Character, often residing in `fit$.internal$formula_jags`.
#' @param rhs_table Output of `get_rhs_table()`
#' @param pars The list that ends up in `fit$pars`
#' @return Character
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_formula_r = function(formula_jags, rhs_table, pars) {
  sim_pars = get_sim_pars(rhs_table, pars)
  rhs_pars = rhs_table$code_name
  cp_pars = setdiff(sim_pars, rhs_pars)

  # Replacements that turns rowwise JAGS code into vectorized R code
  replace_args = c(
    # RHS
    stats::setNames(paste0(", args$", sim_pars), paste0(", ", sim_pars)),
    stats::setNames(paste0("args$", rhs_pars, ", "), paste0(rhs_pars, ", ")),
    stats::setNames(paste0("cbind(args$"), "cbind("),
    stats::setNames("args$", "args$args$"),  # Fix double-inserting args$ above

    # Change points
    stats::setNames(paste0("args$", pars$x, " >="), paste0(pars$x, " >=")),
    stats::setNames(paste0("args$", pars$x, " <"), paste0(pars$x, " <")),

    # General
    stats::setNames("pmin(args$", paste0("pmin("))
  )
  if (length(cp_pars) > 0) {
    replace_args = c(
      replace_args,
      stats::setNames(paste0(" (args$", cp_pars, " + "), paste0(" (", cp_pars, " + ")),  # varying change point
      stats::setNames(paste0(" + args$", cp_pars, ")"), paste0(" + ", cp_pars, ")")),  # varying change point
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


#' Get JAGS code for GARMA residual recursion
#'
#' @aliases get_arma_jagscode get_ar_jagscode
#' @keywords internal
#' @noRd
#' @param ar_order,ma_order Positive integer or `NA` when absent.
#' @param x_name Character. Name of some vector that has the length of the dataset.
#' @return Character JAGS code
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_arma_jagscode = function(ar_order, ma_order, x_name) {
  ar_order = ifelse(is.na(ar_order), 0, ar_order)
  ma_order = ifelse(is.na(ma_order), 0, ma_order)
  checkmate::assert_int(ar_order, lower = 0)
  checkmate::assert_int(ma_order, lower = 0)
  checkmate::assert_string(x_name)
  max_order = max(ar_order, ma_order)
  if (max_order == 0)
    stop_github("get_arma_jagscode() requires a positive AR or MA order.")

  get_terms = function(row, available_ar, available_ma) {
    terms = character()
    if (available_ar > 0) {
      lags = seq_len(available_ar)
      terms = c(terms, paste0("ar", lags, "_[", row, "] * resid_abs_[", row, " - ", lags, "]"))
    }
    if (available_ma > 0) {
      lags = seq_len(available_ma)
      terms = c(terms, paste0("ma", lags, "_[", row, "] * resid_ma_[", row, " - ", lags, "]"))
    }
    paste0(terms, collapse = " +\n              ")
  }

  jagscode = "
  # Apply GARMA recursion to link-scale residuals
  resid_arma_[1] = 0"

  if (max_order >= 2) {
    for (i in 2:max_order) {
      jagscode = paste0(
        jagscode, "\n  resid_arma_[", i, "] = ",
        get_terms(i, min(ar_order, i - 1), min(ma_order, i - 1))
      )
    }
  }

  paste0(
    jagscode,
    "\n  for (i_ in ", max_order + 1, ":length(", x_name, ")) {",
    "\n    resid_arma_[i_] = ", get_terms("i_", ar_order, ma_order),
    "\n  }"
  )
}


# Backwards-compatible internal wrapper for pure AR code generation.
get_ar_jagscode = function(ar_order, x_name) {
  get_arma_jagscode(ar_order, NA, x_name)
}
