# ABOUT: These functions build the change point regression model for R and JAGS
# -----------------


#' Call `get_formula_jags_dpar` for multiple dpars and paste strings
#'
#' @aliases get_formula_jags
#' @keywords internal
#' @inheritParams get_formula_jags_dpar
#' @inheritParams mcp
#' @param dpars A character vector of dpars to including in model building
#' @return A string with JAGS code.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_formula_jags = function(ST, rhs_table, par_x, family) {
  # Add X-helpers which code the X relative to the start of each segment.
  local_x_str = "\n# par_x local to each segment"
  for (i in seq_len(nrow(ST))) {
    segment_start = ifelse(i > 1, yes = paste0(" - ", ST$cp_code_form[i], ""), no = "")  #
    segment_end = ifelse(i < nrow(ST), yes = paste0(ST$cp_code_form[i + 1], ""), no = paste0("cp_", i))  # infinite if last segment.

    local_x_str = paste0(local_x_str, "\n", par_x, "_local_", i, "_[i_] = min(", par_x, "[i_], ", segment_end, ")", segment_start)
  }

  # Build formula for each dpar (note plural "_dpars")
  formula_jags_dpars = rhs_table %>%
    dplyr::select(-matrix_data) %>%  # Throw less data around
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Left-join ST
      this_cp = ST$cp_code_form[[segment]],
      next_cp = ST$cp_code_form[next_intercept],  # replace segment with code string
      form = ST$form[[segment]],

      # One dpar per ar order: (ar, 1) --> ar1
      dpar = paste0(dpar, tidyr::replace_na(order, ""))
    ) %>%

    # Build formula for each dpar
    dplyr::group_by(dpar) %>%
    tidyr::nest() %>%
    dplyr::rowwise() %>%
    dplyr::summarise(
      formula_jags_dpar = get_formula_jags_dpar(data, dpar, par_x)
    ) %>%
    dplyr::pull(formula_jags_dpar) %>%
    paste0(collapse = "\n\n")

  # Concatenate and return
  formula_jags = paste0(local_x_str, "\n\n", formula_jags_dpars)

  # Special case when no terms are present for a given dpar (all ~0): insert "dpar = 0".
  for (dpar in family$dpars) {
    if (!(dpar %in% rhs_table$dpar))
      formula_jags = paste0(formula_jags, "\n\n# All segments are ~ 0 for this par:\n", dpar, "_[i_] = 0")
  }

  class(formula_jags) = c("mcptext", "character")  # Nicer printing
  return(formula_jags)
}


#' Build an R formula (as string) for a dpar
#'
#' @aliases get_formula_jags_dpar
#' @keywords internal
#' @inheritParams mcp
#' @param ST Tibble. Returned by `get_segment_table`.
#' @param dpar_table A rhs_table with only one (dpar, order) combo
#' @return A string with JAGS code.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_formula_jags_dpar = function(dpar_table, dpar, par_x) {
  # Build this! Initiate
  formula_str = paste0("# Formula for ", dpar)
  if (dpar == "sigma") {
    formula_str = paste0(formula_str, "\nsigma_tmp_[i_] =\n")
  } else {
    formula_str = paste0(formula_str, "\n", dpar, "_[i_] =\n")
  }


  df_code_strs = dpar_table %>%
    # Code dpar_data_segment column indices
    dplyr::arrange(segment) %>%  # just to make sure
    dplyr::group_by(segment) %>%

    # Summarise
    dplyr::group_by(segment, x_factor) %>%
    dplyr::summarise(
      # The parts
      indicator_this = paste0("  (", par_x, "[i_] >= ", dplyr::first(this_cp), ")"),
      indicator_next = dplyr::if_else(is.na(dplyr::first(next_cp)) == TRUE, "", paste0(" * (", par_x, "[i_] < ", dplyr::first(next_cp), ")")),
      inprod = paste0(" * inprod(rhs_matrix_[i_, c(", paste0(matrix_col, collapse = ", "), ")], c(", paste0(code_name, collapse = ", "), "))"),
      x_factor = gsub("x", paste0(par_x, "_local_", dplyr::first(segment), "_[i_]"), dplyr::first(x_factor)),

      # All together
      segment_code = paste0(indicator_this, indicator_next, inprod, " * ", x_factor),
      form = dplyr::first(form)
    ) %>%

    # Add title-comment
    dplyr::group_by(segment) %>%
    dplyr::mutate(
      title = dplyr::if_else(dplyr::row_number() == 1, paste0("\n  # Segment ", dplyr::first(segment), ": ", dplyr::first(form), "\n"), ""),
      segment_code = paste0(title, segment_code)
    )

  # Return
  all_predictors = paste0(df_code_strs$segment_code, collapse = " + \n")
  formula_str = paste0(formula_str, all_predictors)
  if (dpar == "sigma")
    formula_str = paste0(formula_str, "\nsigma_[i_] = max(10^-9, sigma_tmp_[i_])  # Count negative sigma as just-above-zero sigma")
  return(formula_str)
}


#' Convert for-looped JAGS code to vectorized R code
#'
#' @aliases get_formula_r
#' @keywords internal
#' @param formula_jags Character, often residing in `fit$.internal$formula_jags`.
#' @param rhs_table Output of `get_rhs_table()`
#' @param pars The list that ends up in `fit$pars`
#' @return Character
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_formula_r = function(formula_jags, rhs_table, pars) {
  param_pars = get_sim_pars(rhs_table, pars)

  # Replacements that turns rowwise JAGS code into vectorized R code
  replace_args = c(
    setNames(paste0(", args$", param_pars), paste0(", ", param_pars)),
    setNames(paste0("args$", param_pars, ", "), paste0(param_pars, ", ")),
    setNames(paste0(" (args$", param_pars, " + "), paste0(" (", param_pars, " + ")),  # varying change point
    setNames(paste0(" + args$", param_pars, ")"), paste0(" + ", param_pars, ")")),  # varying change point
    setNames(paste0("cbind(args$"), "cbind("),
    setNames("args$", "args$args$"),  # Fix double-inserting args$ above
    setNames(paste0("args$", pars$x, " >="), paste0(pars$x, " >=")),
    setNames(paste0("args$", pars$x, " <"), paste0(pars$x, " <")),
    setNames("pmin(args$", paste0("pmin("))
  )
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
  return(formula_r)
}


#' Get JAGS code to model autoregressive effects
#'
#' This is simply code for `resid_arma_`.
#'
#' @aliases get_ar_jagscode
#' @keywords internal
#' @param ar_order Positive integer.
#' @param x_name Character. Name of some vector that has the length of the dataset.
#' @return Character JAGS code
#' @seealso simulate_ar
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_ar_jagscode = function(ar_order, x_name) {
  jagscode = "
  # Apply autoregression to the residuals
  resid_arma_[1] = 0"

  # For data points lower than the full order
  if (ar_order >= 2) {
    for (i in 2:ar_order) {
      jagscode = paste0(jagscode, "
  resid_arma_[", i, "] = ", paste0("ar", 1:(i-1), "_[", i, " - ", 1:(i-1), "] * resid_abs_[", i, " - ", 1:(i-1), "]", collapse = " +\n              "))
    }
  }

  # For full order
  jagscode = paste0(jagscode, "
  for (i_ in ", ar_order + 1, ":length(", x_name, ")) {
    resid_arma_[i_] = 0")
  for (i in seq_len(ar_order)) {
    jagscode = paste0(jagscode, " + \n      ar", i, "_[i_] * resid_abs_[i_ - ", i, "]")
  }

  # Finish up
  jagscode = paste0(jagscode, "
  }")

  return(jagscode)
}
