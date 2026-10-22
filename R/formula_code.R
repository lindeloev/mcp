# ABOUT: The linear predictor of each distributional parameter as a table of terms. 
# - evaluate_formula() is the R version of the mcp model. Used for fitted(), predict(), fit$simulate(), etc.
# - get_formula_jags() generates the JAGS version of the mcp model
# -----------------


# The terms of each dpar's linear predictor: one row per coefficient in each segment where it
# applies (including reused coefficients and group-level deviations) and one row per offset.
# A term applies from change point `start`, which starts its segment, until change point `end`
# (NA: through the last segment). Change points are numbered from cp_0 = min(x) to cp_K = max(x).
# `dpar` names the linear predictor, e.g., "mu", "sigma", or "ar1" for the first AR lag.
# - tables: Model tables, as in `get_fit_model_tables(fit)`.
get_formula_terms = function(tables) {
  columns = c("dpar", "order", "segment", "next_segment", "x_factor", "matrix_col", "name")
  group_effects = tables$group_effects[tables$group_effects$part == "predictor", ]

  # Offsets apply in their own segment only
  offsets = Filter(function(spec) isTRUE(spec$has_offset), tables$design_specs)
  field = function(name, type) unname(vapply(offsets, function(spec) spec[[name]], type))
  offsets = tibble::tibble(
    dpar = field("dpar", character(1)),
    order = field("order", integer(1)),
    segment = field("segment", integer(1)),
    next_segment = ifelse(segment < nrow(tables$segments), segment + 1L, NA_integer_),
    x_factor = "1",
    offset = field("offset_name", character(1))
  )

  terms = dplyr::bind_rows(tables$predictors[columns], group_effects[c(columns, "group_col")], offsets)
  terms$dpar = paste0(terms$dpar, tidyr::replace_na(as.character(terms$order), ""))
  terms$start = terms$segment - 1L
  terms$end = terms$next_segment - 1L
  terms[c("dpar", "segment", "start", "end", "x_factor", "matrix_col", "name", "group_col", "offset")]
}


# Power of x (measured from the change point) in an `x_factor`: "1" is 0, "x" is 1, and "x^3" is 3.
x_degree = function(x_factor) {
  as.integer(ifelse(x_factor == "1", "0", ifelse(x_factor == "x", "1", sub("x^", "", x_factor, fixed = TRUE))))
}


#' Write the JAGS code for the linear predictors of all dpars
#'
#' @param tables Model tables, as in `get_fit_model_tables(fit)`.
#' @param family An mcpfamily object.
#' @param locations Group-level change points sampled as absolute locations (see `is_cp_location()`).
#' @return A string with JAGS code for the body of the likelihood loop.
#' @noRd
get_formula_jags = function(tables, family, locations) {
  par_x = tables$data_columns$par_x
  terms = get_formula_terms(tables)
  K = nrow(tables$segments)

  # Change points cp_0, ..., cp_K. Group-specific change points are sampled locations or include group-level deviations.
  cps = tables$cps
  index = paste0("[", cps$group_col, "_index_[i_]]")
  cp = c(
    "cp_0",
    ifelse(is.na(cps$group_name), cps$name, ifelse(
      cps$group_name %in% locations,
      paste0(cps$group_name, "_location", index),
      paste0("(", cps$name, " + ", cps$group_name, index, ")")
    )),
    paste0("cp_", K)
  )
  cp_at = function(k) cp[k + 1]
  x_at = function(op, k) paste0("(", par_x, "[i_] ", op, " ", cp_at(k), ")")
  applies = function(start, end) paste0(x_at(">=", start), ifelse(is.na(end), "", paste0(" * ", x_at("<", end))))

  # x measured from each segment's change point
  x_local = paste0(
    "\n# ", par_x, " measured from each segment's change point",
    paste0("\nx_local_", seq_len(K), "_[i_] = min(", par_x, "[i_], ", cp_at(seq_len(K)), ")",
           ifelse(seq_len(K) > 1, paste0(" - ", cp_at(seq_len(K) - 1)), ""), collapse = "")
  )

  # One sum per dpar. Coefficients sharing a segment, x-factor, and end are one inprod().
  formula_dpar = function(dpar) {
    dpar_terms = terms[terms$dpar == dpar, ]
    coefficients = dpar_terms[is.na(dpar_terms$offset), ] %>%
      dplyr::mutate(name = ifelse(is.na(.data$group_col), .data$name, paste0(.data$name, "[", .data$group_col, "_index_[i_]]"))) %>%
      dplyr::group_by(.data$segment, .data$x_factor, .data$end) %>%
      dplyr::summarise(
        code = paste0(
          "  ", applies(dplyr::first(.data$start), dplyr::first(.data$end)),
          " * inprod(rhs_matrix_[i_, c(", paste0(.data$matrix_col, collapse = ", "), ")], c(", paste0(.data$name, collapse = ", "), "))",
          " * ", sub("x", paste0("x_local_", dplyr::first(.data$segment), "_[i_]"), dplyr::first(.data$x_factor), fixed = TRUE)
        ),
        .groups = "drop"
      )
    offsets = dpar_terms[!is.na(dpar_terms$offset), ]
    offset_code = if (nrow(offsets) > 0) paste0("  ", applies(offsets$start, offsets$end), " * ", offsets$offset, "[i_]")
    node = if (dpar %in% family$dpar_specs$dpar) paste0("link_", dpar) else dpar
    paste0("# Formula for ", dpar, "\n", node, "_[i_] =\n", paste0(c(coefficients$code, offset_code), collapse = " + \n"))
  }
  formula_jags = paste0(
    x_local, get_garma_threshold_jagscode(tables$predictors), "\n\n",
    paste0(vapply(unique(terms$dpar), formula_dpar, character(1)), collapse = "\n\n")
  )

  # Dpars without terms in any segment (all ~ 0)
  for (dpar in setdiff(family$dpar_specs$dpar, terms$dpar))
    formula_jags = paste0(formula_jags, "\n\n# All segments are ~ 0 for this par:\nlink_", dpar, "_[i_] = 0")

  class(formula_jags) = c("mcptext", "character")  # Nicer printing
  formula_jags
}


#' Evaluate the linear predictors of all dpars in R
#'
#' The R counterpart of `get_formula_jags()`, evaluated for all rows of data and draws of
#' parameter values at once.
#'
#' @param fit An `mcpfit` object.
#' @param data Model data from `get_model_data()`.
#' @param pars Named list of parameter values. Each is one value per draw (a vector) or one
#'   value per row and draw (a rows x draws matrix), e.g., group-level deviations.
#' @param n_draws Number of draws.
#' @return A named list of rows x draws matrices: one for each dpar, e.g., `mu`, `sigma`, and `ar1`.
#' @noRd
evaluate_formula = function(fit, data, pars, n_draws) {
  tables = get_fit_model_tables(fit)
  terms = get_formula_terms(tables)
  x = data$x
  n_rows = length(x)
  # A rows x draws matrix of the parameter `name` times `column`, which has one value per row
  value = function(name, column = rep(1, n_rows)) {
    v = pars[[name]]
    if (is.matrix(v)) column * v else outer(column, v)
  }

  # Change points cp_0, ..., cp_K, including group-level deviations. Outer segments extend beyond the data.
  cps = tables$cps
  cp = c(-Inf, lapply(seq_len(nrow(cps)), function(k) {
    if (is.na(cps$group_name[k])) value(cps$name[k]) else value(cps$name[k]) + value(cps$group_name[k])
  }), Inf)
  cp_at = function(k) cp[[k + 1]]

  # x measured from each segment's change point
  x_local = lapply(seq_len(nrow(tables$segments)), function(s) {
    pmin(x, cp_at(s)) - if (s > 1) cp_at(s - 1) else 0
  })

  # Sum the terms of each dpar. As in the JAGS code, terms sharing a segment, end, and x-factor
  # are summed before multiplying by x and by where they apply. Dpars without terms are zero.
  dpars = unique(c(terms$dpar, fit$family$dpars))
  lapply(stats::setNames(nm = dpars), function(dpar) {
    link = matrix(0, n_rows, n_draws)
    dpar_terms = terms[terms$dpar == dpar, ]
    for (group in split(dpar_terms, paste(dpar_terms$segment, dpar_terms$end, dpar_terms$x_factor))) {
      inprod = 0
      for (i in seq_len(nrow(group)))
        inprod = inprod + if (is.na(group$offset[i])) value(group$name[i], data$matrix[, group$matrix_col[i]]) else data$offsets[[group$offset[i]]]
      if (group$x_factor[1] != "1")
        inprod = inprod * x_local[[group$segment[1]]]^x_degree(group$x_factor[1])
      applies = (x >= cp_at(group$start[1])) * (if (is.na(group$end[1])) 1 else x < cp_at(group$end[1]))
      link = link + applies * inprod
    }
    link
  })
}
