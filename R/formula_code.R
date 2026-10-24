# ABOUT: The linear predictor of each distributional parameter as a table of terms. 
# - evaluate_formula() is the R version of the mcp model. Used for fitted(), predict(), fit$simulate(), etc.
# - get_formula_jags() generates the JAGS version of the mcp model
# -----------------


# The terms of each dpar's linear predictor: one row per coefficient in each segment where it
# applies (including reused coefficients and group-level deviations) and one row per offset.
# A term applies from change point `start`, which starts its segment, until change point `end`
# (NA: through the last segment). Change points are numbered from cp_0 = min(x) to cp_K = max(x).
# `dpar` names the linear predictor, e.g., "mu", "sigma", or "ar1" for the first AR lag.
# - model_tables: Model tables, as in `get_fit_model_tables(fit)`.
get_formula_terms = function(model_tables) {
  columns = c("dpar", "order", "segment", "next_segment", "x_degree", "matrix_col", "name")

  # Offsets apply in their own segment only
  offsets = Filter(function(design_spec) isTRUE(design_spec$has_offset), model_tables$design_specs)
  field = function(name, type) unname(vapply(offsets, function(design_spec) design_spec[[name]], type))
  segment = field("segment", integer(1))
  offsets = tibble::tibble(
    dpar = field("dpar", character(1)),
    order = field("order", integer(1)),
    segment = segment,
    next_segment = ifelse(segment < nrow(model_tables$segments), segment + 1L, NA_integer_),
    x_degree = 0L,
    offset = field("offset_name", character(1))
  )

  terms = dplyr::bind_rows(model_tables$predictors[columns], model_tables$group_predictors[c(columns, "group_col")], offsets)
  terms$dpar = paste0(terms$dpar, tidyr::replace_na(as.character(terms$order), ""))
  terms$start = terms$segment - 1L
  terms$end = terms$next_segment - 1L
  terms[c("dpar", "segment", "start", "end", "x_degree", "matrix_col", "name", "group_col", "offset")]
}


# Exponent for x-terms of degree `x_degree` in JAGS code, e.g., "" for x and "^3" for x^3.
x_power = function(x_degree) ifelse(x_degree > 1, paste0("^", x_degree), "")


#' Write the JAGS code for the linear predictors of all dpars
#'
#' @param model_tables Model tables, as in `get_fit_model_tables(fit)`.
#' @param family An mcpfamily object.
#' @param nodes Named vector of the JAGS node that is sampled for each parameter (see `add_nodes()`).
#' @return A string with JAGS code for the body of the likelihood loop.
#' @noRd
get_formula_jags = function(model_tables, family, nodes) {
  par_x = model_tables$data_columns$par_x
  terms = get_formula_terms(model_tables)
  K = nrow(model_tables$segments)

  # Change points cp_0, ..., cp_K. Group-specific change points are sampled locations or include group-level deviations.
  cps = model_tables$cps
  index = paste0("[", cps$group_col, "_index_[i_]]")
  node = unname(nodes[cps$group_name])
  is_location = !is.na(node) & node != cps$group_name
  cp = c(
    "cp_0",
    ifelse(is.na(cps$group_name), cps$name, ifelse(
      is_location,
      paste0(node, index),
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

  # One sum per dpar. Coefficients sharing a segment, x-degree, and end are one inprod().
  formula_dpar = function(dpar) {
    dpar_terms = terms[terms$dpar == dpar, ]
    coefficients = dpar_terms[is.na(dpar_terms$offset), ] |>
      dplyr::mutate(
        name = ifelse(is.na(.data$group_col), .data$name, paste0(.data$name, "[", .data$group_col, "_index_[i_]]")),
        x_code = ifelse(.data$x_degree == 0, "1", paste0("x_local_", .data$segment, "_[i_]", x_power(.data$x_degree)))
      ) |>
      dplyr::group_by(.data$segment, .data$x_degree, .data$end) |>
      dplyr::summarise(
        code = paste0(
          "  ", applies(dplyr::first(.data$start), dplyr::first(.data$end)),
          " * inprod(rhs_matrix_[i_, c(", paste0(.data$matrix_col, collapse = ", "), ")], c(", paste0(.data$name, collapse = ", "), "))",
          " * ", dplyr::first(.data$x_code)
        ),
        .groups = "drop"
      )
    offsets = dpar_terms[!is.na(dpar_terms$offset), ]
    offset_code = if (nrow(offsets) > 0) paste0("  ", applies(offsets$start, offsets$end), " * ", offsets$offset, "[i_]")
    node = if (dpar %in% family$dpar_specs$dpar) paste0("link_", dpar) else dpar
    paste0("# Formula for ", dpar, "\n", node, "_[i_] =\n", paste0(c(coefficients$code, offset_code), collapse = " + \n"))
  }
  formula_jags = paste0(
    x_local, get_garma_threshold_jags_code(model_tables$predictors), "\n\n",
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
#' @param model_tables Model tables, as in `get_fit_model_tables(fit)`.
#' @param family An mcpfamily object.
#' @param data Model data from `get_model_data()`.
#' @param pars Named list of parameter values. Each is one value per draw (a vector) or one
#'   value per row and draw (a rows x draws matrix), e.g., group-level deviations.
#' @param ndraws Number of draws.
#' @return A named list of rows x draws matrices: one for each dpar, e.g., `mu`, `sigma`, and `ar1`.
#' @noRd
evaluate_formula = function(model_tables, family, data, pars, ndraws) {
  terms = get_formula_terms(model_tables)
  x = data$x
  n_rows = length(x)
  # A rows x draws matrix of the parameter `name` times `column`, which has one value per row
  value = function(name, column = rep(1, n_rows)) {
    par_values = pars[[name]]
    if (is.matrix(par_values)) column * par_values else outer(column, par_values)
  }

  # Change points cp_0, ..., cp_K, including group-level deviations. Outer segments extend beyond the data.
  cps = model_tables$cps
  cp = c(-Inf, lapply(seq_len(nrow(cps)), function(k) {
    if (is.na(cps$group_name[k])) value(cps$name[k]) else value(cps$name[k]) + value(cps$group_name[k])
  }), Inf)
  cp_at = function(k) cp[[k + 1]]

  # x measured from each segment's change point
  x_local = lapply(seq_len(nrow(model_tables$segments)), function(segment) {
    pmin(x, cp_at(segment)) - if (segment > 1) cp_at(segment - 1) else 0
  })

  # Sum the terms of each dpar. As in the JAGS code, terms sharing a segment, end, and x-degree
  # are summed before multiplying by x and by where they apply. Dpars without terms are zero.
  dpars = unique(c(terms$dpar, family$dpar_specs$dpar))
  lapply(stats::setNames(nm = dpars), function(dpar) {
    link = matrix(0, n_rows, ndraws)
    dpar_terms = terms[terms$dpar == dpar, ]
    for (inprod_terms in split(dpar_terms, paste(dpar_terms$segment, dpar_terms$end, dpar_terms$x_degree))) {
      inprod = 0
      for (i in seq_len(nrow(inprod_terms)))
        inprod = inprod + if (is.na(inprod_terms$offset[i])) value(inprod_terms$name[i], data$matrix[, inprod_terms$matrix_col[i]]) else data$offsets[[inprod_terms$offset[i]]]
      if (inprod_terms$x_degree[1] > 0)
        inprod = inprod * x_local[[inprod_terms$segment[1]]]^inprod_terms$x_degree[1]
      applies = (x >= cp_at(inprod_terms$start[1])) * (if (is.na(inprod_terms$end[1])) 1 else x < cp_at(inprod_terms$end[1]))
      link = link + applies * inprod
    }
    link
  })
}
