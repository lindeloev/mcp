#' RowSums of element-wise products.
#' Like an inner product, just vectorized for many y.
#' This ensures R <--> JAGS code compatibility
#'
#' @aliases inprod
#' @param x A matrix
#' @param y A matrix
#' @return A vector
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
inprod = function(x, y) {
  assert_types(x, "matrix")
  assert_types(y, "matrix")
  rowSums(x * y)
}


# Converts logical(0) to null. Returns x otherwise
logical0_to_null = function(x) {
  if (length(x) > 0)
    return(x)
  else
    return(NULL)
}


# if((a %in% b) == FALSE) --> if(a %notin% b)
`%notin%` = Negate(`%in%`)


# List of categorical column names and their unique levels
get_categorical_levels = function(df) {
  assert_types(df, "data.frame", "tibble")
  categorical_cols = colnames(df)[sapply(df, class) %in% c("factor", "logical", "character")]
  lapply(df[, categorical_cols, drop = FALSE], unique)
}


# Ask reminder questions for CRAN export
release_questions = function() {
  c(
    #"TEST: Have you run the extensive tests? options(test_mcp_allmodels = TRUE)",
    "TEST: Have you run the test of fits? (uncomment skip() in test-fits-examples.R and helper-fits.R)",
    "TEST: Have you run `revdepcheck::revdep_check()`?",

    "DOC: Have you built the README plots and checked them? source('vignettes/figures/make_README_plots.R')",
    "DOC: Have you re-built the site using pkgdown::build_site() AFTER deleting caches of articles in 'vignettes/*_cache/'?",
    "DOC: Have you checked all articles and plots after re-building the site?",
    "DOC: Have you run the script to insert the correct logo.png in the HTML meta?"
  )
}


# Returns the AR order or NA if no AR
get_ar_order = function(rhs_table) {
  ar_order = ifelse("ar" %in% rhs_table$dpar, max(rhs_table$order, na.rm = TRUE), NA)
}


#' Remove varying or population terms from a formula
#'
#' WARNING: removes response side from the formula
#'
#' @aliases remove_terms
#' @keywords internal
#' @noRd
#' @param form A formula
#' @param remove Either "varying" or "population". These are removed.
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
remove_terms = function(form, remove) {
  assert_types(form, "formula", len = c(2, 3))
  assert_value(remove, allowed = c("varying", "population"))

  # Find terms with "|"
  attrs = attributes(stats::terms(form))
  term.labels = attrs$term.labels
  varying_bool = stringr::str_detect(term.labels, "\\|")

  # Add parenthesis back to them
  term.labels[varying_bool] = paste0("(", term.labels[varying_bool], ")")

  # Remove non-matching types
  if (remove == "varying") {
    term.labels = term.labels[!varying_bool]
    term.labels = c(attrs$intercept, term.labels)  # Add intercept indicator
  } else if (remove == "population") {
    term.labels = term.labels[varying_bool]
  }

  # Build formula from terms and return
  if (length(term.labels) == 0) {
    return(NULL)
  } else {
    formula_terms = paste0(term.labels, collapse = " + ")
    formula_str = paste0("~", formula_terms)
    return(stats::as.formula(formula_str, env=globalenv()))
  }
}


#' Takes any formula-like input and returns a formula
#' @aliases to_formula
#' @keywords internal
#' @noRd
#' @param form Formula or character (with or without initial tilde/"~")
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
to_formula = function(form) {
  assert_types(form, "character", "formula", len = c(1, 3))
  if (is.character(form)) {
    # Add tilde
    if (!stringr::str_detect(form, "^(\\s|)~")) {
      form = paste0("~", form)
    }
    form = stats::as.formula(form)
  }

  form
}


#' Homogonize enumerating strings in mcp
#'
#' Nice for error messages.
#'
#' @aliases collapse_quote
#' @keywords internal
#' @noRd
#' @param x A character vector
#' @return Character
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
and_collapse = function(x) {
  assert_types(x, "character")
  paste0(x, collapse = " and ")
}


#' Converts formula to string
#'
#' Note: this uses base R and circumvents the length-limitation of `deparse()`
#' and `format()`.
#'
#' @aliases formula_to_char
#' @keywords internal
#' @noRd
#' @param form Any valid formula with any number of tildes.
#' @return A character.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
formula_to_char = function(form) {
  assert_types(form, "formula", len = c(1, 3))
  form_char = as.character(form)
  if (length(form_char) == 2 & form_char[1] == "~") {
    return(paste0(form_char, collapse = " "))
  } else if (length(form_char == 3) & form_char[1] == "~") {
    return(paste0(form_char[c(2, 3)], collapse = " ~ "))
  } else {
    stop_github("Could not decode formula ", deparse(form, width.cutoff = 500))
  }
}


#' Returns the right-hand-side of a formula
#'
#' @aliases get_rhs
#' @keywords internal
#' @noRd
#' @param form Formula, e.g. `~x`, `y ~ x` or `y ~ z ~ x`
#' @return A formula
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs = function(form) {
  assert_types(form, "formula")
  if (length(form) == 2) {
    return(form)
  } else if (length(form) == 3) {
    return(form[-2])
  }
}


#' Returns all vars in the RHS of an mcpmodel
#'
#' @aliases get_rhs_vars
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return Character vector with unique term names
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs_vars = function(model) {
  assert_types(model, "mcpmodel")

  model %>%
    lapply(get_rhs) %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()
}

#' Returns all vars in the RHS of an mcpmodel
#'
#' @aliases get_model_vars
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @return Character vector with unique term names
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_model_vars = function(model) {
  assert_types(model, "mcpmodel")

  model %>%
    lapply(all.vars) %>%
    unlist() %>%
    unique()
}


#' Create model matrix from rhs_table
#'
#' cbinds rhs_table$matrix_data
#' @aliases get_rhs_matrix
#' @keywords internal
#' @noRd
#' @param rhs_table The output of `get_rhs_table()`
#' @return matrix with one column for each row in `rhs_table`
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_rhs_matrix = function(rhs_table) {
  assert_types(rhs_table, "data.frame", "tibble")
  suppressMessages(dplyr::bind_cols(rhs_table$matrix_data, .name_repair = "unique")) %>% # Suppress message about lacking column names
    as.matrix() %>%
    magrittr::set_colnames(rhs_table$code_name)
}


#' Convert from tidy to matrix
#'
#' Converts from the output of `tidy_samples()` or `pp_eval(fit, samples_format = "tidy")`
#' to an `N_draws` X `nrows(newdata)` matrix with fitted/predicted values. This format is
#' used y `brms` and it's useful as `yrep` in `bayesplot::ppc_*` functions.
#'
#' @aliases tidy_to_matrix
#' @keywords internal
#' @noRd
#' @param samples Samples in tidy format
#' @param returnvar An `rlang::sym()` object.
#' @return An  `N_draws` X `nrows(newdata)` matrix.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
tidy_to_matrix = function(samples, returnvar) {
  returnvar = rlang::sym(returnvar)
  samples %>%
    dplyr::select(".draw", "data_row", {{ returnvar }}) %>%
    tidyr::pivot_wider(names_from = "data_row", values_from = {{ returnvar }}) %>%
    dplyr::select(-".draw") %>%
    as.matrix()
}


#' Expand samples with quantiles
#'
#' TO DO: implement using `fitted()` and `predict()` but avoid double-computing the samples? E.g.:
#' `get_quantiles2 = function(fit, quantiles, facet_by = NULL) {`
#'   `fitted(fit, probs = c(0.1, 0.5, 0.9), newdata = data.frame(x = c(11, 50, 100))) %>%`
#'   `tidyr::pivot_longer(tidyselect::starts_with("Q")) %>%`
#'   `dplyr::mutate(quantile = stringr::str_remove(name, "Q") %>% as.numeric() / 100)`
#' `}`
#'
#' @aliases get_quantiles
#' @keywords internal
#' @noRd
#' @inheritParams plot.mcpfit
#' @param samples A tidybayes tibble
#' @param quantiles Vector of quantiles (0.0 to 1.0)
#' @param xvar An rlang::sym() with the name of the x-col in `samples`
#' @param yvar An rlang::sym() with the name of the response col in `samples`
#' @return A tidybayes long format tibble with the column "quantile"
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_quantiles = function(samples, quantiles, xvar, yvar, facet_by = NULL) {
  # Trick to declare no facet = common group for all
  facet_by = logical0_to_null(facet_by)
  if (!is.null(facet_by))
    facet_by = rlang::sym(facet_by)

  # Return data with added quantiles
  samples %>%
    tidyr::expand_grid(quantile = quantiles) %>%

    # Now compute the quantile for each parameter, quantile, and (optionally) facet:
    dplyr::group_by(!!xvar, .data$quantile, !!facet_by) %>%
    #dplyr::group_by(!!rlang::sym(facet_by), .add = TRUE) %>%
    dplyr::summarise(
      y = stats::quantile(!!yvar, probs = .data$quantile[1])
    )
}


#' Print mcplist
#'
#' Shows a list in a more condensed format using `str(list)`.
#' @aliases print.mcplist
#' @inheritParams print.mcpfit
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
print.mcplist = function(x, ...) {
  assert_types(x, "list")
  assert_ellipsis(...)

  # For all-formula list (typically fit$model)
  if (all(sapply(x, is.formula)) == TRUE) {
    cat(paste0("List of ", length(x), "\n"))
    for (i in x) {
      cat(" $ ")
      print(i)
    }
  } else {
    # Other lists
    utils::str(x, vec.len = Inf, give.head = FALSE, give.attr = FALSE)
  }
}


#' Nice Printing of Multiline Texts
#'
#' Useful for `print(fit$jags_code)`, `print(mcp_demo$call)`, etc.
#'
#' @aliases print.mcptext
#' @param x Character, often with newlines.
#' @param ... Currently ignored.
#' @return NULL
#' @export
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @examples
#' mytext = "line1 = 2\n line2 = 'horse'"
#' class(mytext) = "mcptext"
#' print(mytext)
print.mcptext = function(x, ...) {
  assert_types(x, "character", len = 1)
  assert_ellipsis(...)
  cat(x)
}


# Set model environment to parent.frame() for prettier printing
# and because it was created in a different environment than inteded for use.
fix_model_environment = function(model) {
  assert_types(model, "mcpmodel")
  for (i in seq_along(model))
    environment(model[[i]]) = globalenv()
  model
}
