#' Get the JAGS Model, Data, and Initial Values of an mcp Model
#'
#' `r lifecycle::badge("experimental")`
#' 
#' Returns everything needed to sample an `mcp` model directly in JAGS. This gives more flexibility 
#' to modify the JAGS code or to monitor other nodes. The downside is that `mcp` methods (`summary()`, `plot()`, `loo()`, etc.) 
#' cannot use draws sampled this way. Use [mcpfamily()] to fit other likelihoods in `mcp`.
#'
#' @param fit An [mcpfit] object, e.g., from `mcp(..., sample = FALSE)`.
#' @param prior Logical. Return data and initial values for sampling the prior instead of the posterior,
#'   as in `mcp(..., sample = "prior")`.
#' @return An `mcpjags` object, which prints the JAGS code. It is a list with:
#'   * `code`: String. The JAGS model.
#'   * `data`: List with the data.
#'   * `inits`: List with initial values for each chain, including JAGS random-number seeds. These are the initial
#'     values that `mcp()` used to sample `fit`. The seeds are random for fits with `sample = FALSE` and `seed = NULL`.
#'   * `variable.names`: Character vector. The nodes that `mcp()` monitors, i.e., `mcp_pars(fit)$parameter`.
#'   * `modules`: Character vector. The JAGS modules that `mcp()` loads, e.g., `"glm"`.
#'   * `levels`: Named list with the levels of each grouping variable. JAGS indexes them as 1, 2, etc. in
#'     this order, e.g., `x_1_id[2]` is the deviation for `levels$id[2]`.
#' @aliases as_jags mcpjags
#' @export
#' @examples
#' data = mcp_example_data("demo")
#' fit = mcp(list(response ~ 1, ~ 0 + time, ~ 1 + time), data, sample = FALSE)
#' jags = as_jags(fit)
#' jags  # Show the JAGS code. Also in jags$code
#'
#' \donttest{
#' # After making changes to jags$code or other elements, sample it using rjags
#' for (module in jags$modules)
#'   rjags::load.module(module)
#' model = rjags::jags.model(textConnection(jags$code), jags$data, jags$inits, n.chains = length(jags$inits))
#' draws = rjags::coda.samples(model, jags$variable.names, n.iter = 1000)
#' summary(draws)  # This calls coda:::summary.mcmc.list() - not mcp:::summary.mcpfit()
#' }
as_jags = function(fit, prior = FALSE) {
  # Check inputs and extract info
  checkmate::assert_flag(prior)
  model_tables = get_fit_model_tables(fit)
  prior_table = fit$.internal$prior_table
  sampling = fit$.internal$sampling
  if (is.null(sampling))
    stop("This `mcpfit` lacks the settings that `mcp()` sampled it with. Refit it with mcp >= 0.4.0.9000.", call. = FALSE)

  # Build objects for return
  code = .subset2(fit, "jags_code")
  custom_code = has_custom_jags_code(fit)
  data = get_jags_data(fit$data, fit$family, model_tables, code, generated = !custom_code)
  if (prior)
    data = get_jags_prior_data(data, code, model_tables$data_columns)
  sample = if (prior) "prior" else "post"
  seed = if (is.null(sampling$seed)) sample.int(.Machine$integer.max, 1) else sampling$seed  # Empty fits without `seed`
  group_cols = unique(stats::na.omit(model_tables$parameters$group_col))

  # Return
  out = list(
    code = code,
    data = data,
    inits = get_chain_inits(model_tables, prior_table, fit$family, fit$data, sample, sampling$chains, sampling$inits, seed, if (custom_code) code),
    variable.names = prior_table$parameter,
    modules = if (use_jags_glm(fit$family, definition_rows(model_tables$group_coefs))) "glm" else character(),
    levels = stats::setNames(lapply(group_cols, function(group_col) levels(factor(fit$data[[group_col]]))), group_cols)
  )
  class(out) = "mcpjags"
  out
}


#' @export
#' @noRd
print.mcpjags = function(x, ...) {
  rlang::check_dots_empty()
  cat(x$code, "\n\n", sep = "")
  cat("# Also contains $data, $inits (", length(x$inits), " chains), $variable.names, $modules, and $levels. See ?as_jags\n", sep = "")
  invisible(x)
}
