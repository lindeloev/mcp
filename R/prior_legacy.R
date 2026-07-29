# Legacy prior compatibility -----------------------------------------------

# This is the complete translation layer for data constants released before
# readable prior expressions. Keeping token, value, and replacement together
# makes the compatibility contract explicit and prevents the three views from
# drifting apart.
legacy_prior_map = function(context) {
  x = context$data[[context$x_name]]
  y = context$data[[context$y_name]]
  tibble::tibble(
    token = c("MINX", "MAXX", "MEANX", "SDX", "MINY", "MAXY", "MEANY", "SDY", "N_CP"),
    value = c(
      min(x, na.rm = TRUE),
      max(x, na.rm = TRUE),
      mean(x, na.rm = TRUE),
      stats::sd(x, na.rm = TRUE),
      min(y, na.rm = TRUE),
      max(y, na.rm = TRUE),
      mean(y, na.rm = TRUE),
      stats::sd(y, na.rm = TRUE),
      context$n_cp
    ),
    rule = c(
      paste0("min(", context$x_display, ")"),
      paste0("max(", context$x_display, ")"),
      paste0("mean(", context$x_display, ")"),
      paste0("sd(", context$x_display, ")"),
      paste0("min(", context$y_display, ")"),
      paste0("max(", context$y_display, ")"),
      paste0("mean(", context$y_display, ")"),
      paste0("sd(", context$y_display, ")"),
      "n_cp()"
    )
  )
}


warn_legacy_prior_constants = function(prior, context) {
  prior_text = unlist(prior[vapply(prior, is.character, logical(1))], use.names = TRUE)
  if (length(prior_text) == 0)
    return(invisible(NULL))

  map = legacy_prior_map(context)
  pattern = paste0("\\b(?:", paste(map$token, collapse = "|"), ")\\b")
  used = unique(unlist(stringr::str_extract_all(prior_text, pattern)))
  if (length(used) > 0) {
    warning(
      "Deprecated prior data constant(s): ", and_collapse(used), ". Use readable data ",
      "expressions such as `min(", context$x_display, ")`, `sd(", context$y_display,
      ")`, `max(", context$x_display, ") - min(", context$x_display,
      ")`, `segment_width(", context$x_display,
      ")`, or `n_cp()` instead.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


translate_legacy_prior = function(x, context) {
  map = legacy_prior_map(context)
  for (i in seq_len(nrow(map)))
    x = stringr::str_replace_all(x, paste0("\\b", map$token[i], "\\b"), map$rule[i])
  x
}


resolve_legacy_prior_name = function(name, context) {
  map = legacy_prior_map(context)
  match = match(name, map$token)
  if (is.na(match)) NULL else map$value[match]
}


add_legacy_prior_jags_data = function(jags_data, jags_code, context) {
  map = legacy_prior_map(context)
  used = map$token[vapply(map$token, function(token) {
    stringr::str_detect(jags_code, paste0("\\b", token, "\\b"))
  }, logical(1))]
  if (length(used) == 0)
    return(jags_data)

  warning(
    "Deprecated prior data constant(s) in custom `jags_code`: ",
    and_collapse(used), ". Prefer resolved numeric values; see `prior_summary(fit)`.",
    call. = FALSE
  )
  values = stats::setNames(map$value, map$token)
  jags_data[used] = as.list(values[used])
  jags_data
}


legacy_prior_table = function(fit) {
  segments = get_fit_model_tables(fit)$segments
  if (is.null(segments) || is.null(fit$data) || is.null(fit$prior))
    stop("This legacy mcpfit does not contain enough information to reconstruct its priors.")
  context = prior_context(fit$data, segments)
  all_names = names(fit$prior)
  specs = tibble::tibble(
    parameter = all_names,
    code = unname(unlist(fit$prior)),
    description = "Stored prior from a legacy mcpfit",
    source = "legacy"
  )
  compile_prior_specs(specs, all_names, context)
}
