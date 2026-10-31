# ABOUT: Initial values for the sampler: checks of user-supplied `inits` and defaults for the rest.
# ------------------------------------------------------------------------


# Check that `inits` only set parameters that are sampled directly (see add_nodes())
assert_inits = function(inits, pars) {
  invalid = setdiff(unlist(lapply(if (is.list(inits[[1]])) inits else list(inits), names)), pars)
  if (length(invalid) > 0)
    stop("`inits` can only set these parameters of this model: `", paste0(pars, collapse = "`, `"), "`. Got `", paste0(unique(invalid), collapse = "`, `"),
         "`. Parameters with constant priors and parameters that are reparameterized for faster sampling (see `mcp_reparam` in ?mcp-package) cannot be initialized.", call. = FALSE)
}


# Initial values where users supplied none: ordered change points (see get_cp_inits()) and prior medians of
# other population-level parameters, e.g., SDs away from zero. Custom JAGS code may define parameters differently.
get_default_inits = function(prior_table, parameters, x_range, inits, custom_jags_code = NULL) {
  is_cp = parameters$scope == "population" & parameters$role == "change_point"
  cp_rows = prior_table[match(parameters$parameter[is_cp], prior_table$parameter), ]
  is_median = !is.na(prior_table$node) & prior_table$parameter %in% parameters$parameter[parameters$scope == "population" & !is_cp]
  node_priors = stats::setNames(prior_table$node_prior[is_median], prior_table$node[is_median])
  defaults = c(get_cp_inits(cp_rows, x_range, if (!is.list(inits[[1]])) inits), get_median_inits(node_priors))  # Ignores chain-specific inits
  if (is.null(custom_jags_code)) defaults else
    defaults[vapply(names(defaults), function(name) grepl(paste0("\\n\\s*", name, " ~"), custom_jags_code), logical(1))]
}


# Initial values at the prior medians of nodes. Covers priors whose arguments and truncation bounds are numbers.
# - prior: Named list of parsed priors, e.g., `prior_table$node_prior` named by `prior_table$node`.
get_median_inits = function(prior) {
  as_number = function(x) tryCatch(as.numeric(eval(str2lang(x), baseenv())), error = function(e) NA_real_)

  inits = list()
  for (name in names(prior)) {
    parsed = prior[[name]]
    args = vapply(parsed$args, as_number, numeric(1), USE.NAMES = FALSE)
    bounds = c(-Inf, Inf)
    is_bounded = nzchar(parsed$trunc)
    bounds[is_bounded] = vapply(parsed$trunc[is_bounded], as_number, numeric(1), USE.NAMES = FALSE)
    if (!anyNA(args) && !anyNA(bounds))
      inits[[name]] = prior_quantile(parsed$dist, args, bounds, 0.5)
  }
  inits
}


# Initial values for the change points with a prior distribution, except those in `user_inits`:
# equispaced in the observed range (`x_range`) with room for minimum segment lengths, then moved the
# least possible to respect their bounds and order. JAGS fails to initialize if they violate cp_order_.
# - cp_rows: Rows of the prior table for the population-level change points, in order.
get_cp_inits = function(cp_rows, x_range, user_inits = list()) {
  n_cp = nrow(cp_rows)
  if (n_cp == 0)
    return(list())

  # Lower and upper bounds (dunif() or T()) and minimum segment length (see add_nodes()).
  # Constant change points and user-supplied initial values are bounded by themselves.
  env = stats::setNames(as.list(x_range), paste0("cp_", c(0, n_cp + 1)))
  evaluate = function(x, empty) if (is.null(x) || !nzchar(x)) empty else eval(str2lang(x), env, baseenv())
  is_free = !is.na(cp_rows$node)
  lower = upper = min_lengths = numeric(n_cp)
  for (i in seq_len(n_cp)) {
    fixed = c(user_inits[[cp_rows$parameter[i]]], if (cp_rows$kind[i] == "constant") as.numeric(cp_rows$value[i]), NA)[1]
    bounds = if (is_free[i]) get_prior_bounds(cp_rows$node_prior[[i]])
    min_lengths[i] = evaluate(cp_rows$min_length[i], 0)
    lower[i] = if (!is.na(fixed)) fixed else evaluate(bounds[1], -Inf)
    upper[i] = if (!is.na(fixed)) fixed else evaluate(bounds[2], Inf)
  }

  # Keep 1% of the range of x from the bounds and neighboring change points, but within narrow priors
  margin = pmin(0.01 * diff(x_range), (upper - lower) / 3)
  spacing = min_lengths + c(0, pmax(margin[-1], margin[-n_cp]))
  lower = lower + margin
  upper = upper - margin

  # Earliest and latest positions of ordered change points within these bounds
  for (i in seq_len(n_cp)[-1])
    lower[i] = max(lower[i], lower[i - 1] + spacing[i])
  for (i in rev(seq_len(n_cp)[-1]))
    upper[i - 1] = min(upper[i - 1], upper[i] - spacing[i])
  if (any(lower > upper))
    stop("The change point priors (and initial values) leave no room for ordered change points. Check the lower and upper bounds (in dunif() or T()), minimum segment lengths, and initial values of the change points.", call. = FALSE)

  target = x_range[1] + cumsum(min_lengths) + max((diff(x_range) - sum(min_lengths)) / (n_cp + 1), 0) * seq_len(n_cp)
  inits = stats::setNames(as.list(pmin(pmax(target, lower), upper)), cp_rows$parameter)
  inits[is_free & cp_rows$parameter %notin% names(user_inits)]
}
