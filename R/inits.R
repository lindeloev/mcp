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
get_default_inits = function(prior_table, parameters, x_range, inits, custom_jags_code = NULL, cp_deviance = NULL) {
  is_cp = parameters$scope == "population" & parameters$role == "change_point"
  cp_rows = prior_table[match(parameters$parameter[is_cp], prior_table$parameter), ]
  is_median = !is.na(prior_table$node) & prior_table$parameter %in% parameters$parameter[parameters$scope == "population" & !is_cp]
  node_priors = stats::setNames(prior_table$node_prior[is_median], prior_table$node[is_median])
  defaults = c(get_cp_inits(cp_rows, x_range, if (!is.list(inits[[1]])) inits, cp_deviance), get_median_inits(node_priors))  # Ignores chain-specific inits
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
# With `cp_deviance` (see get_cp_deviance()), they are then moved to where the mu model fits best (see scan_cps()),
# within the central 98% of priors with numeric arguments.
# - cp_rows: Rows of the prior table for the population-level change points, in order.
get_cp_inits = function(cp_rows, x_range, user_inits = list(), cp_deviance = NULL) {
  n_cp = nrow(cp_rows)
  if (n_cp == 0)
    return(list())

  # Lower and upper bounds (dunif() or T()) and minimum segment length (see add_nodes()).
  # Constant change points and user-supplied initial values are bounded by themselves.
  env = stats::setNames(as.list(x_range), paste0("cp_", c(0, n_cp + 1)))
  evaluate = function(x, empty) if (is.null(x) || !nzchar(x)) empty else eval(str2lang(x), env, baseenv())
  is_free = !is.na(cp_rows$node)
  lower = upper = min_lengths = numeric(n_cp)
  central = matrix(c(-Inf, Inf), 2, n_cp)
  for (i in seq_len(n_cp)) {
    fixed = c(user_inits[[cp_rows$parameter[i]]], if (cp_rows$kind[i] == "constant") as.numeric(cp_rows$value[i]), NA)[1]
    bounds = if (is_free[i]) get_prior_bounds(cp_rows$node_prior[[i]])
    min_lengths[i] = evaluate(cp_rows$min_length[i], 0)
    lower[i] = if (!is.na(fixed)) fixed else evaluate(bounds[1], -Inf)
    upper[i] = if (!is.na(fixed)) fixed else evaluate(bounds[2], Inf)
    args = if (is.na(fixed) && !is.null(cp_deviance)) tryCatch(vapply(cp_rows$node_prior[[i]]$args, evaluate, numeric(1), empty = NA), error = function(e) NA) else NA  # NA for conditional priors
    quantiles = if (!anyNA(args)) prior_quantile(cp_rows$node_prior[[i]]$dist, args, c(lower[i], upper[i]), c(0.01, 0.99))
    if (!is.null(quantiles))
      central[, i] = quantiles
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
  cp = pmin(pmax(target, lower), upper)
  is_init = is_free & cp_rows$parameter %notin% names(user_inits)
  if (!is.null(cp_deviance)) {
    # Scan within the central priors if the change points still fit there in order
    central_cp = pmin(pmax(cp, central[1, ]), central[2, ])
    if (all(pmax(lower, central[1, ]) <= pmin(upper, central[2, ])) && all(diff(central_cp) >= spacing[-1])) {
      cp = central_cp
      lower = pmax(lower, central[1, ])
      upper = pmin(upper, central[2, ])
    }
    cp = scan_cps(cp, is_init, lower, upper, spacing, cp_deviance)
  }
  stats::setNames(as.list(cp), cp_rows$parameter)[is_init]
}


# Move change points to where the mu model fits best: first anywhere (then re-sorted), which can
# reorder them, then between their neighbors. Moves are kept if they lower the deviance and respect
# the bounds and spacing. Ties keep the current value, so change points that do not affect the deviance stay.
# - cp: Change points within their bounds, e.g., equispaced.
# - is_init: Which change points to move. The others are fixed.
# - lower, upper, spacing: Bounds and minimum distance to the previous change point, from get_cp_inits().
# - cp_deviance: A function of the change points, from get_cp_deviance().
scan_cps = function(cp, is_init, lower, upper, spacing, cp_deviance, n_grid = 40) {
  best = cp_deviance(cp)
  try_grid = function(k, from, to) {
    for (value in from + (to - from) * (seq_len(n_grid) - 0.5) / n_grid) {
      candidate = cp
      candidate[k] = value
      candidate = sort(candidate)
      if (all(candidate >= lower & candidate <= upper) && all(diff(candidate) >= spacing[-1])) {
        deviance = cp_deviance(candidate)
        if (deviance < best - 1e-9 * abs(best)) {
          best <<- deviance
          cp <<- candidate
        }
      }
    }
  }

  if (is.finite(best)) {
    for (sweep in 1:3) {  # Anywhere, until no change point moves
      previous = cp
      for (k in which(is_init))
        try_grid(k, min(lower[is_init]), max(upper[is_init]))
      if (identical(previous, cp))
        break
    }
    for (k in which(is_init))  # Between neighbors
      try_grid(k, max(lower[k], c(-Inf, cp)[k]), min(upper[k], c(cp, Inf)[k + 1]))
  }
  cp
}


# The deviance of the mu model at given change points, with all population-level mu coefficients fitted by
# glm.fit() with `family$glm_family` and other parameters at zero, e.g., group-level deviations and other dpars. NULL
# where the mu model cannot place all change points (AR/MA terms or segments without mu terms), for likelihood weights,
# and for families without a `glm_family`.
# - max_rows: Evaluate at most this many rows, evenly spread over the sorted x.
get_cp_deviance = function(model_tables, family, data, max_rows = 2000) {
  terms = get_formula_terms(model_tables)
  cps = model_tables$cps
  data_columns = model_tables$data_columns
  is_mu = terms$dpar == "mu" & is.na(terms$group_col) & is.na(terms$offset)
  if (nrow(cps) == 0 || is.null(family$glm_family) || !is.null(data_columns$weights) || any(model_tables$coefs$dpar %in% c("ar", "ma")) || !all(seq_len(nrow(cps) + 1) %in% terms$segment[is_mu]))
    return(NULL)

  # Response on a subset of rows. Binomial responses are proportions weighted by trials.
  y = data[[data_columns$response]]
  trials = if (is.null(data_columns$trials)) rep(1, length(y)) else data[[data_columns$trials]]
  x = data[[data_columns$par_x]]
  rows = which(!is.na(y))
  rows = rows[order(x[rows])][unique(round(seq(1, length(rows), length.out = min(length(rows), max_rows))))]
  designs = evaluate_designs(model_tables, data)
  model_data = list(x = x[rows], matrix = designs$matrix[rows, , drop = FALSE], offsets = lapply(designs$offsets, function(offset) offset[rows]))

  # In draw j, coefficient j is 1 and the others are 0, so link_mu[, j] - link_mu[, p + 1] is column j of the design matrix
  coefs = unique(terms$name[is_mu])
  p = length(coefs)
  zeros = setdiff(stats::na.omit(c(terms$name, cps$group_name)), c(coefs, cps$name))
  pars = c(stats::setNames(lapply(seq_len(p), function(j) as.numeric(seq_len(p + 1) == j)), coefs),
           stats::setNames(rep(list(numeric(p + 1)), length(zeros)), zeros))
  function(cp) {
    pars[cps$name] = lapply(cp, rep, p + 1)
    link_mu = evaluate_formula(model_tables, family, model_data, pars, p + 1, terms)$mu
    fit = tryCatch(suppressWarnings(stats::glm.fit(
      link_mu[, seq_len(p), drop = FALSE] - link_mu[, p + 1], y[rows] / trials[rows], trials[rows], offset = link_mu[, p + 1], family = family$glm_family
    )), error = function(e) NULL)
    if (is.null(fit) || !fit$converged) Inf else fit$deviance
  }
}
