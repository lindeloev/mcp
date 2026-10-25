# ABOUT: Reparameterizations that make JAGS sample faster. Each one samples a
# default prior on a better-mixing scale and derives the parameter from it. The
# model, the priors, and hence the posterior are unchanged. options(mcp_reparam = FALSE)
# writes the plain model. See https://lindeloev.github.io/mcp/articles/implementation.html
# -----------------

# Sampling scales for population-level coefficients
#
# Fills the `node` columns of the prior table (see add_nodes()) where a parameter is
# sampled via an auxiliary node:
#
#   <name>_<suffix>_ ~ dist(location * width + shift, scale * width)
#   <name> = (<name>_<suffix>_ - derive) / width
#
# The auxiliary nodes are values that stay in place when a change point moves:
# - Segment 1 (suffix "mid"): the intercept (and e.g. `stateb` in `1 + state * x`) as the
#   value at the expected midpoint of segment 1 under the default change point prior.
# - Disjoined segments ("end"): the intercept as the value at the segment end.
# - Joined segments: the lowest-degree x-term as the value at the segment end ("end"),
#   and further x-terms, e.g., I(x^2), as their rise over the segment ("rise").
#
# Only default dnorm() and dt() priors without truncation qualify (the segment-1 intercept
# also with truncation or dunif()). AR/MA and segments bounded by group-level change
# points are skipped.
add_reparams = function(prior_table, coefs, cps, context, reparam = TRUE) {
  # Sample `name` via the auxiliary node <name>_<suffix>_ (see the formulas above)
  sampled_as = c(
    mid = "the value at the expected midpoint of segment 1",
    end = "the value at the segment end",
    rise = "the rise over the segment"
  )
  set_aux = function(name, suffix, shift = NA, width = NA, derive = shift) {
    row = match(name, prior_table$parameter)
    node = paste0(name, "_", suffix, "_")
    parsed = prior_table$parsed[[row]]
    if (!is.na(width)) {
      parsed$args[1] = if (parsed$args[1] == "0") "0" else paste0(parsed$args[1], " * ", width)
      parsed$args[2] = paste0(parsed$args[2], " * ", width)
    }
    # Bounded priors are only shifted (see is_reparam_prior())
    if (!is.na(shift))
      parsed = shift_prior(parsed, function(x) if (x == "0") shift else paste0(x, " + ", shift))
    value = if (is.na(derive)) node else paste0(node, " - (", derive, ")")
    if (!is.na(width))
      value = paste0(if (is.na(derive)) value else paste0("(", value, ")"), " / ", width)

    prior_table$node[row] <<- node
    prior_table$node_prior[[row]] <<- parsed
    prior_table$derive[row] <<- value
    prior_table$sampled_as[row] <<- sampled_as[[suffix]]
  }

  # Only default location-scale priors can be shifted and scaled exactly. The segment-1
  # intercept is shifted without change points, so bounded priors (truncation or dunif) work too.
  is_reparam_prior = function(name, segment_1 = FALSE) {
    row = match(name, prior_table$parameter)
    if (!reparam || is.na(row))
      return(FALSE)
    parsed = prior_table$parsed[[row]]
    prior_table$source[row] == "default" && parsed$dist %in% c("dnorm", "dt", if (segment_1) "dunif") && (segment_1 || !any(nzchar(parsed$trunc)))
  }
  # Terms without other variables (a model-matrix column of ones)
  is_pure = vapply(coefs$matrix_data, function(values) all(values == 1, na.rm = TRUE), logical(1))

  # Helpers to write x-terms, e.g., "xE2_2 * (cp_2 - cp_1)^2"
  width = function(k, x_degree) paste0("(cp_", k, " - cp_", k - 1, ")", x_power(x_degree))
  K = nrow(cps) + 1

  # Each dpar separately. AR/MA are not anchored at change points.
  for (dpar in setdiff(unique(coefs$dpar), c("ar", "ma"))) {
    own = coefs[coefs$dpar == dpar, ]
    own_pure = is_pure[coefs$dpar == dpar]
    # x-terms without other variables, in the segment that defines them
    slopes = own[own$x_degree > 0 & own_pure & own$segment == own$definition_segment, ]

    # Segment 1: the intercept and x-term partners (see get_x_partners()) at the expected midpoint
    for (coef in own$name[own$segment == 1 & own$x_degree == 0]) {
      terms = get_x_partners(own, coef)
      if (nrow(terms) > 0 && is_reparam_prior(coef, segment_1 = TRUE)) {
        x_max = max(context$data[[context$par_x]], na.rm = TRUE)
        at = paste0("(", sprintf("%.15g", context$x_min + (x_max - context$x_min) / (2 * K)), ")")
        shift = paste0(x_partners_at(terms, at), collapse = " + ")
        set_aux(coef, "mid", shift = shift, derive = shift)
      }
    }

    # Segments 2+ with x-terms
    for (k in setdiff(unique(slopes$segment), 1)) {
      # Skip segments bounded by group-level change points
      if (any(!is.na(cps$group_col[cps$name %in% paste0("cp_", c(k - 1, k))])))
        next
      terms = slopes[slopes$segment == k, ]
      terms = terms[order(terms$x_degree), ]  # Lowest degree first
      intercept = own$name[own$segment == k & own$par_type == "Intercept"]

      # Disjoined segment: the intercept as the value at the segment end
      if (length(intercept) == 1) {
        if (is_reparam_prior(intercept))
          set_aux(intercept, "end", shift = paste0(terms$name, " * ", width(k, terms$x_degree), collapse = " + "))
        next
      }

      # Joined segment: the first x-term as the value at the segment end, others as their rise
      start = get_value_at_cp(own, own_pure, k)
      for (i in seq_len(nrow(terms))) {
        if (!is_reparam_prior(terms$name[i]))
          next
        if (i == 1 && !is.null(start)) {
          set_aux(terms$name[i], "end", shift = start, width = width(k, terms$x_degree[i]))
        } else {
          set_aux(terms$name[i], "rise", width = width(k, terms$x_degree[i]))
        }
      }
    }
  }
  prior_table
}


# JAGS expression for a dpar's value at cp_(k-1), where joined segment k starts:
# the terms from earlier segments that are still included in segment k, with
# x-terms at their full width. Segment 1 measures x from 0. NULL if it is not one
# number, e.g., when other variables are included.
get_value_at_cp = function(own, own_pure, k) {
  # Terms from earlier segments that are still included in segment k
  included = own$segment < k & (is.na(own$next_segment) | own$next_segment > k)
  if (any(!own_pure[included]))
    return(NULL)
  rows = own[included, ]
  if (nrow(rows) == 0)
    return("0")

  # Intercepts as-is; x-terms at their full width
  width = ifelse(rows$segment == 1, "cp_1", paste0("(cp_", rows$segment, " - cp_", rows$segment - 1, ")"))
  terms = ifelse(rows$x_degree == 0, rows$name, paste0(rows$name, " * ", width, x_power(rows$x_degree)))
  paste(terms, collapse = " + ")
}
