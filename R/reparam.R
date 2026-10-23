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
# The following are sampled as values that stay in place when a change point moves.
# All dpars except ar and ma follow the same rules.
#
# - Segment 1: The intercept as the value at a fixed point, the expected midpoint of
#   segment 1 under the default change point prior, min(x) + (max(x) - min(x)) / (2 * K).
#   Likewise for other coefficients with x-terms, e.g., `stateb` in `1 + state * x`.
# - Joined segments: The lowest-degree x-term as the value at the segment end, and
#   further x-terms (e.g., I(x^2)) as their rise over the segment.
# - Disjoined segments: The intercept as the value at the segment end.
#
# Only default dnorm() and dt() priors without truncation are reparameterized, except
# that the segment-1 intercept may also have a truncated or dunif() prior.
# Segments bounded by group-level change points are skipped because their data use
# group-specific widths.
add_reparams = function(prior_table, predictors, cps, context, reparam = TRUE) {
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
  is_pure = vapply(predictors$matrix_data, function(values) all(values == 1, na.rm = TRUE), logical(1))

  # Helpers to write x-terms, e.g., "xE2_2 * (cp_2 - cp_1)^2"
  power = function(x_factor) ifelse(x_factor == "x", "", sub("x", "", x_factor, fixed = TRUE))
  width = function(k, x_factor) paste0("(cp_", k, " - cp_", k - 1, ")", power(x_factor))
  x_terms_at = function(terms, at) paste0(terms$name, " * ", at, power(terms$x_factor), collapse = " + ")
  K = nrow(cps) + 1

  # Each dpar separately. AR/MA are not anchored at change points.
  for (dpar in setdiff(unique(predictors$dpar), c("ar", "ma"))) {
    own = predictors[predictors$dpar == dpar, ]
    own_pure = is_pure[predictors$dpar == dpar]
    # x-terms without other variables, in the segment that defines them
    slopes = own[own$x_factor != "1" & own_pure & own$segment == own$definition_segment, ]

    # Segment 1: the intercept as the value at a fixed point. Likewise for other coefficients
    # with x-terms, e.g., `stateb` with `stateb.x` in `1 + state * x` (see get_x_partners()).
    for (coef in own$name[own$segment == 1 & own$x_factor == "1"]) {
      terms = get_x_partners(own, coef)
      if (nrow(terms) > 0 && is_reparam_prior(coef, segment_1 = TRUE)) {
        x_max = max(context$data[[context$x_name]], na.rm = TRUE)
        at = paste0("(", sprintf("%.15g", context$x_min + (x_max - context$x_min) / (2 * K)), ")")
        set_aux(coef, "mid", shift = x_terms_at(terms, at), derive = x_terms_at(terms, at))
      }
    }

    # Segments 2+ with x-terms
    for (k in setdiff(unique(slopes$segment), 1)) {
      # Skip segments bounded by group-level change points
      if (any(!is.na(cps$group_col[cps$name %in% paste0("cp_", c(k - 1, k))])))
        next
      terms = slopes[slopes$segment == k, ]
      terms = terms[order(nchar(power(terms$x_factor))), ]  # Lowest degree first
      intercept = own$name[own$segment == k & own$par_type == "Intercept"]

      # Disjoined segment: the intercept as the value at the segment end
      if (length(intercept) == 1) {
        if (is_reparam_prior(intercept))
          set_aux(intercept, "end", shift = paste0(terms$name, " * ", width(k, terms$x_factor), collapse = " + "))
        next
      }

      # Joined segment: the first x-term as the value at the segment end, others as their rise
      start = get_value_at_cp(own, own_pure, k)
      for (i in seq_len(nrow(terms))) {
        if (!is_reparam_prior(terms$name[i]))
          next
        if (i == 1 && !is.null(start)) {
          set_aux(terms$name[i], "end", shift = start, width = width(k, terms$x_factor[i]))
        } else {
          set_aux(terms$name[i], "rise", width = width(k, terms$x_factor[i]))
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
  power = ifelse(rows$x_factor %in% c("1", "x"), "", sub("x", "", rows$x_factor, fixed = TRUE))
  terms = ifelse(rows$x_factor == "1", rows$name, paste0(rows$name, " * ", width, power))
  paste(terms, collapse = " + ")
}
