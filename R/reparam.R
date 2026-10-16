# ABOUT: Reparameterizations that make JAGS sample faster. Each one samples a
# default prior on a better-mixing scale and derives the parameter from it. The
# model, the priors, and hence the posterior are unchanged. options(mcp_reparam = FALSE)
# writes the plain model. See https://lindeloev.github.io/mcp/articles/implementation.html
# -----------------

# Sampling scales for population-level coefficients
#
# Adds the columns `aux_suffix`, `aux_shift`, `aux_width`, and `aux_derive` to the
# prior table. Where `aux_suffix` is not NA, the parameter is sampled via an
# auxiliary node (see get_aux_prior_str()):
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
add_reparams = function(table, predictors, cps, context, reparam = TRUE) {
  # Empty columns, filled by set_aux() for reparameterized parameters
  aux_cols = c("aux_suffix", "aux_shift", "aux_width", "aux_derive")
  table[aux_cols] = NA_character_
  set_aux = function(name, suffix, shift = NA, width = NA, derive = shift) {
    table[match(name, table$parameter), aux_cols] <<- list(suffix, shift, width, derive)
  }

  # Only default location-scale priors can be shifted and scaled exactly. The segment-1
  # intercept is shifted without change points, so bounded priors (truncation or dunif) work too.
  is_reparam_prior = function(name, segment_1 = FALSE) {
    row = match(name, table$parameter)
    call = parse_prior_call(split_prior_truncation(table$value[row])$distribution)
    reparam && !is.na(row) && table$source[row] == "default" && !is.null(call) &&
      call$name %in% c("dnorm", "dt", if (segment_1) "dunif") && (segment_1 || !grepl("T(", table$value[row], fixed = TRUE))
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
  table
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


# JAGS code for a parameter sampled via an auxiliary node (see add_reparams()).
# `aux` is the parameter's row in the prior table.
get_aux_prior_str = function(name, value, description, aux) {
  # Prior on the auxiliary node: dist(location * width + shift, scale * width)
  aux_name = paste0(name, "_", aux$aux_suffix, "_")
  if (!is.na(aux$aux_shift) || !is.na(aux$aux_width)) {
    parts = split_prior_truncation(value)
    call = parse_prior_call(parts$distribution)
    location = call$args[1]
    if (!is.na(aux$aux_width)) {
      location = if (location == "0") "0" else paste0(location, " * ", aux$aux_width)
      call$args[2] = paste0(call$args[2], " * ", aux$aux_width)
    }
    if (!is.na(aux$aux_shift))
      location = if (location == "0") aux$aux_shift else paste0(location, " + ", aux$aux_shift)
    call$args[1] = location

    # Bounded priors are only shifted (see is_reparam_prior())
    shift = function(bound) if (bound == "") "" else paste0(bound, " + ", aux$aux_shift)
    if (call$name == "dunif")
      call$args[2] = shift(call$args[2])
    value = paste0(call$name, "(", paste(call$args, collapse = ", "), ")")
    if (!is.null(parts$truncation))
      value = paste0(value, " T(", paste(vapply(parse_prior_call(parts$truncation)$args, shift, character(1)), collapse = ", "), ")")
  }

  # The parameter derived from the auxiliary node: (aux - derive) / width
  derive = if (is.na(aux$aux_derive)) aux_name else paste0(aux_name, " - (", aux$aux_derive, ")")
  if (!is.na(aux$aux_width))
    derive = paste0(if (is.na(aux$aux_derive)) derive else paste0("(", derive, ")"), " / ", aux$aux_width)
  # JAGS comment and code
  sampled_as = c(
    mid = "; sampled as the value at the expected midpoint of segment 1",
    end = "; sampled as the value at the segment end",
    rise = "; sampled as the rise over the segment"
  )
  paste0(
    "  ", aux_name, " ~ ", prior_to_jags(value), "  # ", description, sampled_as[[aux$aux_suffix]], "\n",
    "  ", name, " = ", derive, "\n"
  )
}
