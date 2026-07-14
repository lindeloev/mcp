# ABOUT: These functions "pad" the regression model from get_formula.R
# resulting in a full JAGS model
# -----------------

#' Make JAGS code for Multiple Change Point model
#'
#' @aliases get_jags_code
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param formula_jags String. The formula string returned by `get_formula_jags()`.
#' @param ST Segment table. Returned by `get_segment_table()`.
#' @param ar_order,ma_order NA or positive integer. The GARMA component orders.
#' @param prior_table Resolved prior metadata from `get_prior()`.
#' @param prior_context Data summaries used while resolving priors.
#' @return String. A JAGS model.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
get_jags_code = function(prior, ST, formula_jags, ar_order, ma_order, family, par_x,
                          prior_table = NULL, prior_context = NULL) {
  prior_description = if (is.null(prior_table)) {
    stats::setNames(rep("Prior", length(prior)), names(prior))
  } else {
    stats::setNames(prior_table$description, prior_table$parameter)
  }
  prior_kind_ = if (is.null(prior_table)) {
    NULL
  } else {
    stats::setNames(prior_table$kind, prior_table$parameter)
  }
  prior_context_ = prior_context

  # Begin building JAGS model. `mm` is short for "mcp model".
  # Add fixed variables.
  mm = paste0("model {")

  ####################################
  # DIRICHLET PRIOR ON CHANGE POINTS #
  ####################################
  # Get change point priors and check if they are Dirichlet
  cps = prior[stringr::str_detect(names(prior), "^cp_[1-9]+$")]
  is_dirichlet = stringr::str_detect(cps, "^dirichlet\\([1-9]+\\)$")
  if (any(is_dirichlet)) {
    if (!all(is_dirichlet))
      stop("All or none of the change point priors can be 'dirichlet(N)' and all N > 0.")

    # Build JAGS code. cp_betas is a simplex. cp_i is scaled to the observed range of x.
    mm = paste0(mm, "
  # Scaled Dirichlet prior on change points
  cp_betas ~ ddirch(c(", paste0(stringr::str_extract(cps, "[0-9]+"), collapse = ", "), ", 1))  # Scaled Dirichlet prior on change points")  # OBS: adds an extra 1
    for (i in seq_along(cps)) {
      mm = paste0(mm, "
  cp_", i, " = ", format_prior_number(prior_context_$x_min),
                  " + sum(cp_betas[1:", i, "]) * ", format_prior_number(prior_context_$x_span),
                  "  # Within the observed change-point span")
    }

    # Clean up. Remove any dirichlet priors from the list of priors
    is_dirichlet2 = stringr::str_detect(prior, "^dirichlet\\([1-9]+\\)$")
    prior[is_dirichlet2] = NULL
  }

  ################
  # OTHER PRIORS #
  ################
  # ... also handles non-Dirichlet priors

  # Split up priors into population and varying
  prior_pop = prior[!names(prior) %in% ST$cp_group]
  prior_varying = prior[names(prior) %in% ST$cp_group]

  # Use get_prior_str() to add population-level priors
  mm = paste0(mm, "
  # mcp helper values\n")

  # Helpers for change points:
  mm = paste0(mm, "  cp_0 = ", format_prior_number(prior_context_$x_min), "\n")
  mm = paste0(mm, "  cp_", max(ST$segment), " = ", format_prior_number(prior_context_$x_max), "

  # Priors for population-level effects\n")
  for (i in seq_along(prior_pop)) {
    name = names(prior_pop)[i]
    mm = paste0(mm, get_prior_str(
      prior_pop, i,
      description = prior_description[[name]],
      kind = if (is.null(prior_kind_)) NULL else prior_kind_[[name]]
    ))
  }


  # Use get_prior_str() to add varying priors
  if (length(prior_varying) > 0) {
    mm = paste0(mm, "\n  # Priors for varying effects\n")
    for (i in 1:length(prior_varying)) {
      mm = paste0(mm, get_prior_str(
        prior = prior_varying,
        i = i,
        varying_group = stats::na.omit(ST$cp_group_col[ST$cp_group == names(prior_varying[i])]),
        description = prior_description[[names(prior_varying)[i]]],
        kind = if (is.null(prior_kind_)) NULL else prior_kind_[[names(prior_varying)[i]]]
      ))
    }
  }


  #########
  # GARMA #
  #########
  has_arma = !is.na(ar_order) || !is.na(ma_order)
  if (has_arma)
    mm = paste0(mm, get_arma_jagscode(ar_order, ma_order, par_x))



  ###########
  # FORMULA #
  ###########
  # Transform formula_jags into JAGS format. Insert par_x and varying indices
  for (i in seq_len(max(ST$segment))) {
    formula_jags = gsub(paste0("CP_", i, "_INDEX"), paste0("[", ST$cp_group_col[i], "[i_]]"), formula_jags)
  }

  # Insert formula_jags
  mm = paste0(mm, "
  # Model and likelihood
  for (i_ in 1:length(", par_x, ")) {")

  # Add JAGS code for fitted values and indent it
  mm = paste0(mm, gsub("\n", "\n    ", formula_jags))


  ##############
  # LIKELIHOOD #
  ##############
  # Transform link-scale predictors to distribution-scale parameters. These are
  # deterministic JAGS nodes and are not monitored by default.
  for (dpar in family$dpar_specs$dpar) {
    spec = get_dpar_spec(family, dpar)
    link_code = paste0("link_", dpar, "_[i_]")
    if (dpar == "mu" && has_arma)
      link_code = paste0(link_code, " + resid_arma_[i_]")

    linkinv_str = get_link_str(spec$link, inverse = TRUE)
    response_code = ifelse(
      linkinv_str == "",
      link_code,
      paste0(linkinv_str, "(", link_code, ")")
    )
    if (!is.na(spec$lower))
      response_code = paste0("max(", format(spec$lower, scientific = TRUE), ", ", response_code, ")")

    mm = paste0(mm, "\n    ", dpar, "_[i_] = ", response_code)
  }

  mu_code = "mu_[i_]"

  # Prepare variance code
  has_weights = !all(is.na(ST$weights))
  weights = ifelse(has_weights, yes = paste0(ST$weights[1], "[i_]"), no = "1")

  # Family- and link-dependent likelihood
  mm = paste0(mm, "\n\n    # Likelihood and log-density for family = ", family$family, "()
    ")

  if (family$family == "gaussian") {
    mm = paste0(mm, ST$y[1], "[i_] ~ dnorm(", mu_code, ", ", weights, " / sigma_[i_]^2)  # SD as precision")
  } else if (family$family == "binomial") {
    mm = paste0(mm, ST$y[1], "[i_] ~ dbin(", mu_code, ", ", ST$trials[1], "[i_])")
  } else if (family$family == "bernoulli") {
    mm = paste0(mm, ST$y[1], "[i_] ~ dbern(", mu_code, ")")
  } else if (family$family == "poisson") {
    mm = paste0(mm, ST$y[1], "[i_] ~ dpois(", mu_code, ")")
  } else if (family$family == "negbinomial") {
    mm = paste0(
      mm,
      "nb_prob_[i_] = shape_[i_] / (shape_[i_] + ", mu_code, ")\n    ",
      ST$y[1], "[i_] ~ dnegbin(nb_prob_[i_], shape_[i_])"
    )
  }

  # Compute link-scale residuals for GARMA
  if (has_arma) {
    if (family$family == "binomial") {
      mm = paste0(
        mm,
        "\n    garma_y_[i_] = min(max(", ST$y[1], "[i_], garma_boundary_[i_]), ", ST$trials[1], "[i_] - garma_boundary_[i_]) / ", ST$trials[1], "[i_]",
        "\n    garma_link_y_[i_] = ", family$linkfun_str, "(garma_y_[i_])"
      )
    } else if (family$family %in% c("poisson", "negbinomial")) {
      mm = paste0(
        mm,
        "\n    garma_y_[i_] = max(", ST$y[1], "[i_], garma_boundary_[i_])",
        "\n    garma_link_y_[i_] = ", family$linkfun_str, "(garma_y_[i_])"
      )
    } else {
      mm = paste0(mm, "\n    garma_link_y_[i_] = ", ST$y[1], "[i_]")
    }
    mm = paste0(mm, "\n    resid_abs_[i_] = garma_link_y_[i_] - link_mu_[i_]")
    if (!is.na(ma_order))
      mm = paste0(mm, "\n    resid_ma_[i_] = garma_link_y_[i_] - link_mu_[i_] - resid_arma_[i_]")
  }


  ###############
  # OTHER STUFF #
  ###############

  # Finish up
  mm = paste0(mm, "
  }
}")

  # Return the model string
  mm
}


#' Get JAGS code for a prior
#'
#' @aliases get_prior_str
#' @keywords internal
#' @noRd
#' @inheritParams mcp
#' @param i The index in `prior` to get code for
#' @param varying_group String or NULL. Null indicates a population-
#'   level prior. String indicates a varying-effects prior (one for each group
#'   level).
#' @param description Short comment to include in generated JAGS code.
#' @param kind One of distribution, alias, expression, or constant.
#' @return A string
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @encoding UTF-8
get_prior_str = function(prior, i, varying_group = NULL,
                          description = "Prior", kind = NULL) {
  # Helpers
  value = prior[[i]]
  name = names(prior[i])
  description = gsub("[\r\n]+", " ", description)
  if (is.null(kind))
    kind = prior_kind(value, names(prior))

  # JAGS does not support dinvgamma or dloginvgamma. Write it manually.
  # If inverse_shape ~ Gamma(a, b), then
  # shape = exp(-log(inverse_shape)) ~ InvGamma(a, b).
  if (stringr::str_detect(value, "^dloginvgamma\\(")) {
    if (!is.null(varying_group))
      stop("dloginvgamma() is currently only supported for population-level parameters.")

    call = parse_prior_call(value)
    if (is.null(call) || call$name != "dloginvgamma" || length(call$args) != 2)
      stop("Expected dloginvgamma(shape, scale). Got '", value, "'.")

    inverse_name = paste0(name, "_inverse")
    return(paste0(
      "  ", inverse_name, " ~ dgamma(", call$args[1], ", ", call$args[2], ")  # ", description, "\n",
      "  ", name, " = -log(", inverse_name, ")  # Log of inverse-gamma parameter\n"
    ))
  }

  if (kind == "distribution") {
    # Convert to precision
    value = sd_to_prec(value)

    # ... and this is a population-level effect
    if (is.null(varying_group)) {
      return(paste0("  ", name, " ~ ", value, "  # ", description, "\n"))
    } else {
      # It is a varying effect!
      return(paste0("  for (", varying_group, "_ in 1:n_unique_", varying_group, ") {
    ", name, "_uncentered[", varying_group, "_] ~ ", value, "  # ", description, "
  }
  ", name, " = ", name, "_uncentered - mean(", name, "_uncentered)  # vectorized zero-centering\n"))
    }
  }

  paste0("  ", name, " = ", value, "  # ", description, "\n")
}


#' Transform a JAGS Prior from SD to Precision.
#'
#' JAGS uses precision rather than SD. This function converts
#' `dnorm(4.2, 1.3)` into `dnorm(4.2, 1/1.3^2)`. It allows users to specify
#' priors using SD and then it's transformed for the JAGS code. It works for the
#' following distributions: dnorm|dt|dcauchy|ddexp|dlogis|dlnorm. In all of
#' these,
#' tau/sd is the second parameter.
#'
#' @aliases sd_to_prec
#' @param prior_str String. A JAGS prior. Can be truncated, e.g.
#'   `dt(3, 2, 1) T(my_var, )`.
#' @return A string
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @encoding UTF-8
#' @export
sd_to_prec = function(prior_str) {
  parts = split_prior_truncation(prior_str)
  call = parse_prior_call(parts$distribution)
  if (is.null(call) || call$name %notin% c("dnorm", "dt", "dcauchy", "ddexp", "dlogis", "dlnorm"))
    return(prior_str)
  if (length(call$args) < 2)
    stop("Expected a scale as the second argument of '", prior_str, "'.")

  call$args[2] = paste0("1/(", gsub(" ", "", call$args[2]), ")^2")
  converted = paste0(call$name, "(", paste(call$args, collapse = ", "), ") ")
  if (!is.null(parts$truncation))
    converted = paste0(converted, gsub("\\s+", "", parts$truncation))
  converted
}
