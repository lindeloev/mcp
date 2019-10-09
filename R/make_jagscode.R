#' Make JAGS code for Multiple Change Point model
#'
#' @param data A data.frame or tibble containing the variables expressed in `segments` in long format.
#' @param segments A vector of formulas - one for each segment. The left-hand side specifies the form of the change point (on x). The right-hand side specifices the form of intercepts and slopes. See `mcp` examples for details.
#' @param prior A named list with parameter names as names and a JAGS distribution as value, e.g., `list(int_1 = "dunif(10, 30)")`.
#' @param param_x A string. Only relevant if no segments contains slope (no hint at what x is). Set this, e.g., param_x = "time".
#' @keywords jags, model
#' @export
#' @examples
#'
make_jagscode = function(data, segments, prior, param_x = NULL) {

  # General handy values
  nsegments = length(segments)
  if(nsegments == 0) stop("At least one segment is needed")

  all_pars = list()
  for(i in 1:nsegments) all_pars[paste0("segment_", i)] = list()

  # Figure out the name of the y-axis an the x-axis
  param_name_y = as.character(segments[[1]][[2]])

  all_slopes = c()
  for(i in 1:nsegments) {
    slope_name = attr(terms(segments[[i]]), "term.labels")
    if(!is.null(slope_name)) all_slopes = c(all_slopes, slope_name)
  }
  if(!length(all_slopes)) {
    if(!is.null(param_x)) param_name_x = param_x
    else stop("This is a plateau-only model, but `param_x` is missing: mcp(..., param_x = \"name\")")
  } else param_name_x = all_slopes[1]

  # Check that only one slope has been used
  if(length(unique(all_slopes)) > 1) stop(paste0("More than one slope predictor in the segments: ", paste0(all_slopes, collapse = ", ")))


  # Begin model. `mm` is short for "mcp model".
  # Make it global so that mm_add can manipulate it.
  mm <<- ""

  # Add fixed variables.
  mm_add("
  data {
    # X values
    MINX = ", min(data[, param_name_x]), "
    MAXX = ", max(data[, param_name_x]), "
    MEANX = ", mean(unlist(data[, param_name_x])), "
    SDX = ", sd(unlist(data[, param_name_x])), "

    # Y values
    MINY = ", min(data[, param_name_y]), "
    MAXY = ", max(data[, param_name_y]), "
    MEANY = ", mean(unlist(data[, param_name_y])), "
    SDY = ", sd(unlist(data[, param_name_y])), "
  }

  model {"
  )


  # Add prior on sigma
  if("sigma" %in% names(prior)) {
    par_value = prior[["sigma"]]
  } else {
    par_value = "dnorm(0, 1 / SDY^2) T(0, )  # mcp default"
  }
  mm_add("\n
    # Prior on sigma (standard deviation of likelihood function)
    sigma ~ ", par_value)
  all_pars["sigma"] = "sigma"


  # Add priors on change points
  mm_add("\n\n  # Priors on change points. Enforces order.\n")
  mm_add("  cp_0 = -10^100  # mcp helper value; minus infinity\n")
  if(nsegments == 1) cp_loop = c()  # No change point for 1-segment models
  else cp_loop = 1:(nsegments - 1)
  for(i in cp_loop) {
    param_name = paste0("cp_", i)
    all_pars[param_name] = param_name

    # If prior is manually specified
    if(param_name %in% names(prior)) {
      if(i > 1) trunc = paste0(" T(cp_", i-1, ", )")
      else trunc = ""
      par_value = paste0(prior[[param_name]], trunc)
    }

    # Otherwise use default priors
    else {
      if(i == 1) par_value = "dunif(MINX, MAXX)"
      else if(i == nsegments) par_value = paste0("dunif(cp_", i-1, ", MAXX)")
      else par_value = paste0("dunif(cp_", i-1, ", MAXX)")
      par_value = paste0(par_value, "  # mcp default")  # Comment
    }

    # Add new line to mm
    new_prior = paste0("  ", param_name, " ~ ", par_value, "\n")
    mm_add(new_prior)
  }
  mm_add("  cp_", i+1, " = 10^100  # mcp helper value; plus infinity\n")


  # Population-priors on intercepts
  mm_add("\n  # Priors on intercepts\n")

  for(i in 1:nsegments) {
    # If there is an intercept
    if(attr(terms(segments[[i]]), "intercept")) {
      param_name = paste0("int_", i)
      all_pars[[paste0("segment_", i)]]["int_free"] = param_name

      # If the prior is manually specified
      if(param_name %in% names(prior)) {
        par_value = prior[[param_name]]
      }

      # Otherwise use default priors
      else {
        par_value = "dt(0, 1 / (SDY * 3)^2, 1)  # mcp default"
      }

      # Add it to the model
      new_prior = paste0("  ", param_name, " ~ ", par_value, "\n")
      mm_add(new_prior)
    }
  }


  # Population priors on slopes
  mm_add("\n  # Priors on slopes\n")

  for(i in 1:nsegments) {
    # If there is a slope
    if(length(attr(terms(segments[[i]]), "term.labels"))) {

      # Set parameter name
      param_name = paste0(param_name_x, "_", i)
      all_pars[[paste0("segment_", i)]]["slope_free"] = param_name

      # If there is an explicit prior...
      if(param_name %in% names(prior)) {
        par_value = prior[[param_name]]
      }

      # Otherwise use default priors
      else {
        par_value = "dt(0, 1 / (SDY / (MAXX - MINX) * 3)^2, 1)  # mcp default"
      }

      # Add it to the model
      new_prior = paste0("  ", param_name, " ~ ", par_value, "\n")
      mm_add(new_prior)
    }
  }



  # TO DO
  # * Build list of parameters involved in each segment
  mm_add("
    # Model and likelihood
    for(i_ in 1:length(", param_name_x, ")) {

      # Additive segments
      y_[i_] = \n")

  for(i in 1:nsegments) {
    segment = all_pars[[paste0("segment_", i)]]
    segment_pars = names(segment)

    if(!is.null(segment_pars)) {  # If there are parameters for this segment
      # Indicator whether this segment is included
      mm_add("             (", param_name_x, "[i_] > cp_", i - 1, ") * (")

      # Add intercept of this segment
      has_int = "int_free" %in% segment_pars
      has_slope = "slope_free" %in% segment_pars
      if(has_int) {
        mm_add("int_", i)
      }

      # Add int and slope if both are in this segment
      if(has_int & has_slope) mm_add(" + ")

      # Add slope of this segment
      if(has_slope) {
        mm_add(param_name_x, "_", i, " * (min(", param_name_x, "[i_], cp_", i, ") ")
        if(i > 1) mm_add("- cp_", i - 1)
        mm_add(")")
      }

      # Close line
      if(i < nsegments) mm_add(") +  #segment ", i, "\n")  # More terms to come
      else mm_add(")  # segment ", i)  # This was the last
    }
    else {
      if(i < nsegments) mm_add("             0 +  #segment ", i, "\n")  # More terms to come
      else mm_add("             0  # segment ", i)  # This was the last
    }
  }

  # Finally the likelihood
  mm_add("\n\n    # Likelihood
      ", param_name_y, "[i_] ~ dnorm(y_[i_], 1 / sigma^2)")

  # Log-density.
  mm_add("\n\n    # Log-density for WAIC computation
      loglik_[i_] = logdensity.norm(", param_name_y, "[i_], y_[i_], 1 / sigma^2)")

  # Finish up
  mm_add("
    }
  }")

  return(list(
    model = mm,
    all_pars = all_pars,
    param_name_x = param_name_x,
    param_name_y = param_name_y
  ))
}




#' Quick paste
#'
#' @param ... Terms to be paste0'ed to the string `mm`.
#' @export
#' @examples
#' mm = "model{\n"
#' mm_add('  sigma = dnorm(', mu, ', ' 1 /', sigma, '^2)\n')
#' mm_add('\}')
#' cat(mm)
#'
mm_add = function(...) {
  mm <<- paste0(mm, ...)
}
