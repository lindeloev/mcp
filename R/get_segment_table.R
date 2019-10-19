# TO DO:
# * Nest cp varying effects too?

source("R/lme4_utils.R")


# Takes a formula and returns a string representation of y, cp, and rhs
unpack_tildes = function(segment, i) {
  if(attributes(terms(segment))$response == 0)
    stop("Empty left-hand side in segment ", i, ": ", segment)

  # Get it as one string
  form_str = as.character(segment)
  form_str = paste(form_str[2], form_str[1], form_str[3])

  # Check for rel(0)
  if("rel(0)" %in% attributes(terms(segment))$term.labels)
    stop("rel(0) is not (currently) supported in segment formulas.")

  # List of strings for each section
  chunks = stringr::str_trim(strsplit(form_str, "~")[[1]])

  if(length(chunks) == 2) {
    # Only one tild. This is the first segment or y is implicit from earlier segment(s)
    return(tibble::tibble(
      form = form_str,
      form_y = ifelse(i == 1, chunks[1], NA),
      #form_cp_str = ifelse(i == 1, NA, chunks[1]),
      form_cp = ifelse(i == 1, NA, paste0(" ~ ", chunks[1])),
      #form_cp = as.formula(paste0(" ~ ", chunks[2])),
      form_rhs = paste0(" ~ ", chunks[2])
      #form_rhs = as.formula(paste0(" ~ ", chunks[2]))
    ))
  } else if(length(chunks) == 3) {
    return(tibble::tibble(
      form = form_str,
      form_y = chunks[1],
      form_cp = paste0(" ~ ", chunks[2]),
      #form_cp = as.formula(paste0(" ~ ", chunks[2])),
      form_rhs = paste0(" ~ ", chunks[3])
      #form_rhs = as.formula(paste0(" ~ ", chunks[3]))
    ))
  } else {
    stop("Error in segment ", i, ": None or more than two tildes in a segment formula.")
  }
}

# Checks if all terms are in the data
check_terms_in_data = function(form, data, i) {
  form_terms = all.vars(form)
  found = form_terms %in% colnames(data)
  if(!all(found))
    stop("Error in segment ", i, ": Term '", paste0(form_terms[!found], collapse="' and '"), "' found in formula but not in data.")
}


# Unpacks y variable name
unpack_y = function(form_y, i) {
  if(!is.na(form_y)) {
    if(!grepl("^[A-Za-z0-9]+$", form_y))
      stop("Error in segment ", i, ": Invalid format for response variable. Only a single column name is (currently) allowed")
  }
  tibble::tibble(y = form_y)
}


#' Takes a cp formula (as astring) and returns its properties
#' @param form_cp Segment formula as string.
#' @param i Positive integer. Segment number.
unpack_cp = function(form_cp, i) {
  if(is.na(form_cp)) {
    return(tibble::tibble(
      cp_int = FALSE,
      cp_int_rel = FALSE,
      cp_ran_int = FALSE,
      cp_ran_group = NA
    ))
  }
  form_cp = as.formula(form_cp)

  # Varying effects
  varying = findbars(form_cp)
  if(!is.null(varying)) {
    if(length(varying) > 1)
      stop("Error in segment", i, " (change point): only one varying effect allowed. Found ", form_cp)

    varying_parts = as.character(varying[[1]])
    if(!varying_parts[2] %in% "1")
      stop("Error in segment ", i, " (change point): Only plain intercepts are allowed in varying effects, e.g., (1|id).", i)

    if(!grepl("^[A-Za-z0-9]+$", varying_parts[2]))
      stop("Error in segment ", i, " (change point): Grouping variable in varying effects.")
  }

  # Fixed effects
  population = attributes(terms(nobars(form_cp)))
  is_int_rel = population$term.labels == "rel(1)"
  if(any(is_int_rel))
    population$term.labels = population$term.labels[-is_int_rel]  # code as no term
  #if(any(is_int_rel))
  #  population$intercept = 1  # There was an intercept. I.e. "~ 0 + rel(1)" is still taken to mean a (relative) intercept.

  if(length(population$term.labels) > 0)
    stop("Error in segment ", i, " (change point): Only intercepts (1 or rel(1)) are allowed in population-level effects.")

  if(is.null(varying) & population$intercept == 0)
    stop("Error in segment ", i, " (change point): no intercept or varying effect. You can do e.g., ~ 1 or ~ (1 |id).")

  # Return as list.
  if(!is.null(varying)) {
    # If there is a varying effect
    return(tibble::tibble(
      cp_int = population$intercept == 1,
      cp_int_rel = any(is_int_rel),  # the intercept is relative
      cp_ran_int = ifelse(varying_parts[2] == "1", TRUE, NA),  # placeholder for later
      cp_ran_group = varying_parts[3]
    ))
  } else {
    # If there is no varying effect
    return(tibble::tibble(
      cp_int = population$intercept == 1,
      cp_int_rel = any(is_int_rel),
      cp_ran_int = FALSE,
      cp_ran_group = NA
    ))
  }
}

unpack_rhs = function(form_rhs, i) {
  form_rhs = formula(form_rhs)

  # Varying effects
  varying = findbars(form_rhs)
  if(!is.null(varying)) {
    # For each varying effect...
    V = tibble::tibble()
    for(term in varying) {
      #varying = lapply(varying, function(x) {
      parts = as.character(term)

      # Check that there is just one grouping term
      if(!grepl("^[A-Za-z0-9]+$", parts[3]))
        stop("Error in segment ", i, " (linear): Grouping variable in varying effects for change points.")

      # Check that nothing is relative
      if(any(stringr::str_detect(parts, "rel\\(")))
        stop("Error in segment ", i, " (linear): rel() not supported in varying effects.")

      # LHS: Split intercepts and variable
      preds = strsplit(gsub(" ", "", parts[2]), "\\+")[[1]]
      slope = preds[!preds %in% c("0", "1")]
      if(length(slope) > 1)
        stop("Error in segment ", i, " (linear): More than one variable in LHS of varying effect.")
      else if(length(slope) == 0)
        # If not slope
        slope = NA

      # Return
      new_row = tibble::tibble(
        int = !"0" %in% preds | parts[2] == "",  # bool. Is integer present?
        slope = slope,
        group = parts[3])
      V = dplyr::bind_rows(V, new_row)
    }
  } else V = NA

  # Population-level intercepts
  population = attributes(terms(nobars(form_rhs)))
  int_rel = population$term.labels == "rel(1)"
  population$term.labels = population$term.labels[!population$term.labels %in% "rel(1)"]  # code as no term
  if(any(int_rel))
    population$intercept = 1


  # Population-level slopes
  n_slopes = length(population$term.labels)
  if(n_slopes > 1) {
    stop("Error in segment ", i, " (linear): Only one slope allowed in population-level effects but '", paste0(population$term.labels, collapse="' and '"), "' found.")
  } else if(n_slopes == 1) {
    # One slope. Code as slope_rel = TRUE/FALSE and x = str
    slope_rel = stringr::str_detect(population$term.labels, "rel\\(")
    population$term.labels = gsub("rel\\(|\\)", "", population$term.labels)
    slope = population$term.labels
  } else if(n_slopes == 0) {
    slope = NA
    slope_rel = FALSE
  }

  # Return as list.
  return(tibble::tibble(
    int = population$intercept == 1,
    int_rel = any(int_rel),
    slope = slope,
    slope_rel = slope_rel,
    varying = list(V)
  ))
}


#' Build a table describing a list of segments
#'
#' Used internally for most mcp functions.
#'
#' @aliases get_segment_table
#' @inheritParams mcp
#' @importFrom stats terms
#' @importFrom utils modifyList
#' @import dplyr
#' @return A tibble.
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @export
#' @examples
#' segments = list(
#'   y ~ 1 + x,
#'   1 ~ 1
#' )
#' data = data.frame(y = 1:5, x=1:5)
#' get_segment_table(segments, NULL, NULL)

get_segment_table = function(segments, data = NULL, param_x = NULL) {
  #############################################################
  # PART 1: BUILD "SEGMENT TABLE (ST)" FROM ISOLATED SEGMENTS #
  #############################################################
  ST = tibble::tibble()
  for(i in 1:length(segments)) {
    # Get ready...
    segment = segments[[i]]
    row = tibble::tibble(segment = i)

    if(!is.null(data))
      check_terms_in_data(segment, data, i)

    # Go! Unpack this segment
    row = dplyr::bind_cols(row, unpack_tildes(segment, i))
    row = dplyr::bind_cols(row, unpack_y(row$form_y, i))
    row = dplyr::bind_cols(row, unpack_cp(row$form_cp, i))
    row = dplyr::bind_cols(row, unpack_rhs(row$form_rhs, i))

    ST = dplyr::bind_rows(ST, row)
  }

  # Return the cols we need
  ST = dplyr::select(ST, -form_y, -form_cp, -form_rhs)


  ################################
  # PART 2: LOOK ACROSS SEGMENTS #
  ################################

  # Check segment 1: rel() not possible here.
  if(any(ST[1, c("cp_int", "cp_int_rel", "cp_ran_int", "cp_ran_group")] != FALSE, na.rm = T))
    stop("Change point defined in first segment. This should not be possible. Submit bug report in the GitHub repo.")
  if(any(ST[1, c("int_rel", "slope_rel")] != FALSE, na.rm = T))
    stop("rel() cannot be used in segment 1. There is nothing to be relative to.")

  # rel() in segment 2+
  rel_slope_after_plateau = lag(is.na(ST$slope), 1) & ST$slope_rel != 0
  if(any(rel_slope_after_plateau))
    stop("rel(slope) is not meaningful after a plateau segment (without a slope). Use absolute slope to get the same behavior. Found in segment ", which(rel_slope_after_plateau))
  if(nrow(ST) > 1) {
    if(ST$cp_int_rel[2] == TRUE)
      stop("rel() cannot be used for change points in segment 2. There are no earlier change points to be relative to. Relative changepoints work from segment 3 and on.")
  }


  # Set ST$x (what is the x-axis dimension?)
  derived_x = unique(na.omit(ST$slope))
  if(length(derived_x) == 1) {
    # One x derived from segments
    if(is.null(param_x))  # param_x not provided. Rely on derived
      ST$x = derived_x
    else if(param_x == derived_x)  # param_x provided and matches devided.
      ST$x = derived_x
    else  # no x info or contradicting info
      stop("param_x provided but it does not match the predictor found in segment RHS.")
  } else if(length(derived_x) == 0) {
    # Zero x derived from segments. Rely on param_x?
    if(all(is.na(ST$slope) & is.character(param_x)))
      ST$x = param_x
    else
      stop("This is a plateau-only model so no x-axis variable could be derived from the segment formulas. Use argument 'param_x' to set it explicitly")
  } else if(length(derived_x) > 1)
    # More than one...
    stop("More than one predictor found: '", paste0(unique(na.omit(ST$slope)), collapse = "' and '"), "'")

  # Only one response variable allowed
  derived_y = unique(na.omit(ST$y))
  if(length(derived_y) != 1)
    stop("There should be exactly one response variable. Found '", paste0(derived_y, collapse="' and '", "'."))

  # Recode relative columns so 0 = not relative. N > 0 is the number of consecutive "relatives"
  ST = ST %>% dplyr::mutate(
    cp_int_rel = sequence(rle(cp_int_rel)$lengths) * cp_int_rel,
    int_rel = sequence(rle(int_rel)$lengths) * int_rel,
    slope_rel = sequence(rle(slope_rel)$lengths) * slope_rel
  ) %>%
    # Replace NAs with the latest y value
    tidyr::fill(y) %>%

    # Add variable names
    mutate(
      int_name = ifelse(int, yes = paste0("int_", segment), no = NA),
      slope_name = ifelse(!is.na(slope), yes = paste0(slope, "_", segment), no = NA),
      slope_code = slope_name,  # Will be modified in next step
      cp_name = paste0("cp_", segment - 1),
      cp_code = cp_name  # Will be modified in next step
    )

  # Add column "cp_code" to do the relative codings, i.e. cp_3 = (cp_2 + cp_3)
  # Add parentheses in the end so it doesn't go crazy.
  # TO DO: there is probably a prettier way to do this.
  for(i in 1:nrow(ST)) {
    ST$cp_code[i] = ifelse(ST$cp_int_rel[i] != 0,
                                 yes = paste0(ST$cp_code[i], "+", ST$cp_code[i-1]),
                                 no = ST$cp_code[i])
  }
  ST$cp_code = ifelse(ST$cp_int_rel != 0,
                            yes = paste0("(", ST$cp_code, ")"),
                            no = ST$cp_code)

  # Do the same for slopes
  for(i in 1:nrow(ST)) {
    ST$slope_code[i] = ifelse(ST$slope_rel[i] != 0,
                                    yes = paste0(ST$slope_code[i], "+", ST$slope_code[i-1]),
                                    no = ST$slope_code[i])
  }
  ST$slope_code = ifelse(ST$slope_rel != 0,
                               yes = paste0("(", ST$slope_code, ")"),
                               no = ST$slope_code)

  # Return
  ST
}
