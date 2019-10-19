##' From the right hand side of a formula for a mixed-effects model,
##' determine the pairs of expressions that are separated by the
##' vertical bar operator.  Also expand the slash operator in grouping
##' factor expressions and expand terms with the double vertical bar operator
##' into separate, independent random effect terms.
##'
##' @title Determine random-effects expressions from a formula
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @param term a mixed-model formula
##' @return pairs of expressions that were separated by vertical bars
##' @section Note: This function is called recursively on individual
##' terms in the model, which is why the argument is called \code{term} and not
##' a name like \code{form}, indicating a formula.
##' @author From package \code{lme4}.
##' @importFrom methods is
##' @examples
##' \dontrun{
##' findbars(f1 <- Reaction ~ Days + (Days|Subject))
##' ## => list( Days | Subject )
##' findbars(y ~ Days + (1|Subject) + (0+Days|Subject))
##' ## => list of length 2:  list ( 1 | Subject ,  0+Days|Subject)
##' findbars(~ 1 + (1|batch/cask))
##' ## => list of length 2:  list ( 1 | cask:batch ,  1 | batch)
##' identical(findbars(~ 1 + (Days || Subject)),
##'     findbars(~ 1 + (1|Subject) + (0+Days|Subject)))
##' stopifnot(identical(findbars(f1),
##'                     list(expression(Days | Subject)[[1]])))
##' }
##' @family utilities
##' @keywords models utilities
## @export
findbars <- function(term)
{
  ## Recursive function applied to individual terms
  fb <- function(term)
  {
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(fb(term[[2]]))
    stopifnot(is.call(term))
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(fb(term[[2]]))
    c(fb(term[[2]]), fb(term[[3]]))
  }
  ## Expand any slashes in the grouping factors returned by fb
  expandSlash <- function(bb)
  {
    ## Create the interaction terms for nested effects
    makeInteraction <- function(x)
    {
      if (length(x) < 2) return(x)
      trm1 <- makeInteraction(x[[1]])
      trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
      list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
    }
    ## Return the list of '/'-separated terms
    slashTerms <- function(x)
    {
      if (!("/" %in% all.names(x))) return(x)
      if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor",call.=FALSE)
      list(slashTerms(x[[2]]), slashTerms(x[[3]]))
    }

    if (!is.list(bb))
      expandSlash(list(bb))
    else
      unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
          ## lapply(unlist(...)) - unlist returns a flattened list
          lapply(unlist(makeInteraction(trms)),
                 function(trm) substitute(foo|bar, list(foo = x[[2]], bar = trm)))
        else x
      }))
  }## {expandSlash}

  modterm <- expandDoubleVerts(
    if(is(term, "formula")) term[[length(term)]] else term)
  expandSlash(fb(modterm))
}





##' Remove the random-effects terms from a mixed-effects formula,
##' thereby producing the fixed-effects formula.
##'
##' @title Omit terms separated by vertical bars in a formula
##' @param term the right-hand side of a mixed-model formula
##' @return the fixed-effects part of the formula
##' @section Note: This function is called recursively on individual
##' terms in the model, which is why the argument is called \code{term} and not
##' a name like \code{form}, indicating a formula.
##' @author From package \code{lme4}
##' @importFrom stats reformulate
##' @examples
##' \dontrun{
##' nobars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days
##' }
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @family utilities
##' @keywords models utilities
## @export
nobars <- function(term) {
  nb <- nobars_(term)  ## call recursive version
  if (is(term,"formula") && length(term)==3 && is.symbol(nb)) {
    ## called with two-sided RE-only formula:
    ##    construct response~1 formula
    nb <- reformulate("1",response=deparse(nb))
  }
  ## called with one-sided RE-only formula, or RHS alone
  if (is.null(nb)) {
    nb <- if (is(term,"formula")) ~1 else 1
  }
  nb
}

nobars_ <- function(term)
{
  if (!anyBars(term)) return(term)
  if (isBar(term)) return(NULL)
  if (isAnyArgBar(term)) return(NULL)
  if (length(term) == 2) {
    nb <- nobars_(term[[2]])
    if(is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars_(term[[2]])
  nb3 <- nobars_(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

isBar <- function(term) {
  if(is.call(term)) {
    if((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
      return(TRUE)
    }
  }
  FALSE
}

isAnyArgBar <- function(term) {
  if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
    for(i in seq_along(term)) {
      if(isBar(term[[i]])) return(TRUE)
    }
  }
  FALSE
}

anyBars <- function(term) {
  any(c('|','||') %in% all.names(term))
}


##' From the right hand side of a formula for a mixed-effects model,
##' expand terms with the double vertical bar operator
##' into separate, independent random effect terms.
##'
##' @title Expand terms with \code{'||'} notation into separate \code{'|'} terms
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @param term a mixed-model formula
##' @return the modified term
##' @family utilities
##' @importFrom stats formula terms as.formula
##' @keywords models utilities
## @export
expandDoubleVerts <- function(term)
{
  expandDoubleVert <- function(term) {
    frml <- formula(substitute(~x,list(x=term[[2]])))
    ## FIXME: do this without paste and deparse if possible!
    ## need term.labels not all.vars to capture interactions too:
    newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
    if(attr(terms(frml), "intercept")!=0)
      newtrms <- c("1", newtrms)

    as.formula(paste("~(",
                     paste(vapply(newtrms, function(trm)
                       paste0(trm, "|", deparse(term[[3]])), ""),
                       collapse=")+("), ")"))[[2]]
  }

  if (!is.name(term) && is.language(term)) {
    if (term[[1]] == as.name("(")) {
      term[[2]] <- expandDoubleVerts(term[[2]])
    }
    stopifnot(is.call(term))
    if (term[[1]] == as.name('||'))
      return( expandDoubleVert(term) )
    ## else :
    term[[2]] <- expandDoubleVerts(term[[2]])
    if (length(term) != 2) {
      if(length(term) == 3)
        term[[3]] <- expandDoubleVerts(term[[3]])
    }
  }
  term
}
