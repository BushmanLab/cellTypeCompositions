##' Combine duplicate rows of a table of counts
##'
##' When a table of counts has many duplicate rows, computations are
##' sometimes faster if duplicate rows are combined and a vector of
##' case weights is created.
##' @title uniTab Combine Duplicate Rows
##' @param tab matrix or table of counts.  Rows correspond to lineages
##'     and columns to cell types.
##' @return \code{list} with elements \code{tab} - the rows of
##'     \code{unique(tab)}, \code{n} - a vector of row counts, and
##'     \code{row.index} - a mapping of the rows to the \code{tab}
##'     argument to those of \code{uniTab(tab)[["tab"]]}.
##' @export
##' @author Charles Berry
uniTab <- function(tab){
    tab <- unclass(as.matrix(tab))
    utab <- unique(tab)
    tab.index <- match(
        do.call(paste,as.data.frame(tab)),
        do.call(paste,as.data.frame(utab)))
    list(tab=utab,n=as.vector(table(tab.index)),
         row.index = tab.index)}
