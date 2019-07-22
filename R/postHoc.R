##' posterior samples from gibbsScan objects
##'
##' Sample clones and cell numbers in the parent sample from a
##' \code{"gibbsScan"} object.
##'
##' \code{rInvisibles} uses the fitted values of
##' \code{lapply(object[elts],"[[","eta")} and
##' \code{lapply(object[elts],"[[","dataToEta")} and the numbers of
##' cells observed in \code{wtab} to sample the total number of clones
##' with corresponding values of \sQuote{eta} in the parent sample
##' that were not observed and the number of cells for those clones.
##'
##' @title Clones and cells in the parent sample.
##' @param object a \code{"gibbsScan"} object
##' @param wtab the object produced by \code{\link{uniTab}}
##' @param plist the argument of \code{\link{exOGS}}
##' @param elts a vector of elements to select from \code{object} or
##'     \code{NULL} (the default) to select all elements.
##' @return list with each element a list with vector \code{R} of
##'     counts of cells and vector \code{n.invis} a count of clones
##'     not observed.
##' @author Charles Berry
##' @importFrom stats rnbinom
##' @export
rInvisibles <- function(object,wtab,plist, elts=NULL){
    if (is.null(elts)) elts <- seq_along(object)
    om <- exOGS(plist)
    r <- rowSums(wtab[["tab"]])
    n <- wtab[["n"]]
    di <- wtab[["data.index"]]
    rdi <- r[di]
    lapply(elts, function(ielt) {
            elt <- object[[ielt]]
            d2e <- elt[["dataToEta"]]
            rho <- colSums(crossprod(om, elt[["eta"]]))
            R.samp <- rdi + rnbinom(length(di), rdi+1, rho[d2e])
            n.invis <- rnbinom(length(R.samp), 1, 1 - rho[d2e]^R.samp)
            list(R=R.samp, n.invis = n.invis )})
}
