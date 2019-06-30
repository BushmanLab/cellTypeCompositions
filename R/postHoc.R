## post hoc adjustments to gibbsScan objects


##' Impute clones or cells in the parent sample from a \code{"gibbsScan"} object.
##'
##' \code{etaNparent} uses the fitted values of
##' \code{lapply(object[elts],"[[","eta")} and
##' \code{lapply(object[elts],"[[","dataToEta")} and the numbers of
##' cells observed in \code{wtab} to impute the total number of clones
##' with corresponding values of \sQuote{eta} in the parent sample for
##' each observed clone.
##'
##'\code{etaRparent} uses the fitted values of
##' \code{lapply(object[elts],"[[","eta")} and
##' \code{lapply(object[elts],"[[","dataToEta")} and the numbers of
##' cells observed in \code{wtab} to impute the total number of cells
##' with corresponding values of \sQuote{eta} in the parent sample for
##' each observed clone.
##' 
##' @title Clones or cells in the parent sample.
##' @param object a \code{"gibbsScan"} object
##' @param wtab the object produced by \code{\link{uniTab}}
##' @param plist the argument of \code{\link{exOGS}}
##' @param elts a vector of elements to select from \code{object} or
##'     \code{NULL} (the default) to select all elements.
##' @return list of vectors of counts sof cells or clones
##' @importFrom stats dnbinom pbinom qnbinom xtabs
##' @export
##' @author Charles Berry
etaNparent <- function(object,wtab,plist, elts=NULL){
    prDisc <- function(r,rho){
        1.0 - (1-rho)^r * ((rho) /(1-(1-rho)^2))^(r+1)
    }
    if (is.null(elts)) elts <- seq_along(object)
    om <- exOGS(plist)
    r <- rowSums(wtab[["tab"]])
    di <- wtab[["data.index"]]
    lapply(elts, function(ielt){
        elt <- object[[ielt]]
        d2e <- elt[["dataToEta"]]
        rho <- colSums(crossprod(om, elt[["eta"]]))
        pdiscinv <- 1 / prDisc(r[di], rho[d2e] )
        disctab <- 
            xtabs( pdiscinv ~ d2e)
        as.vector(disctab)
    })
}

##' @rdname etaNparent
##' @param eps The tolerance to use in selecting plausible values of a
##'     negative binomial
##' @export
etaRparent <- function(object,wtab,plist, elts=NULL,eps=1e-5){
    if (is.null(elts)) elts <- seq_along(object)
    om <- exOGS(plist)
    r <- rowSums(wtab[["tab"]])
    n <- wtab[["n"]]
    di <- wtab[["data.index"]]
    uniq.r <- unique(r)
    r.index <- match(r,uniq.r)
    lapply(elts, function(ielt){
        elt <- object[[ielt]]
        d2e <- elt[["dataToEta"]]
        rho <- colSums(crossprod(om, elt[["eta"]]))    
        index <- cbind(r.index[di],d2e)
         uind <- unique(index)
        totals <- apply(uind,1,function(z){
            r <- uniq.r[z[1]]
            rho <- rho[z[2]]
            lowerx <- qnbinom(eps,r,rho,lower.tail=TRUE)
            upperx <- qnbinom(eps,r,rho,lower.tail=FALSE)
            x <- lowerx:upperx
            prx <- dnbinom(x,r,rho)
            N  <- r + x
            adjN <- prx  / pbinom(0, N , rho, lower.tail=FALSE)
            sum( adjN * N )
        })
        totmat <- array(0.0,c( length(uniq.r), elt[["etaM"]]))
        totmat[ uind ] <- totals
        totmat[ index ]
    })
    }
