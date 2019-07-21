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
##' @param sample logical. If \code{TRUE}, sample the unobserved
##'     clones (\code{etaNparent}) or cells (\code{etaRparent}) else
##'     estimate the number.
##' @param perEta logical. If \code{TRUE} accumulate the clones
##'     acording to the columns of \code{object[[elt]][["eta"]]} else
##'     return the vector of counts.
##' @return list of vectors of counts of cells or clones
##' @importFrom stats dnbinom pbinom qnbinom rnbinom xtabs
##' @export
##' @author Charles Berry
etaNparent <- function(object,wtab,plist, elts=NULL,sample=FALSE,perEta=TRUE){
    prDisc <- function(r.vec,rho.vec,eps=1e-8){
        stopifnot(length(r.vec)== length(rho.vec))
        sapply(seq_along(r.vec),function(i){
            r <- r.vec[i]
            rho <- rho.vec[i]
            lower <- qnbinom(eps,r+1,rho)
            if (pbinom(0,r+lower,rho, lower.tail=FALSE,log.p=TRUE) > -eps)
                return(1.0)
            upper <- qnbinom(eps,r+1,rho,lower=FALSE) 
            sum( dnbinom(lower:upper,r+1,rho) /
                 pbinom(0, r+(lower:upper), rho, lower.tail=FALSE))
        })}
    if (is.null(elts)) elts <- seq_along(object)
    om <- exOGS(plist)
    r <- rowSums(wtab[["tab"]])
    di <- wtab[["data.index"]]
    lapply(elts, function(ielt){
        elt <- object[[ielt]]
        d2e <- elt[["dataToEta"]]
        rho <- colSums(crossprod(om, elt[["eta"]]))
        prd <- prDisc(r[di], rho[d2e])
        discwt <- if (sample){
                      1+rnbinom(length(d2e), 1L, prd)
                  } else {
                      1 / prd
                  }
        
        if (perEta) {
            disctab <- xtabs( discwt ~ d2e )
            as.vector(disctab)
        } else {
            discwt
        }
        })       
}

    
##' @rdname etaNparent
##' @param eps The tolerance to use in selecting plausible values of a
##'     negative binomial
##' @param inflateN \code{NULL} or a list who canonical element is a
##'     vector of \code{length(wtab[["data.index"]])} integers to
##'     indicate how many times to repeat each observed clone if
##'     \code{isTRUE(sample)}.
##' @importFrom stats rnbinom dnbinom qnbinom
##' @export
etaRparent <- function(object,wtab,plist, elts=NULL,eps=1e-5,
                       sample=FALSE,inflateN=NULL){
    if (is.null(elts)) elts <- seq_along(object)
    om <- exOGS(plist)
    r <- rowSums(wtab[["tab"]])
    n <- wtab[["n"]]
    di <- wtab[["data.index"]]
    if (sample){
        do_inflate <- !is.null(inflateN)
        lapply(elts, function(ielt) {
            elt <- object[[ielt]]
            d2e <- elt[["dataToEta"]]
            if (do_inflate) {
                stopifnot(length(inflateN[[ielt]]) == length(wtab[["data.index"]]))
                di <- rep(di,inflateN[[ielt]])
                d2e  <- rep(d2e,inflateN[[ielt]])
            }
            rho <- colSums(crossprod(om, elt[["eta"]]))
            r[di]+rnbinom(length(di), r[di]+1, rho[d2e] )
        })
    } else {
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
                lowerx <- qnbinom(eps,r+1,rho,lower.tail=TRUE)
                upperx <- qnbinom(eps,r+1,rho,lower.tail=FALSE)
                x <- lowerx:upperx
                prx <- dnbinom(x,r+1,rho)
                N  <- r + x
                adjN <- prx  / pbinom(0, N , rho, lower.tail=FALSE)
                sum( adjN * N )
            })
            totmat <- array(0.0,c( length(uniq.r), elt[["etaM"]]))
            totmat[ uind ] <- totals
            totmat[ index ]
        })
    }
}

