##' Log Likelihood of Poisson observations given rates, proportions
##' lost, and mis-sorted
##'
##' 
##' @title Log Likelohood of Poisson Mixture 
##' @param ex.Sample matrix of expected sample counts - one column for each clone
##' @param ex.OGS expected cell counts observed given sample counts of one
##' @param obs matrix of observed counts with one row for each clone
##' @return \code{nrow(obs)} loglikelihoods
##' @author Charles Berry
##' @export
##' @importFrom stats dpois
poislogpost <- function(ex.Sample,ex.OGS,obs){
    stopifnot( length(obs) == nrow(ex.Sample) )
    lambda <- t(ex.OGS) %*% ex.Sample
    res <-   dpois(obs,lambda,log=TRUE)
    dim(res) <- dim(lambda)
    colSums( res )
  }

##' Compute the logposterior from an object produced by \code{ctSampler}
##'
##' The details of the setup are inferred from
##' \code{attr(obj,"call")}, then \code{\link{poislogpost}} is
##' called on each sample to find its log posterior.  A typical use is
##' to wrap a call to \code{ctSampler} with this function.  If the
##' objects used to construct the object \code{obj} are not available,
##' the function will fail.
##' @title log posterior of Gibbs samples
##' @param obj see \code{\link{ctSampler}}
##' @return a list of vectors of log posteriors
##' @author Charles Berry
##' @export
logposterior <- function(obj){
    objcall <- attr( obj, "call" )
    wtab <- eval.parent(objcall$gMat)
    uop <- eval.parent(objcall$uop)
    lapply(1:nrow(wtab), function(i) poislogpost(obj[[i]],uop,wtab[i,]))
}


##' Convert Expected Sample Counts to Expected Observation Counts
##'
##' 
##' @title  Expected Observation Counts
##' @param pl list with elements \code{upsilon}, \code{omega}, and \code{psi} 
##' @return matrix of multipliers
##' @author Charles Berry
##' @export
exOGS <- function(pl) {
    with(pl,diag(upsilon)%*%omega%*%diag(psi))
}
##' Detection Probabilities for Clones Given Their Sample Expectations
##'
##' 
##' @title Probability of Detection
##' @param ex.Sample matrix of expected sample counts - one column for each clone
##' @param pl list with elements \code{upsilon}, \code{omega}, and \code{psi} 
##' @param subset vector of rows of \code{oxOGS(pl)} represented in \code{ex.Sample}
##' @return vector of detection probabilities
##' @author Charles Berry
##' @xport
##' @importFrom stats ppois
pod <- function(ex.Sample,pl,subset=1:nrow(ex.Sample)){
    om <- exOGS(pl)
    lambda <- rowSums(om[subset,])%*%ex.Sample
    ppois(0,lambda,lower.tail=FALSE)
}

