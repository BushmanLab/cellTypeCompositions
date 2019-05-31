##' Log Likelihood of Poisson observations given rates, proportions
##' lost, and mis-sorted
##'
##' \code{poislogpost} gives the log likelihood for a collection of
##' posterior samples from a Poisson parent for the expected numbers
##' of observations. \code{marglogpost} implicitly integrates over a
##' uniform proper prior for the expected total in a sample and
##' returns values that are proportional to the posterior.
##' @title Log Likelohood of Poisson Mixture
##' @param ex.Sample matrix of expected sample counts - one column for
##'     each clone. For \code{marglogpost} can be proportions.
##' @param ex.OGS expected cell counts observed given sample counts of
##'     one
##' @param obs matrix of observed counts with one row for each clone
##' @return \code{nrow(obs)} loglikelihoods
##' @author Charles Berry
##' @export
##' @importFrom stats dpois
poislogpost <- function(ex.Sample,ex.OGS,obs){
    stopifnot( length(obs) == ncol(ex.OGS) )
    lambda <- t(ex.OGS) %*% ex.Sample
    res <-   dpois(obs,lambda,log=TRUE)
    if (is.matrix(lambda)) {
        dim(res) <- dim(lambda)
        colSums( res )
    } else {
        sum(res)
    }
    
}


## internal multinomial mass funciton
dmulti <- function(x,prob,log=FALSE){
    stopifnot(all(prob>0))
    prob <- as.matrix(prob)
    cs <- colSums(prob) 
    stopifnot(cs <= 1.0 + 1e-7)
    xs <- colSums(x*log(prob) - lgamma(x+1))+lgamma( colSums(as.matrix(x)) + 1 )
    if (log) xs else exp(xs)
}

##' @rdname poislogpost
##' @export
marglogpost <- function(ex.Sample,ex.OGS,obs){
    ## obs can be a vector - recycled to match ex.Sample
    ## or it can be a matrix with dim(obs) == rev( dim(ex.Sample) )
    rho <- t( ex.OGS ) %*% ex.Sample
    psum <- colSums(rho)
    prop <- prop.table(rho,2)
    if (is.matrix(obs)) obs <- t( obs )
    - log( psum ) + dmulti( obs, prop, log=TRUE)
  }


##' Compute the logposterior from an object produced by \code{ctSampler}
##'
##' The details of the setup are inferred from
##' \code{attr(obj,"call")}, then \code{\link{poislogpost}} is called
##' on each sample to find its log posterior.  A typical use is to
##' wrap a call to \code{ctSampler} with this function.  If the
##' objects used to construct the object \code{obj} are not available,
##' the function will fail.
##' @title log posterior of Gibbs samples
##' @param obj see \code{\link{ctSampler}}
##' @param integrate logical, use marginal posterior over a locally
##'     uniform prior (in which case the result is proportional to the
##'     posterior).
##' @return a list of vectors of log posteriors
##' @author Charles Berry
##' @export
logposterior <- function(obj,integrate=FALSE){
    objcall <- attr( obj, "call" )
    wtab <- eval.parent(objcall$gMat)
    if (!is.matrix(wtab))
            dim(wtab) <- c(1, length(wtab))
    uop <- eval.parent(objcall$uop)
    pfun  <- if (integrate) marglogpost else poislogpost
    lapply(1:nrow(wtab), function(i) pfun(obj[[i]],uop,wtab[i,]))
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

