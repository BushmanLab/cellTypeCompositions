##' Counts from the Parent Specimen are Sampled and Confused 
##'
##' Given \code{R[i]} cells with cell type probabilities
##' \code{eta[,i]}, \code{n[i]} draws are made.  Subsequently, these
##' are misclassified and possibly omitted acording to the
##' probabilities in \code{uop}.
##' @title Subsample Counts
##' @param n integer vector. How many draws are to be made for each
##'     column of \code{eta} and element of \code{R}.
##' @param eta matrix with rows giving the probability distribution of
##'     cell type
##' @param R vector of cell numbers
##' @param uop matrix of classification probabilities whose rows sum
##'     to 1.0 or less (when downsampling is in effect).
##' @return matrix of counts with \code{sum(n)} rows
##' @export
##' @importFrom stats rbinom
##' @author Charles Berry
sampleParent <- function(n,eta,R,uop){
    stopifnot(length(n) == ncol(eta), length(n) == length(R))
    stopifnot( uop <= 1.0)
    k <- nrow(eta)
    rho <- t(uop) %*% eta
    if (any(colSums(rho)>1.0)) stop("eta and/or uop misspecified") 
    rhoMiss <- 1 - colSums(rho)
    pr.done <- apply(rho,2,cumsum)
    res <- array(0L,c(sum(n),nrow(eta)))
    R.left <- rep(R, n)
    p.left  <- rep(1.0, sum(n))
    for (i in 1:k){
        rhovec <- rep(rho[ i, ],n)
        res[,i] <- rbinom( sum(n), R.left, pmin(1.0, rhovec/p.left))
        p.left <- p.left - rhovec
        R.left <- R.left - res[ , i ]
    }
    res
}
