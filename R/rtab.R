##' Sample Cell Type Table
##'
##' Given the \code{ctScan} fit of a dataset, create a sample table
##' from the posterior.
##' @title rtab
##'
##' @param scan the object produced by \code{\link{gibbsScan}} or a
##'     list whose last element has components \dQuote{dataToEta},
##'     \dQuote{dataToLambda}. \dQuote{eta}, and \dQuote{lambda} akin
##'     to those of a \code{ctScan} object.
##' 
##' @param om a filtration matrix such as made by \code{\link{omega}}
##' @param elt which element of the scan to use in sampling.
##' @param tol lower limit on probability of observing a clone
##' @param impute.unseen if \code{TRUE} unseen clones will be
##'     sampled. This might be \code{FALSE} if \code{scan} represents
##'     a population that is to be sampled.
##' @return table of counts
##' @importFrom stats rpois rnbinom ppois
##' @examples
##' ascan <- list()
##' ascan$dataToEta <- rep(1:7,each=3)
##' ascan$dataToLambda <- rep(1:3, 7)
##' ascan$lambda <- c( 1.0, 5.0, 20.0)
##' ascan$eta <- cbind( diag(3) * 0.7 + 0.1, 1/3, -0.3 * diag(3)+0.4)
##' ascan <- list( ascan )
##' om <-  diag(9:7)/10 + 0.025
##' atab <- rtab(ascan, om)
##' @export
##' @author Charles Berry
  rtab <- function(scan, om, elt=length(scan),tol=1e-3, impute.unseen=TRUE){
    elt <- scan[[elt]]
    etaLambdaN <- table(elt[["dataToEta"]], elt[["dataToLambda"]])
    rho <- t(elt[["eta"]])%*%om
    rhoPlusLambda <- rowSums(rho)%o%elt[["lambda"]]
    is.etaLambda <- etaLambdaN != 0
    if( any(rhoPlusLambda[is.etaLambda] < tol) )
        stop("probability of observing clone(s) < tol")
    etaLambdaN0 <-
        if (impute.unseen){
            rnbinom(sum( is.etaLambda ),
                    etaLambdaN[is.etaLambda],
                    ppois(0,rhoPlusLambda[is.etaLambda],lower.tail=FALSE))
        } else {
            0
        }
    etaLambdaN[is.etaLambda] <- etaLambdaN[is.etaLambda] + etaLambdaN0
    rhoLambda <- matrix(t(rho) %o% elt[["lambda"]],nrow=ncol(om))
    tab <- rpois( sum(etaLambdaN)*ncol(om),
		 rhoLambda[, rep(1:length(etaLambdaN),
				 etaLambdaN)])
    dim(tab) <- c(ncol(om),sum(etaLambdaN))
    t(tab[,colSums(tab) != 0L])
  }
