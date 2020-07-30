##' Sample Cell Type Table
##'
##' Given the \code{ctScan} fit of a dataset, create a sample table
##' from the posterior.
##' @title rtab
##' @param scan the object produced by \code{\link{gibbsScan}}
##' @param uop a filtration matrix such as made by \code{\link{exOGS}}
##' @param elt which element of the scan to use in sampling.
##' @param tol lower limit on probability of observing a clone
##' @return table of counts 
##' @importFrom stats rpois
##' @export
##' @author Charles Berry
  rtab <- function(scan, uop, elt=length(scan),tol=1e-3){
    elt <- scan[[elt]]
    etaLambdaN <- table(elt[["dataToEta"]], elt[["dataToLambda"]])
    rho <- t(elt[["eta"]])%*%uop
    rhoPlusLambda <- rowSums(rho)%o%elt[["lambda"]]
    is.etaLambda <- etaLambdaN != 0
    if( any(rhoPlusLambda[is.etaLambda] < tol) )
        stop("probability of observing clone(s) < tol")
    etaLambdaN0 <-
      rnbinom(sum( is.etaLambda ),
	      etaLambdaN[is.etaLambda],
	      ppois(0,rhoPlusLambda[is.etaLambda],lower.tail=FALSE))
      etaLambdaN[is.etaLambda] <- etaLambdaN[is.etaLambda] + etaLambdaN0
    rhoLambda <- matrix(t(rho) %o% elt[["lambda"]],nrow=ncol(uop))
    tab <- rpois( sum(etaLambdaN)*ncol(uop),
		 rhoLambda[, rep(1:length(etaLambdaN),
				 etaLambdaN)])
    dim(tab) <- c(ncol(uop),sum(etaLambdaN))
    t(tab[,colSums(tab) != 0L])
  }
