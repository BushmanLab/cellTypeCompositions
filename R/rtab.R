##' Sample Cell Type Table
##'
##' Given the \code{ctScan} fit of a dataset, create a sample table
##' from the posterior.
##' @title rtab
##' @param scan the object produced by \code{\link{gibbsScan}}
##' @param om a filtration matrix such as made by \code{\link{omega}}
##' @param elt which element of the scan to use in sampling.
##' @param tol lower limit on probability of observing a clone
##' @return table of counts 
##' @importFrom stats rpois rnbinom ppois
##' @export
##' @author Charles Berry
  rtab <- function(scan, om, elt=length(scan),tol=1e-3){
    elt <- scan[[elt]]
    etaLambdaN <- table(elt[["dataToEta"]], elt[["dataToLambda"]])
    rho <- t(elt[["eta"]])%*%om
    rhoPlusLambda <- rowSums(rho)%o%elt[["lambda"]]
    is.etaLambda <- etaLambdaN != 0
    if( any(rhoPlusLambda[is.etaLambda] < tol) )
        stop("probability of observing clone(s) < tol")
    etaLambdaN0 <-
      rnbinom(sum( is.etaLambda ),
	      etaLambdaN[is.etaLambda],
	      ppois(0,rhoPlusLambda[is.etaLambda],lower.tail=FALSE))
      etaLambdaN[is.etaLambda] <- etaLambdaN[is.etaLambda] + etaLambdaN0
    rhoLambda <- matrix(t(rho) %o% elt[["lambda"]],nrow=ncol(om))
    tab <- rpois( sum(etaLambdaN)*ncol(om),
		 rhoLambda[, rep(1:length(etaLambdaN),
				 etaLambdaN)])
    dim(tab) <- c(ncol(om),sum(etaLambdaN))
    t(tab[,colSums(tab) != 0L])
  }
