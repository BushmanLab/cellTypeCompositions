##' Auxiliary Gibbs Sampler for Dirichlet Process Prior 
##'
##' The auxiliary Gibbs sampler for a Dirichlet Process prior is
##' implemented.
##' @title Dirichlet Process Prior Gibbs Sampler
##' @param wtab the result of \code{\link{uniTab}}
##' @param om the object returned by \code{\link{exOGS}(paramList)}
##' @param eta (optional) initialization and workspace for the sampled
##'     parameters
##' @param etaN (optional) initialization and workspace for the counts
##'     of sampled parameter vectors
##' @param dataToEta index from the original data to the value of eta
##'     sampled
##' @param etaM number of samples provided in the intialization
##' @param auxM how many auxiliary paramters to draw
##' @param alpha how much weight to place on new samples
##' @param verbose whether to report intermediate actions
##' @param etaCols how large to make the workspace (too small a value
##'     will result in r=termination with an error)
##' @return list with elements \code{eta}, \code{etaN], \code{etaM},
##'     and \code{dataToEta} which are updates to the correspondingly
##'     named inputs for \code{gibbsDPP} and a list of two elements,
##'     \code{last} and \code{launch}, given the last two such lists
##'     for \code{gibbsCan}
##' @export
##' @author Charles Berry
gibbsDPP <- function(
                     wtab,
                     om, 
                     eta=NULL,
                     etaN=NULL,
                     dataToEta=NULL,
                     etaM = 0L,
                     auxM = 5L,
                     alpha = 100.0,
                     verbose = 0L,
                     etaCols = 500L
                     )
{
  if (is.null(eta))
    eta <- array(0.0,c(nrow(om), etaCols))
  if (is.null(etaN)) etaN <- rep(0L, ncol(eta))
  if (is.null(dataToEta)) dataToEta <- 
                            rep(-1L,length(wtab$data.index))
  stopifnot( all( etaM >= dataToEta + 1L) )
  auxGibbs(wtab, om, eta, etaN, dataToEta, etaM,
           auxM , alpha , verbose )
}


##' @rdname gibbsDPP
##' @param nscans integer with value \code{nscans>=1L}. How many scans
##'     to do.
##' @param ctParms list of values for \code{nburn}, \code{ncore}, and
##'     \code{method} to pass to \code{\link{ctSampler}}
##' @param ijvals \code{0L} by default, but for splitting and merging
##'     a value of 2 is needed. Users should not usually change the
##'     default.
##' @param ... unused at present
##' @export
gibbsScan <- function(wtab,
		      om,
		      eta=NULL,
		      etaN=NULL,
		      dataToEta=NULL, etaM=0L,
		      alpha=1.0,auxM=10L,
		      ctParms=NULL,
		      nscans = 1L,
		      etaCols=500L,
		      ijvals=0L,
		      verbose=FALSE,
		      ...){

  if (is.null(eta))  eta <- array(0.0,c(nrow(om),etaCols))
  if (is.null(dataToEta)) dataToEta <- rep(-1L,length(wtab[["data.index"]]))
  if (is.null(eta))
    eta <- array(0.0,c(nrow(om), etaCols))
  if (is.null(etaN)) etaN <- rep(0L, ncol(eta))
  if (is.null(dataToEta)) dataToEta <- 
			    rep(-1L,length(wtab$data.index))
  stopifnot( all( etaM >= dataToEta + 1L) )
  stopifnot(length(etaN)==ncol(eta),
	    nrow(om) == ncol(wtab[["tab"]]))

  if (ijvals!=0){
    if (etaM!=ijvals) warning("ijvals != etaM seems wrong")
    if (ncol(eta)!=ijvals) warning("ncol(eta) != ijvals is usually an error")
    if (is.null(dataToEta) || any(dataToEta[1:ijvals]<0) )
      stop("dataToEta[1:ijvals] must be given as non-negative integers")
  } else {
    if (auxM == 0L) warning("auxM == 0 & ijvals == 0 is usually a mistake")
  }
  
  
  ctp <- list(nburn=10L, ncore=1L, method="sampleX")
  ctp[names(ctParms)] <- ctParms

  logpriorC <-
    function(etaN,alpha){
      if (alpha==0.0) return( 0.0 )
      n <- sum(etaN)
      length(etaN)*log(alpha) +
	sum(lgamma(etaN)) -
	sum(log(alpha+0:(n-1)))
    }
  logpost <- function(wtab,om,eta,etaN,dataToEta,etaM,alpha){
    etaN <- etaN[1:etaM]
    eta <- eta[, 1:etaM, drop=FALSE]
    uniq.indexes <- unique(cbind(wtab$data.index, dataToEta))
    uniq.sums <- table(match(paste(wtab$data.index, dataToEta),
			     paste(uniq.indexes[,1, drop=FALSE],
                                   uniq.indexes[,2, drop=FALSE])))
    marglik <-
      marglogpost(eta[,1L+uniq.indexes[,2], drop=FALSE],
                  om,wtab[["tab"]][uniq.indexes[,1], , drop=FALSE])
    margliksum <- sum( marglik * uniq.sums)
    ## Equations 6 and 10 Jain and Neal, 2007 yield: 
    margliksum+logpriorC(etaN,alpha)+etaM*lgamma(nrow(om))
  }


  pass2 <- list(eta=eta,etaN=etaN,dataToEta=dataToEta,etaM=etaM)

  logpst <- rep(0.0,nscans)

  for (iscan in seq_len(nscans)){
    ## if there are updates, do them:
    pass1 <- pass2   
    pass2  <- auxGibbs(wtab,om,
		       pass1[["eta"]],
		       pass1[["etaN"]],
		       pass1[["dataToEta"]],
		       etaM=pass1[["etaM"]],
		       auxM=auxM,
		       alpha=alpha,
		       ijvals=ijvals,
		       verbose=verbose)
    

    eta.by.ct <- table(pass2[["dataToEta"]],
                         wtab[["data.index"]]) %*% wtab[["tab"]]

    eta.tuned <- sapply(ctSampler(eta.by.ct,om,
                                  pass2[["eta"]][, 1:pass2[["etaM"]], drop=FALSE ],
                                  nkeep=1, # makes no sense otherwise
                                  nburn=ctp[["nburn"]],
                                  ncores=ctp[["ncore"]],
                                  method=ctp[["method"]]),
                        prop.table)
    
    pass2[["eta"]][, 1:pass2[["etaM"]] ] <- eta.tuned
    logpst[ iscan ] <-
      logpost(wtab, om, pass2[["eta"]], pass2[["etaN"]],
	      pass2[["dataToEta"]],
	      pass2[["etaM"]],
	      alpha)
    
  }
  
  list(last=pass2, launch=pass1, logpost=logpst)
  
}



