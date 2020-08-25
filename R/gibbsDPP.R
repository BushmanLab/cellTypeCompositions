##' Auxiliary Gibbs Sampler for Dirichlet Process Prior 
##'
##' The auxiliary Gibbs sampler for a Dirichlet Process prior is
##' implemented. The usual call is to \code{gibbsScan}, but
##' \code{update} can be used to continue a chain (with no arguments)
##' possibly changing the parameters.
##'
##' \code{gibbsDPP} only does an initial Gibbs scan with no following
##' tuning of parameters or subsequent iteration.  It may be helpful
##' for low-level access to the scanning routine, but rarely is used.
##'
##' \code{gibbsScan} is the usual entry point.  It conducts the Gibbs
##' scan for \code{eta} and \code{lambda} and subseqeunt tuning alone the
##' lines of Algorithm 8 in Neal, 2000, but the scan for \code{lambda}
##' uses the integrated posterior from the base distribution as in
##' Algorithm 1 of Neal, 2000 to decide whether to add a new value.
##'
##' The method for \code{update} allows iterations under
##' \code{gibbsScan(..., keep = TRUE)} to resume. If another value is
##' used for \code{keep} some care is needed to mimic the return value
##' as described next.
##'
##' \code{keep.default} is used internally under the default
##' (\code{gibbsScan(..., keep = TRUE)}) to return values from each
##' iteration. Users can provide their own function to alter the
##' collection of values from each iterate. The function should have
##' the same number of arguments as \code{keep.default} and may pass
##' those arguments to \code{keep.default} when \code{i == nkeep} to
##' allow \code{update} to resume iteration.
##' @title Dirichlet Process Prior Gibbs Sampler
##' @param wtab the result of \code{\link{uniTab}}
##' @param om the object returned by \code{\link{exOGS}(paramList)}
##' @param eta (optional) initialization and workspace for the sampled
##'     parameters
##' @param etaN (optional) initialization and workspace for the counts
##'     of sampled parameter vectors
##' @param dataToEta index from the original data to the value of eta
##'     sampled
##' @param lambda (optional) initialization and workspace for the
##'     sampled abundance parameters
##' @param lambdaN (optional) initialization and workspace for the
##'     counts of sampled parameters
##' @param dataToLambda index from the original data to the value of
##'     lambda sampled
##' @param etaM number of samples provided in the intialization
##' @param auxEtaM how many auxiliary paramters to draw
##' @param alphaEta how much weight to place on new samples
##' @param lambdaM how many auxiliary paramters to draw
##'
##' @param auxLambdaM how much space to allot for new samples. Only
##'     one new sample will be drawn if greater than zero
##' @param alphaLambda how much weight to place on new samples
##' @param verbose whether to report intermediate actions
##' @param etaCols how large to make the workspace (too small a value
##'     will result in termination with an error)
##' @param lambdaSize dprior prior value for draws of eta
##' @param dprior symmetric Dirichlet prior value
##' @param lambdaShape base gamma distribution prior shape
##' @param lambdaRate base distribution prior rate
##' @param ijvals \code{0L} by default, but for (not yet implemented)
##'     splitting and merging a value of 2 is needed. Users should not
##'     usually change the default.
##' @param ... currently unused
##' 
##' @return list which by default has elements \code{eta}, \code{etaN}, \code{etaM},
##'     \code{lambda}, \code{lambdaN}, \code{lambdaM}, and
##'     \code{dataToEta} and \code{dataToLambda} which are updates to
##'     the correspondingly named inputs for \code{gibbsDPP} and a
##'     list of lists of such elements for \code{gibbsScan}. The
##'     element \code{logLik} is the loglikelihood for \code{lambda}
##'     and \code{eta} conditioning on the values of all other
##'     parameters.
##'
##' @references Radford M. Neal (2000) \emph{Markov Chain Sampling
##'     Methods for Dirichlet Process Mixture Models}, Journal of
##'     Computational and Graphical Statistics, 9:2, 249-265, DOI:
##'     10.1080/10618600.2000.10474879

##' @export
##' @author Charles Berry
gibbsDPP <- function(
                     wtab,
                     om, 
                     eta = NULL,
                     etaN = NULL,
                     dataToEta = NULL,
                     lambda = NULL,
                     lambdaN = NULL,
                     dataToLambda = NULL,
                     etaM = 0L,
                     auxEtaM = 5L,
                     alphaEta = 100.0,
                     lambdaM = 0L,
                     auxLambdaM = 5L,
                     alphaLambda = 1.0,
                     verbose = 0L,
                     etaCols = 500L,
                     lambdaSize = 10L,
                     dprior=1.0,
                     lambdaShape=1.0,
                     lambdaRate=0.01,
                     ijvals=0L,
                     ...
                     )
{
  mc <- match.call()
  if (is.null(eta))
    eta <- array(0.0,c(nrow(om), etaCols))
  if (is.null(etaN)) etaN <- rep(0L, ncol(eta))
  if (is.null(dataToEta)) dataToEta <- 
                            rep(0L,length(wtab$data.index))
  stopifnot( all( etaM >= dataToEta) )

  if (is.null(lambda))
    lambda <- rep(as.double(NA), lambdaSize)
  if (is.null(lambdaN)) lambdaN <- rep(0L, length(lambda))
  if (is.null(dataToLambda)) dataToLambda <-
                               rep(0L,length(wtab$data.index))
  stopifnot( all( lambdaM >= dataToLambda ) )

  tab <- I( wtab[["tab"]] ) # force copy
  di <- wtab[["data.index"]] - 1L
  dataToEta <- dataToEta - 1L
  dataToLambda <- dataToLambda - 1L

  res <- auxGibbs(tab, di,
                  om,
                  eta,
                  etaN,
                  dataToEta,
                  lambda,
                  lambdaN,
                  dataToLambda,
                  etaM,
                  auxEtaM,
                  alphaEta,
                  lambdaM,
                  auxLambdaM,
                  alphaLambda,
                  ijvals,
                  verbose,
                  dprior,
                  lambdaShape,
                  lambdaRate)

  
  dim(res[["etaN"]]) <- dim(res[["dataToEta"]]) <- dim(res[["lambda"]]) <-
    dim(res[["lambdaN"]]) <- dim(res[["dataToLambda"]]) <- NULL
  res[["dataToEta"]] <-  res[["dataToEta"]] + 1L
  res[["dataToLambda"]] <-  res[["dataToLambda"]] + 1L
  res[["call"]] <- mc
  class(res) <- "gibbsDPP"
  res
}


##' @rdname gibbsDPP
##' @param nkeep integer with value \code{nkeep>=1L}. How many scans
##'     to do.
##' @param nthin how many iterations per saved value
##' @param nburn discard this many iterations before saving any values
##' @param keep vector of elements of the result in each iteration to
##'   retain or a function called with \code{pass}, the intermediate results, 
##'
##' @param ... currently unused
##' @param abEta if \code{NULL} treat \code{alphaEta} as a fixed
##'     value. Otherwise the first two elements are the shape and
##'     rate parameters of Gamma prior for \code{alphaEta}.
##' @param abLambda akin to \code{abEta}
##' @param niter.tune how many cycles of tuning for each gibbs scan
##'
##' @importFrom stats rbeta rgamma dgamma runif
##' @export
##' @examples
##' tab <- diag(10, nrow = 3)[ rep(1:3,2), ] +
##'   rbinom(18,2,0.5)
##' uop <- prop.table(diag(3)+0.05,1)
##' gibbsScan(uniTab(tab), uop)
gibbsScan <- function(wtab,
		      om,
		      eta = NULL,
		      etaN = NULL,
		      dataToEta = NULL, etaM = 0L,
		      alphaEta = 1.0,auxEtaM = 10L,
                      lambda = NULL,
                      lambdaN = NULL,
                      dataToLambda = NULL, lambdaM = 0L,
                      alphaLambda = 1.0, auxLambdaM = 10L,
		      nkeep  =  1L, nthin = 1L, nburn = 0L,
		      etaCols = 50L,
                      lambdaSize = 20L,
		      ijvals = 0L,
                      abEta = c(0.0001,0.0001),
                      abLambda = c(0.01,0.01),
		      verbose = FALSE, dprior = 1.0,
                      lambdaShape = 1.0, lambdaRate=0.01,
                      keep = TRUE, niter.tune=10L,
		      ...){
  ## define helper functions:

  ## Escobar and West, Page 10  (Chapter in Dey et al, 2012)
  ralpha <- function(alpha,a,b,Istar,I){
    eta <- rbeta(1,alpha+1,I)
    shape <- a + Istar - 1
    rate <- (b - log(eta))
    plus1 <- runif(1, -shape, I*rate) < 0.0
    rgamma(1,shape+plus1, rate = rate)
  }
    
  mc <- match.call()

  stopifnot(all(om>=0.0))
  if (any(om==0.0)){
    om[om==0.0] <- .Machine[["double.eps"]]
    warning("Converted zeroes in om to machine epsilon")
  }
  if (is.null(eta))  eta <- array(0.0,c(nrow(om),etaCols))
    if (is.null(dataToEta)){
        dataToEta <- rep(-1L,length(wtab[["data.index"]]))
    } else {
        stopifnot(length(wtab[["data.index"]])==length(dataToEta))
    }    
  if (is.null(eta))
    eta <- array(0.0,c(nrow(om), etaCols))
  if (is.null(etaN)) etaN <- rep(0L, ncol(eta))
  stopifnot( all( etaM >= dataToEta) )
  stopifnot(length(etaN)==ncol(eta),
            nrow(om) == ncol(wtab[["tab"]]))
  stopifnot(nthin >= 1L)
  stopifnot(all(abEta>0), length(abEta)==0 || length(abEta==2))
  stopifnot(all(abLambda>0), length(abLambda)==0 || length(abLambda==2))
  if (ijvals!=0){
    if (etaM!=ijvals) stop("ijvals != etaM or 0L is not permitted")
      if (is.null(dataToEta) || any(dataToEta[1:ijvals]<1) )
      stop("dataToEta[1:ijvals] must be given as positive integers")
  }
  
  if (auxEtaM != 0L && etaM + auxEtaM > ncol(eta)){
    ## need to pad eta and etaN
    padby  <- max(etaM+auxEtaM, etaCols) - ncol(eta)
    eta  <- cbind(eta,array(0.0, c(nrow(eta), padby)))
    etaN <- c(etaN, rep(0.0, padby))
  }

  if (is.null(lambda))
    lambda <- rep(as.double(NA), lambdaSize)
  if (is.null(lambdaN)) lambdaN <- rep(0L, length(lambda))
    if (is.null(dataToLambda)) {
        dataToLambda <- rep(0L,length(wtab[["data.index"]]))
    } else {
        stopifnot(length(dataToLambda)==length(wtab[["data.index"]]))
    } 
  stopifnot( all( lambdaM >= dataToLambda ) )
  if (auxLambdaM != 0L && lambdaM + auxLambdaM > length(lambda)){
    ## need to pad lambda and lambdaN
    padby  <- max(lambdaM+auxLambdaM, lambdaSize) - length(lambda)
    lambda  <- c(lambda,rep(0.0, padby))
    lambdaN <- c(lambdaN, rep(0.0, padby))
  }

  if (is.function(keep)){
    keep.fun <- keep
  } else {
    keep.fun <- keep.default
  }
  
  tab <- I( wtab[["tab"]] ) # force copy
  di <- wtab[["data.index"]] - 1L
  dataToEta <- dataToEta - 1L
  dataToLambda <- dataToLambda - 1L

  N <- sum(wtab[["n"]])
  log.dalphaEta <- 0.0
  
  pass2 <- list(eta = eta,
                etaN = etaN,
                dataToEta = dataToEta,
                lambda = lambda,
                lambdaN = lambdaN,
                dataToLambda = dataToLambda, 
                etaM = etaM,
                lambdaM = lambdaM,
                alphaEta = alphaEta,
                alphaLambda = alphaLambda)

  keepers <- list()
  isamp <- 0L

  iters.per.keeper <- rep(as.integer(nthin), nkeep)
  iters.per.keeper[ 1L ]  <-   iters.per.keeper[ 1L ] + nburn

   ## browser()  
  for (i in 1:nkeep){
    for (j in 1:iters.per.keeper[ i ]){
      ## if there are updates, do them:
      pass1 <- pass2   
      pass2  <- auxGibbs(tab,di, om,
                         pass1[["eta"]],
                         pass1[["etaN"]],
                         pass1[["dataToEta"]],
                         pass1[["lambda"]],
                         pass1[["lambdaN"]],
                         pass1[["dataToLambda"]],
                         pass1[["etaM"]],
                         auxEtaM,
                         alpha = pass1[["alphaEta"]],
                         pass1[["lambdaM"]],
                         auxLambdaM,
                         pass1[["alphaLambda"]],
                         ijvals = ijvals,
                         verbose = verbose,
                         dprior = dprior)
      

        sparm <- 
            sampleParms(tab, di, om,
                        pass2[["dataToEta"]],
                        pass2[["dataToLambda"]],
                        pass2[["eta"]],
                        pass2[["etaM"]],
                        pass2[["lambda"]],
                        pass2[["lambdaM"]],
                        dprior, lambdaShape, lambdaRate, niter.tune, verbose)

        pass2[["eta"]] <- sparm[["eta"]]
        pass2[["lambda"]] <- sparm[["lambda"]]
        
      ## update alpha*
      if (!is.null(abEta)) {
        pass2[["alphaEta"]] <-
          ralpha(pass1[["alphaEta"]],abEta[1],abEta[2],
                 pass2[["etaM"]], N)
        log.dalphaEta <- dgamma(pass2[["alphaEta"]], abEta[1],
                                abEta[2], log = TRUE)
      }
      if (!is.null(abLambda)) {
        pass2[["alphaLambda"]] <-
          ralpha(pass1[["alphaLambda"]], abLambda[1],abLambda[2],
                 pass2[["lambdaM"]], N)
        log.dalphaLambda <- dgamma(pass2[["alphaLambda"]],
                                   abLambda[1], abLambda[2], log = TRUE)
      } else {
        pass2[["alphaLambda"]] <- alphaLambda
        log.dalphaLambda <- 0.0
        }
    }
    

    log.lik <- sparm[["logLik"]] - sum(wtab[["n"]]*lfactorial(wtab[["tab"]]))
         
    keepers[[ i ]] <- keep.fun(pass2, log.lik, i, nkeep )             
  }
  attr(keepers,"call") <- mc
  class(keepers) <- "ctScan"
  keepers
}


##' @rdname gibbsDPP
##' @param object a \dQuote{ctScan} object
##' @param elt integer to select a sample to use as the starting point
##' @param ... parameters to revise using \code{update} 
##' @method update ctScan
##' @export
update.ctScan <- function(object, elt = length(object), ...){
  mc <- match.call()
  nameOnly <- c("eta","etaN","dataToEta","etaM","alphaEta",
                "lambda","lambdaN","dataToLambda","lambdaM","alphaLambda")
  objcall <- attr(object,"call")
  newparms <- list(...)
  newnames <- names(newparms)
  newind <- match(newnames,names(objcall),0L)
  if (any(newnames == 0L)) stop( "couldn't match ", newnames[newind==0L])
  objcall[newnames] <- newparms
  nameOnly <- setdiff(nameOnly, newnames)
  if (is.null(elt)) elt <- length(object)
  elt  <- object[[elt]]
  for ( i in nameOnly) objcall[[i]] <- as.name(i)
  objcall[["update.call"]] <- mc
  eval(objcall,elt, parent.frame())
}

##' @rdname gibbsDPP
##' @param pass a list of parameter values following an iteration, viz.
##' the values in elements of the default return value of \code{gibbsScan}.
##' 
##' @param log.posterior the loglikelihood evaluated at the current
##'     values of \code{pass}
##' @param i the iteration number
##' @param nkeep the maximum number of iterations
##' @export
keep.default <-
    function(pass, log.posterior, i, nkeep){
        with(pass,
             list(eta = eta[, 1:etaM],
                  etaN = etaN[ 1:etaM],
                  dataToEta = as.vector(dataToEta + 1L),
                  etaM = etaM,
                  lambda = lambda[1:lambdaM],
                  lambdaN = lambdaN[ 1:lambdaM],
                  dataToLambda = as.vector(dataToLambda + 1L),
                  lambdaM = lambdaM,
                  alphaEta = alphaEta,
                  alphaLambda = alphaLambda,
                  logLik = log.posterior
                  ))
    }
