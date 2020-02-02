##' Auxiliary Gibbs Sampler for Dirichlet Process Prior 
##'
##' The auxiliary Gibbs sampler for a Dirichlet Process prior is
##' implemented. The usual call is to \code{gibbsScan}, but
##' \code{update} can be used to continue a chain (with no arguments)
##' possibly changing the parameters.
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
##'     will result in termination with an error)
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
                     eta = NULL,
                     etaN = NULL,
                     dataToEta = NULL,
                     lambda = NULL,
                     lambdaN = NULL,
                     dataToLambda = NULL,
                     etaM = 0L,
                     auxM = 5L,
                     alpha = 100.0,
                     lambdaM = 0L,
                     auxLambdaM = 5L,
                     alphaLambda = 1.0,
                     verbose = 0L,
                     etaCols = 500L,
                     lambdaSize = 10L,
                     dprior=1.0,
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
                  auxM,
                  alpha,
                  lambdaM,
                  auxLambdaM,
                  alphaLambda,
                  ijvals,
                  verbose,
                  dprior)
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
##' @param ctParms list of values for \code{nburn}, \code{ncore}, and
##'     \code{method} to pass to \code{\link{ctSampler}}
##' @param ijvals \code{0L} by default, but for splitting and merging
##'     a value of 2 is needed. Users should not usually change the
##'     default.
##' @param ab if \code{NULL} treat \code{alpha} as a fixed
##'     value. Otherwise the first two elemenbts are rthe shape and
##'     rate parameters of Gamma prior for \code{alpha}.
##' @param dpriors The default value is \code{c(1.0,1.0)}. A smaller
##'     value for \code{dpriors[1]} causes draws from the prior for
##'     \code{eta} to be overdispersed. A larger value for
##'     \code{dpriors[2]} causes posterior samples of \code{eta} to
##'     shrink towards equality. Non-default value can be used to
##'     influence the behavior of the sampler. The computations of the
##'     posterior do not account for non-default values of
##'     \code{dpriors} and should be ignored when non-default values
##'     are used.
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
		      alpha = 1.0,auxM = 10L,
                      lambda = NULL,
                      lambdaN = NULL,
                      dataToLambda = NULL, lambdaM = 0L,
                      alphaLambda = 1.0, auxLambdaM = 10L,
		      ctParms = NULL,
		      nkeep  =  1L, nthin = 1L, nburn = 0L,
		      etaCols = 500L,
                      lambdaSize = 20L,
		      ijvals = 0L,
                      ab = c(0.0001,0.0001),
		      verbose = FALSE, dpriors = c(1.0,1.0),
		      ...){
  ## define helper functions:

  ## JN and Liue et al have this as the prior
  logPriorC <-
    function(etaN,alpha){
      if (alpha==0.0) return( 0.0 )
      n <- sum(etaN)
      length(etaN)*log(alpha) +
        sum(lgamma(etaN)) -
        sum(log(alpha+0:(n-1)))
    }

  ## Antoniak has this as the prior
  logPriorAnt <- function(etaN,alpha){
    m <- table(etaN)
    i <- unique(sort(etaN))
    stopifnot(i>0L)
    n <- sum(i*m)
    num <- lfactorial( n ) + sum(m * log(alpha))
    denom <- 
      sum( m *log(i) ) + sum(lfactorial(m)) +
      sum( log(alpha + 1:n - 1) )
    num - denom
  }


  ## Escobar and West, Page 10  (Chapter in Dey et al, 2012)
  ralpha <- function(alpha,a,b,Istar,I){
    eta <- rbeta(1,alpha+1,I)
    shape <- a + Istar - 1
    rate <- (b - log(eta))
    plus1 <- runif(1, -shape, I*rate) < 0.0
    rgamma(1,shape+plus1, rate)
  }
  

  logpost <- function(wtab,om,eta,etaN,dataToEta,etaM,alpha){
    etaN <- etaN[1:etaM]
    eta <- eta[, 1:etaM, drop = FALSE]
    uniq.indexes <- unique(cbind(wtab$data.index, dataToEta))
    uniq.sums <- table(match(paste(wtab$data.index, dataToEta),
                             paste(uniq.indexes[,1, drop = FALSE],
                                   uniq.indexes[,2, drop = FALSE])))
    marglik <-
      marglogpost(eta[,1L+uniq.indexes[,2], drop = FALSE],
                  om,wtab[["tab"]][uniq.indexes[,1], , drop = FALSE])
    margliksum <- sum( marglik * uniq.sums)
    ## Equations 6 and 10 Jain and Neal, 2007 yield the sum of: 
    c(loglik = margliksum, logprior = logPriorC(etaN,alpha)+etaM*lgamma(nrow(om)))
  }

  mc <- match.call()

  stopifnot(all(om>=0.0))
  if (any(om==0.0)){
    om[om==0.0] <- .Machine[["double.eps"]]
    warning("Converted zeroes in om to machine epsilon")
  }
  stopifnot( length( dpriors ) == 2L )
  dprior <- dpriors[ 1L ]
  dpriorLambda <- dpriors[ 2L ]
  if (is.null(eta))  eta <- array(0.0,c(nrow(om),etaCols))
  if (is.null(dataToEta)) dataToEta <- rep(-1L,length(wtab[["data.index"]]))
  if (is.null(eta))
    eta <- array(0.0,c(nrow(om), etaCols))
  if (is.null(etaN)) etaN <- rep(0L, ncol(eta))
  stopifnot( all( etaM >= dataToEta) )
  stopifnot(length(etaN)==ncol(eta),
            nrow(om) == ncol(wtab[["tab"]]))
  stopifnot(nthin >= 1L)
  stopifnot(all(ab>0), length(ab)==0 || length(ab==2))
  if (ijvals!=0){
    if (etaM!=ijvals) stop("ijvals != etaM (0L, by default) is not permitted")
    if (ncol(eta)!=ijvals) warning("ncol(eta) != ijvals is usually an error")
    if (is.null(dataToEta) || any(dataToEta[1:ijvals]<1) )
      stop("dataToEta[1:ijvals] must be given as positive integers")
  } else {
    if (auxM == 0L) warning("auxM == 0 & ijvals == 0 is usually a mistake")
  }

  if (auxM != 0L && etaM + auxM > ncol(eta)){
    ## need to pad eta and etaN
    padby  <- max(etaM+auxM, etaCols) - ncol(eta)
    eta  <- cbind(eta,array(0.0, c(nrow(eta), padby)))
    etaN <- c(etaN, rep(0.0, padby))
  }

  if (is.null(lambda))
    lambda <- rep(as.double(NA), lambdaSize)
  if (is.null(lambdaN)) lambdaN <- rep(0L, length(lambda))
  if (is.null(dataToLambda)) dataToLambda <-
                               rep(0L,length(wtab$data.index))
  stopifnot( all( lambdaM >= dataToLambda ) )
  if (auxLambdaM != 0L && lambdaM + auxLambdaM > length(lambda)){
    ## need to pad lambda and lambdaN
    padby  <- max(lambdaM+auxLambdaM, lambdaSize) - length(lambda)
    lambda  <- c(lambda,rep(0.0, padby))
    lambdaN <- c(lambdaN, rep(0.0, padby))
  }

  
  tab <- I( wtab[["tab"]] ) # force copy
  di <- wtab[["data.index"]] - 1L
  dataToEta <- dataToEta - 1L
  dataToLambda <- dataToLambda - 1L


  
  ctp <- list(nburn=10L, ncore = 1L, method="sampleX")
  ctp[names(ctParms)] <- ctParms


  N <- sum(wtab[["n"]])
  log.dalpha <- 0.0
  
  pass2 <- list(eta = eta,
                etaN = etaN,
                dataToEta = dataToEta,
                lambda = lambda,
                lambdaN = lambdaN,
                dataToLambda = dataToLambda, 
                etaM = etaM,
                lambdaM = lambdaM)

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
                         auxM = auxM,
                         alpha = alpha,
                         pass1[["lambdaM"]],
                         auxLambdaM,
                         alphaLambda,
                         ijvals = ijvals,
                         verbose = verbose,
                         dprior = dpriors[1])
      
      pass2[c("eta", "lambda")] <-
        sampleParms(tab, di, om,
                    pass2[["dataToEta"]],
                    pass2[["dataToLambda"]],
                    pass2[["eta"]],
                    pass2[["etaM"]],
                    pass2[["lambda"]],
                    pass2[["lambdaN"]],
                    pass2[["lambdaM"]],
                    dprior, dpriorLambda, verbose)
      
      if (!is.null(ab)) {
        alpha <- ralpha(alpha,ab[1],ab[2],pass2[["etaM"]], N)
        log.dalpha <- dgamma(alpha, ab[1], ab[2], log = TRUE)
        
        alphaLambda <- ralpha(alphaLambda,ab[1],ab[2],pass2[["lambdaM"]], N)
        log.dalphaLambda <- dgamma(alphaLambda, ab[1], ab[2], log = TRUE)
      }
    }
    


    
    
    keepers[[ i ]] <-
      with( pass2,
           list(eta = eta[, 1:etaM],
                etaN = etaN[ 1:etaM],
                dataToEta = as.vector(dataToEta + 1L),
                etaM = etaM,
                lambda = lambda[1:lambdaM],
                lambdaN = lambdaN[ 1:lambdaM],
                dataToLambda = as.vector(dataToLambda + 1L),
                lambdaM = lambdaM,
                alpha = alpha,
                alphaLambda = alphaLambda
                ## logpost = logpost(wtab, om, pass2[["eta"]], pass2[["etaN"]],
                ##                   pass2[["dataToEta"]],
                ##                   pass2[["etaM"]],
                ##                   alpha) + log.dalpha
                ))
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
##' @S3method update ctScan
update.ctScan <- function(object, elt = length(object), ...){
  mc <- match.call()
  nameOnly <- c("eta","etaN","dataToEta","etaM","alpha",
                "lambda","lambdaN","dataToLambda","lambdaM","alphaLambda")
  objcall <- attr(object,"call")
  newparms <- list(...)
  newnames <- names(newparms)
  newind <- match(newnames,names(objcall),0L)
  if (any(newnames == 0L)) stop( "couldn't match ", newnames[newind==0L])
  objcall[newnames] <- newparms
  if (is.null(elt)) elt <- length(object)
  elt  <- object[[elt]]
  for ( i in nameOnly) objcall[[i]] <- as.name(i)
  objcall[["update.call"]] <- mc
  eval(objcall,elt, parent.frame())
}

