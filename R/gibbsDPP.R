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
##' @param lambda (optional) initialization and workspace for the sampled
##'     abundance parameters
##' @param lambdaN (optional) initialization and workspace for the counts
##'     of sampled parameters
##' @param dataToLambda index from the original data to the value of lambda
##'     sampled 
##' @param etaM number of samples provided in the intialization
##' @param auxEtaM how many auxiliary paramters to draw
##' @param alphaEta how much weight to place on new samples
##' @param lambdaM how many auxiliary paramters to draw
##' @param auxLambdaM how much weight to place on new samples
##' @param alphaLambda how much weight to place on new samples
##' @param verbose whether to report intermediate actions
##' @param etaCols how large to make the workspace (too small a value
##'     will result in termination with an error)
##' @param lambdaSize 
##' ##' @param dprior prior value for draws of eta
##' @param ijvals \code{0L} by default, but for splitting and merging
##'     a value of 2 is needed. Users should not usually change the
##'     default.
##' @param dprior when calling \code{gibbsDPP} same as
##'   \code{dpriors[1]} when calling \code{gibbsScan}
##' @param ... currently unused
##' @return list with elements \code{eta}, \code{etaN}, \code{etaM},
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
                     auxEtaM = 5L,
                     alphaEta = 100.0,
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
                  auxEtaM,
                  alphaEta,
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
##' @param dpriors The default value is \code{c(1.0,1.0)}. A smaller
##'     value for \code{dpriors[1]} causes draws from the prior for
##'     \code{eta} to be overdispersed. A larger value for
##'     \code{dpriors[2]} causes posterior samples of \code{eta} to
##'     shrink towards equality. Non-default value can be used to
##'     influence the behavior of the sampler. The computations of the
##'     posterior do not account for non-default values of
##'     \code{dpriors} and should be ignored when non-default values
##'     are used.
##' 
##' @param keep vector of elements of the result in each iteration to
##'   retain or a function called with \code{pass}, the intermediate results, 
##'
##' @param ... currently unused
##' @param abEta if \code{NULL} treat \code{alphaEta} as a fixed
##'     value. Otherwise the first two elemenbts are rthe shape and
##'     rate parameters of Gamma prior for \code{alphaEta}.
##' @param abLambda akin to \code{abEta}
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
		      etaCols = 500L,
                      lambdaSize = 20L,
		      ijvals = 0L,
                      abEta = c(0.0001,0.0001),
                      abLambda = c(1.0,1.0),
		      verbose = FALSE, dpriors = c(1.0,1.0),
                      keep = TRUE,
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
    rgamma(1,shape+plus1, rate = rate)
  }
  

  logpost <- function(tab,di,om,eta,etaN,dataToEta,etaM,alpha,
                      lambda,lambdaN,dataToLambda,
                      lambdaM,alphaLambda){
    ## multinomial
    etaN <- etaN[1:etaM]
    eta <- eta[, 1:etaM, drop = FALSE]
    lambdaN <- lambdaN[1:lambdaM]
    lambda <- lambda[1:lambdaM]
    ## data by Eta by Lambda sums
    iel <- di * (etaM * lambdaM) + dataToEta * lambdaM + dataToLambda
    uiel <- unique(iel)
    ui <- uiel %/% (etaM * lambdaM)
    ue <- (uiel %/% lambdaM) %% etaM
    ul <- (uiel %% lambdaM)
    uind <- match(iel, uiel)
    iel.sums <- tabulate(uind)
    ## data by Eta sums
    uie <- ui * etaM + ue
    uuie <- unique(uie)
    uuind <- match(uie, uuie)
    uui <- uuie%/%etaM
    uue <-   uuie %% etaM
    ie.sums <- tapply(iel.sums, uuind, FUN=sum)

    marglik <-
      marglogpost(eta[, 1L+uue, drop = FALSE],
                  om,tab[ 1L+uui, , drop = FALSE])
    margliksum <- sum(marglik * ie.sums )
    ## Equations 6 and 10 Jain and Neal, 2007 yield 
    logPriEta <- logPriorC(etaN,alphaEta)+etaM*lgamma(nrow(om))
    logPriLambda <- logPriorC( lambdaN, alphaLambda ) # + lambdaM*log(1.0)
    ## rho
    r <- rowSums(tab)
    p <- rowSums(t(om)) %*% eta
    logprp <- 
        mapply(logprobp, #( r.elt, p.elt, lambda.elt)
               r[ 1L+ui ], p[ 1L+ue ], lambda[ 1L+ul ] )
    logprho <- sum( iel.sums * logprp )
    c(multinom=margliksum, n.obs=logprho,eta=logPriEta, lambda=logPriLambda)
  }

  keep.default <-
    function(pass, log.posterior, i, nkeep, keep){
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
                logpost = log.posterior
                ))[ keep ]
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
  stopifnot(all(abEta>0), length(abEta)==0 || length(abEta==2))
  stopifnot(all(abLambda>0), length(abLambda)==0 || length(abLambda==2))
  if (ijvals!=0){
    if (etaM!=ijvals) stop("ijvals != etaM (0L, by default) is not permitted")
    if (ncol(eta)!=ijvals) warning("ncol(eta) != ijvals is usually an error")
    if (is.null(dataToEta) || any(dataToEta[1:ijvals]<1) )
      stop("dataToEta[1:ijvals] must be given as positive integers")
  } else {
    if (auxLambdaM == 0L) warning("auxLambdaM == 0 & ijvals == 0 is usually a mistake")
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
  if (is.null(dataToLambda)) dataToLambda <-
                               rep(0L,length(wtab$data.index))
  stopifnot( all( lambdaM >= dataToLambda ) )
  if (auxLambdaM != 0L && lambdaM + auxLambdaM > length(lambda)){
    ## need to pad lambda and lambdaN
    padby  <- max(lambdaM+auxLambdaM, lambdaSize) - length(lambda)
    lambda  <- c(lambda,rep(0.0, padby))
    lambdaN <- c(lambdaN, rep(0.0, padby))
  }

  if (is.function(keep)){
    keep.fun <- keep
    keep <- TRUE
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
    
    
    log.posterior <-
      with(pass2,
           c(logpost(
             tab, di, om, eta, etaN, dataToEta, etaM, alphaEta,
             lambda, lambdaN, dataToLambda,
             lambdaM, alphaLambda),
             alpha.eta = log.dalphaEta, alpha.lambda = log.dalphaLambda))
    
    keepers[[ i ]] <-
      keep.fun(pass2,log.posterior,i,nkeep,keep)             
             


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

