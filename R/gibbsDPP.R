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
##'     named inputs.
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
    auxGibbs(wtab, om, eta, etaN, dataToEta, etaM,
             auxM , alpha , verbose )
}
