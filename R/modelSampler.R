##' Sample Modelor Estimate Posteriors from Model Specific Posteriors
##' and a Model Prior
##'
##' See \code{\link{modelSampler}}, for which these functions are helpers. 
##' @title Model Posterior Gibbs Samples
##' 
##' @param model.post matrix of clone and model specific posterior
##'     means, e.g. take the mean antilog of each element of \code{
##'     logposterior( ctSampler( ... )) } to get each column of
##'     \code{model.post}
##' @param model.probs vector of model probabilities
##' @param clone.n number of replications of each set of clone counts,
##'     such as \code{\link{uniTab}(...)$n}
##' @return vector of sampled model probabilities for
##'     \code{modePostSampler} or a vector of posterior expectations
##'     for \code{modelPostEM}
##' @author Charles Berry
##' @export
##' @importFrom stats rgamma
modelPostSampler <- function(model.post,model.probs,clone.n){
    totals <- 
        colSums(
            prop.table(
                sweep(model.post, 2, model.probs, "*"),
                1) * clone.n)
    rg <- rgamma( length(totals), 1+totals)
    prop.table(rg)
  }


##' @rdname modelPostSampler
##' @export
modelPostEM <- function(model.post,model.probs,clone.n){
    prop.table(1.0+
               colSums(
                   prop.table(
                       sweep(model.post,2,model.probs,"*"),1)*clone.n))
}


##' Sampler or Estimate Model Posteriors
##'
##' Given a uniform Dirichlet prior on the models and the expected
##' likelihoods or Bayes Factors for a sample of clones, take a sample
##' from the posterior model probability vector or estimate the mean.
##'
##' The default produces Gibbs samples from the model vector
##' posterior. However,
##' @title Model Posteriors
##' @param modelLike matrix of clone and model specific mean
##'     likelihoods, e.g. take the mean antilog of each element of
##'     \code{ logposterior( ctSampler( ... )) } to get each column of
##'     \code{model.post}
##' @param clone.n the multiplicity of clones having each value of
##'     \code{modelLike}
##' @param modelPrior a vector of \code{ncol(modeLike)} probabilities
##' @param nkeep how many samples to retain
##' @param nthin the spacing of samples retained among those collected
##' @param nburn how many samples to discard before the first retained
##'     sample
##' @param sampfun the function use to perform
##'     sampling. \code{modelPostSampler} and \code{modePostEM} are
##'     built in choices.
##' @return
##' @export
##' @author Charles Berry
modelSampler <- function(modelLike, clone.n=1.0, modelPrior=NULL,
			   nkeep=1L, nthin=1L, nburn=0L, sampfun=modelPostSampler){
    if (is.null(modelPrior)) modelPrior <- prop.table(rep(1,ncol(modelLike)))
    nc <- ncol(modelLike)
    stopifnot( length(clone.n)==1 || length(clone.n) == nrow(modelLike),
	      nc == length(modelPrior),
	      nkeep >= 1L)
    prmat <- matrix(0.0,nrow=nc,ncol=nkeep)
    prworking <- modelPrior
    for (i in seq_len(nburn)) prworking <-
				sampfun(modelLike, prworking, clone.n)
    prmat[,1L] <- prworking <- sampfun(modelLike, prworking, clone.n)
    for (i in seq_len(nkeep-1L)){
      for (j in seq_len(nthin)) prworking <-
				 sampfun(modelLike, prworking, clone.n)
      prmat[, i+1]  <- prworking
    }
    prmat
  }
