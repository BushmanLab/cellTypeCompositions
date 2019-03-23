##' Sample Model Posteriors from Model Specific Posteriors and a Model Prior
##'
##' Given a uniform Dirichlet prior on the
##' @title Model Posterior Samples
##' @param model.post matrix of clone and model specific posterior
##'     means, e.g. take the mean antilog of each element of \code{
##'     logposterior( ctSampler( ... )) } to get each column of
##'     \code{model.post}
##' @param model.probs vector of model probabilities
##' @param clone.n number of replications of each set of clone counts,
##'     such as \code{\link{uniTab}(...)$n}
##' @return vector of sampled model probabilities
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
