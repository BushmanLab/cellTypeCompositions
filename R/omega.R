##' Convert Expected Sample Counts to Expected Observation Counts
##'
##' 
##' @title  Expected Observation Counts
##' @param pl list with elements \code{upsilon}, \code{xi}, and \code{psi} 
##' @return matrix of multipliers
##' @author Charles Berry
##' @export
omega <- function(pl) {
    with(pl,diag(upsilon)%*%xi%*%diag(psi))
}
