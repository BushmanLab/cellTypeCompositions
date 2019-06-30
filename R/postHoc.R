## post hoc adjustments to gibbsScan objects

etaNparent <- function(object,wtab,plist, elts=NULL){
    prDisc <- function(r,rho){
        1.0 - (1-rho)^r * ((rho) /(1-(1-rho)^2))^(r+1)
    }
    if (is.null(elts)) elts <- seq_along(object)
    om <- exOGS(plist)
    r <- rowSums(wtab[["tab"]])
    di <- wtab[["data.index"]]
    lapply(elts, function(ielt){
        elt <- object[[ielt]]
        d2e <- elt[["dataToEta"]]
        rho <- colSums(crossprod(om, elt[["eta"]]))
        pdiscinv <- 1 / prDisc(r[di], rho[d2e] )
        disctab <- 
            xtabs( pdiscinv ~ d2e)
        as.vector(disctab)
    })
}

etaRparent <- function(object,wtab,plist, elts=NULL){
    if (is.null(elts)) elts <- seq_along(object)
    om <- exOGS(plist)
    r <- rowSums(wtab[["tab"]])
    di <- wtab[["data.index"]]
    lapply(elts, function(ielt){
        elt <- object[[ielt]]
        d2e <- elt[["dataToEta"]]
        rho <- colSums(crossprod(om, elt[["eta"]]))    
        
