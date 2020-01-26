
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

/*
  auxGibbs.cpp

  Auxiliary Gibbs Sampler for negative multinomial sampler of cell
  type proportions.

  Author: Charles C. Berry
  Date: 24-01-2020
        10-01-2020
        22-04-2019
        16-060-2019
*/

#define MORECOLS 10L
#define MAXCOLS 1000L
#define MORESIZE 10L
#define MAXSIZE 10L

// uniform dirichlet random numbers

inline vec rdirich( int n, double dprior=1.0 ){
  vec rg( n );
  for (int i = 0; i<n; i++) {
    rg(i)=Rf_rgamma(dprior, 1.0);
  }
  return rg/sum(rg);
}


// log probability of r given p and lambda
double logprobp( int r, double p, double lambda ){
  double logpr;
  if ( ISNA(lambda) ){
    logpr = -log( p );
  }
  else
    {
      if (lambda<DOUBLE_XMIN) lambda = DOUBLE_XMIN;
      logpr =
	((double) r - 1.0) * (log( p ) + log( lambda )) -
	(double) r * log1p( - lambda * (1.0 - p )) + log1p( -lambda );
    }
  return logpr;
}

// vectorized log probability of r given p and lambda

inline rowvec logprob_p( int r, double p, rowvec lambda ){
  rowvec result(lambda.size());
  for (int i = 0L; i<lambda.size(); i++) 
    result[i] = logprobp( r, p, lambda[i] );
  return result;
}    

// log probability vector (sans multinomial coefficient)

inline rowvec logprob(irowvec& tabrow, mat& om, mat& eta,
	       rowvec& etaN, int etaLast, double lambda){
  rowvec logpr(etaLast);
  int J = eta.n_rows;
  int K = om.n_cols;
  for (int rc = 0; rc<etaLast; rc++){
    double rhosum = 0.0;
    int tabsum = 0L;
    double logprc = 0.0;

    for (int k = 0; k<K; k++){
      tabsum+=tabrow(k);
      double rhoelt = 0.0;
      for (int j=0;j<J;j++) rhoelt+=om(j,k)*eta(j,rc);
      logprc+= tabrow(k) * log(rhoelt);
      rhosum+= rhoelt;
    }
    logprc-= (double) tabsum * log(rhosum);
    logprc+= logprobp( tabsum, rhosum, lambda );
    logprc+= log( etaN(rc));
    logpr(rc) = logprc;
  }
  return logpr;
}


// sample one index
inline int newIndex(rowvec logpr){
  double maxlogpr = max(logpr);
  double prcum = 0.0;
  for (int i =0;i<logpr.size(); i++){
    prcum += exp(logpr(i)-maxlogpr);
    logpr(i) = prcum;
  }  
  double ur = Rf_runif(0.0,prcum);
  int index=0L;
  for (; logpr[index] < ur && index<logpr.size(); index++);
  return index;
}

// [[Rcpp::export]]
List auxGibbs(List wtab, arma::mat& om, 
              arma::mat eta,
              arma::rowvec etaN,
              arma::ivec diToEta,
	      arma::rowvec lambda,
	      arma::rowvec lambdaN,
	      arma::ivec diToLambda,
              int etaM = 0L,
              int auxM = 5L, double alpha = 100.0,
	      int lambdaM = 0L,
	      int auxLambdaM = 5L, double alphaLambda = 5.0,
              int ijvals = 0L,
              int verbose = 0L,
	      double dprior=1.0) {
  // we get a list from R
  // pull std::vector<double> from R list
  // this is achieved through an implicit
  // call to Rcpp::as
  
  int etaCols = eta.n_cols;
  int lambdaSize = lambda.size();
  int J = om.n_rows;
  imat tab = wtab["tab"];
  ivec di = wtab["data.index"];
  di = di - 1L;
  int ndat = di.size();

  int decN = 0L;
  int incNnew = 0L;
  int incNold = 0L;

  int decLambdaN = 0L;
  int incLambdaNnew = 0L;
  int incLambdaNold = 0L;

  for (int i=ijvals;
       i<ndat && etaM+auxM <= etaCols && lambdaM+auxLambdaM <= lambdaSize;
       i++){
  
    if (verbose>1L) Rprintf("i = %d\n",i);
    
    int di2e = diToEta[ i ];
    double etaN1; // singletons need one less 
    if (di2e >= 0L && etaN( di2e ) == 1.0){
      etaN1 = 1;
      etaN( di2e ) = alpha/auxM;
    }
    else
      { etaN1 = 0; 
	if (di2e>=0L) etaN( di2e )--;
      }
    
    int di2lam = diToLambda[ i ];
    double lambdaN1; // singletons need one less 
    if (di2lam >= 0L && lambdaN( di2lam ) == 1.0){
      lambdaN1 = 1.0;
      lambdaN( di2lam ) = alphaLambda/auxLambdaM;
    }
    else
      { lambdaN1 = 0; 
	if (di2lam>=0L) lambdaN( di2lam )--;
      }
    
    // sample auxM from prior
    for (int j = 0; j < auxM-etaN1; j++){
      eta.col(j + etaM ) = rdirich(J, dprior);
      etaN( j+etaM ) = alpha/auxM;
    }
    
    // rho and logprob
    // initially use lambdaVal = NA_REAL;
    double lambdaVal = (di2lam < 0 ) ? NA_REAL : lambda( di2lam );
    
    irowvec tr = tab.row(di( i ));
    int newind =
      newIndex(logprob( tr, om, eta, etaN,
			etaM + auxM - (int) etaN1, lambdaVal));
    
    // update-eta
    
    if (etaN1){
      //singleton case
      if (newind == di2e)
	{
	  // retain di2e
	  etaN( di2e ) = 1;
	}
      else if (newind < etaM)
	{
	  //move di2e
	  etaN( newind )++;
	  diToEta( i ) = newind;
	  // shift left
	  etaN.subvec(di2e, etaM-1L) = etaN.subvec(di2e+1L, etaM);
	  eta.cols(di2e, etaM-1L) = eta.cols(di2e+1L, etaM);
	  etaM--;
	  decN++;
	  for (int idi=0; idi<ndat; idi++) 
	    if (diToEta[ idi ] >= di2e) diToEta[ idi ]--;
	  if (verbose > 1L) Rprintf("diToEta=%d\n",diToEta[ i ]);
	}
      else
	{
	  // copy to di2e
	  etaN( di2e ) = 1; incNnew++;
	  eta.col(di2e) = eta.col(newind);
	} 
    }
    else
      // initial run or etaN[ di2e ] >= 2
      {
	if (newind >= etaM)
	  {
	    etaN( etaM ) =1;
	    if (newind>etaM) eta.col(etaM) = eta.col(newind);
	    diToEta( i ) = etaM;
	    etaM++;
	    incNnew++;
	  }
	else
	  {
	    etaN( newind )++;
	    diToEta( i ) = newind;
	    incNold++;
	  }
      
      }
    if (etaM+auxM > etaCols){
      int addCols = MORECOLS;
      if (etaCols + addCols <= MAXCOLS){
        etaCols += addCols;
        eta.resize(J, etaCols); 
        etaN.resize( etaCols );
      } else {
	stop("Cannot resize eta");
	  }
      if (verbose) Rprintf("eta has %d Columns\n", eta.n_cols);
    }


    // update-lambda

    // sample lambdaM from prior
    for (int j = 0; j < auxLambdaM-lambdaN1; j++){
      lambda(j + lambdaM ) = Rf_rbeta(1.0,1.0);
      lambdaN( j+lambdaM ) = alphaLambda/auxLambdaM;
    }

    double rhosum = (double) sum( trans(eta.col(newind))*om );
    newind = newIndex(
		      logprob_p( arma::sum(tr), rhosum, lambda.head(lambdaM+auxLambdaM) ) +
		      log( lambdaN.head( lambdaM + auxLambdaM))
		      );

    if (lambdaN1){
      //singleton case
      if (newind == di2lam)
	{
	  // rlambdain di2e
	  lambdaN( di2lam ) = 1;
	}
      else if (newind < lambdaM)
	{
	  //move di2e
	  lambdaN( newind )++;
	  diToLambda( i ) = newind;
	  // shift left
	  lambdaN.subvec(di2lam, lambdaM-1L) = lambdaN.subvec(di2lam+1L, lambdaM);
	  lambda.subvec(di2e, lambdaM-1L) = lambda.subvec(di2e+1L, lambdaM);
	  lambdaM--;
	  decLambdaN++;
	  for (int idi=0; idi<ndat; idi++) 
	    if (diToLambda[ idi ] >= di2lam) diToLambda[ idi ]--;
	  if (verbose > 1L) Rprintf("diToLambda=%d\n",diToLambda[ i ]);
	}
      else
	{
	  // copy to di2lam
	  lambdaN( di2lam ) = 1; incLambdaNnew++;
	  lambda( di2lam ) = lambda( newind );
	} 
    }
    else
      // initial run or lambdaN[ di2lam ] >= 2
      {
	if (newind >= lambdaM)
	  {
	    lambdaN( lambdaM ) =1;
	    if (newind>lambdaM) lambda(lambdaM) = lambda(newind);
	    diToLambda( i ) = lambdaM;
	    lambdaM++;
	    incLambdaNnew++;
	  }
	else
	  {
	    lambdaN( newind )++;
	    diToLambda( i ) = newind;
	    incLambdaNold++;
	  }
      
      }
    if (lambdaM+auxLambdaM > lambdaSize){
      if (verbose) {
	Rprintf("lambda has %d Elts ", lambda.n_elem);
	Rprintf("lambdaSize = %d lambdaM = %d auxLambdaM = %d\n",
		lambdaSize, lambdaM, auxLambdaM);
      }
      int addSize = MORESIZE;
      if (lambdaSize + addSize <= MAXSIZE){
        lambdaSize += addSize;
        lambda.resize( lambdaSize ); 
        lambdaN.resize( lambdaSize );
      } else {
	stop("Cannot resize lambda");
	  }
      if (verbose) Rprintf("lambda has %d Elts\n", lambda.n_elem);
    }



    
  }
  
  if (etaM+auxM > etaCols)
    warning("etaM+auxM = %d would have exceeded %d (maximum)", etaM+auxM, MAXCOLS);
  
  if (lambdaM+auxLambdaM > lambdaSize)
    warning("lambdaM+auxLambdaM = %d would have exceeded %d (maximum)",
	    lambdaM+auxLambdaM, MAXSIZE);
  
  if (verbose)  {
    Rprintf("delete Eta= %d add = %d use existing = %d ",
	    decN, incNnew, incNold);
    Rprintf("delete Lambda= %d add = %d use existing = %d\n",
	    decLambdaN, incLambdaNnew, incLambdaNold);
  }
  
  // return an R list; this is achieved
  // through an implicit call to Rcpp::wrap
  return List::create(_["eta"] = eta,
                      _["etaN"] = etaN,
                      _["dataToEta"] = diToEta,
                      _["etaM"] = etaM,
		      _["lambda"] = lambda,
                      _["lambdaN"] = lambdaN,
                      _["dataToLambda"] = diToLambda,
                      _["lambdaM"] = lambdaM
		      );
}

