
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

/*
  auxGibbs.cpp

  Auxiliary Gibbs Sampler for negative multinomial sampler of cell
  type proportions.

  Author: Charles C. Berry
  Date: 22-04-2019

*/

// uniform dirichlet random numbers

inline NumericVector rdirich( int n ){
  NumericVector rg( n );
  double rgsum = 0.0;
  for (int i = 0; i<n; i++) {
    rg[i]=Rf_rgamma(1.0,1.0);
    rgsum+=rg[i];
  }
  for (int i = 0; i<n; i++) {
    rg[i]/=rgsum;
  }
  return rg;
}


// log probability vector 
inline rowvec logprob(NumericVector tabrow, NumericMatrix& om, NumericMatrix& eta,
			    NumericVector& etaN, int etaLast, mat& rho){

  mat om2 = Rcpp::as<mat>(om);
  rho.cols(0L,etaLast) = om2.t() * as<mat>(eta).cols(0L,etaLast);
  rowvec rhoSum = sum( rho.cols(0L, etaLast), 0L);
  int tabsum = sum(Rcpp::as<vec>(tabrow));
  rowvec logpr = Rcpp::as<rowvec>(tabrow) * log(rho.cols(0L,etaLast));
  logpr = logpr - (tabsum -1L)*log(rhoSum) +
      log(Rcpp::as<rowvec>(etaN).subvec(0L,etaLast) );
  return logpr;
}

// sample one index
inline int newIndex(rowvec logpr){
  rowvec pr = cumsum(exp(logpr));
  double prsum = as_scalar(pr.tail(1L));
  double ur = Rf_runif(0.0,prsum);
  int index=0L;
  for (; pr[index] < ur && index<pr.size(); index++);
  return index;
}

// [[Rcpp::export]]
List auxGibbs(List wtab, NumericMatrix om, 
	      NumericMatrix eta_orig,
	      NumericVector etaN_orig,
	      IntegerVector diToEta_orig,
	      int etaM = 0L,
	      int auxM = 5L, double alpha = 100.0,
	      int ijvals = 0L,
	      int verbose = 0L) {
  // we get a list from R
  // pull std::vector<double> from R list
  // this is achieved through an implicit
  // call to Rcpp::as

  NumericMatrix eta(clone(eta_orig));
  NumericVector etaN(clone(etaN_orig));
  IntegerVector diToEta(clone(diToEta_orig));

  int etaCols = eta.ncol();
  int J = om.nrow();
  int K = om.ncol();
  mat rho(K, etaCols);
  
  
  // NumericMatrix rho( K, etaCols );
  // NumericVector rhoSum( etaCols );

  double *etaNpt = REAL(etaN);
  double *etapt = REAL(eta);

  NumericMatrix tab = wtab["tab"];
  IntegerVector di_orig = wtab["data.index"];
  IntegerVector di(clone(di_orig));
  di = di - 1L;
  int ndat = di.size();

  int decN = 0L;
  int incNnew = 0L;
  int incNold = 0L;
  for (int i=ijvals; i<ndat; i++){
  
    int di2e = diToEta[ i ];


    double etaN1; // singletons need one less 
    if (di2e >= 0L && etaNpt[ di2e] == 1.0){
      etaN1 = 1;
      etaNpt[ di2e ] = alpha/auxM;
    }
    else
      { etaN1 = 0; 
	if (di2e>=0L) etaNpt[ di2e ]--;
      }

    
    // sample auxM from prior
    for (int j = 0; j < auxM-etaN1; j++){
      eta( _, j + etaM) = rdirich(J);
      etaNpt[ j+etaM ] = alpha/auxM;
    }
    
    // rho and logprob

    int newind =
      newIndex(logprob(tab(di[ i ], _),om,eta,etaN,etaM+auxM-etaN1-1L,rho));

    // update-eta

    if (etaN1){
      //singleton case
      if (newind == di2e)
	{
	  // retain di2e
	  etaNpt[ di2e ] = 1;
	}
      else if (newind < etaM)
	{
	  //move di2e
	  etaNpt[ newind ]++;
	  diToEta[ i ] = newind;
	  // shift left
	  std::copy(etaNpt+1L+di2e, etaNpt+etaM, etaNpt+di2e);
	  std::copy(etapt + J * (1L + di2e ),
		    etapt + J * ( etaM ),
		    etapt + J * di2e );
	  etaM--;
	  decN++;
	  for (int idi=0; idi<ndat; idi++) 
	    if (diToEta[ idi ] >= di2e) diToEta[ idi ]--;
	  // Rprintf("diToEta=%d\n",diToEta[ i ]);
	}
      else
	{
	  // copy to di2e
	  etaNpt[ di2e ] = 1;incNnew++;
	  eta(_,di2e) = eta(_,newind);
	} 
    }
    else
      // initial run or etaN[ di2e ] >= 2
      {
	if (newind >= etaM)
	  {
	    etaNpt[ etaM ] =1;
	    if (newind>etaM) eta(_,etaM) = eta(_,newind);
	    diToEta[ i ] = etaM;
	    etaM++;
	    incNnew++;
	  }
	else
	  {
	    etaNpt[ newind ]++;
	    diToEta[ i ] = newind;
	    incNold++;
	  }
	
      }
  }

  if (verbose)  Rprintf("delete = %d add = %d use existing = %d\n",
			decN, incNnew, incNold);
  // return an R list; this is achieved
  // through an implicit call to Rcpp::wrap
  return List::create(_["eta"] = eta,
		      _["etaN"] = etaN,
		      _["dataToEta"] = diToEta,
		      _["etaM"] = etaM);
}

