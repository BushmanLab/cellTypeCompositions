
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/*
  auxGibbs.cpp

  Auxiliary Gibbs Sampler for negative multinomial sampler of cell
  type proportions.

  Author: Charles C. Berry
  Date: 22-04-2019

*/;
NumericVector rdirich( int n ){
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
};


// [[Rcpp::export]]
arma::rowvec logprob(NumericVector tabrow, NumericMatrix om, NumericMatrix eta,
		     NumericVector etaN, int etaLast){

  arma::mat om2 = Rcpp::as<arma::mat>(om);
  arma::mat rho = om2.t() * as<arma::mat>(eta).cols(0L,etaLast);
  arma::rowvec rhoSum = arma::sum( rho, 0L);
  int tabsum = arma::sum(Rcpp::as<arma::vec>(tabrow));
  arma::rowvec logpr = Rcpp::as<arma::rowvec>(tabrow) * log(rho);
  logpr = logpr - (tabsum -1L)*log(rhoSum) +
      log(Rcpp::as<arma::rowvec>(etaN).subvec(0L,etaLast) );
  return logpr;
}

// [[Rcpp::export]]
int newIndex(arma::rowvec logpr){
  arma::rowvec pr = cumsum(exp(logpr));
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
	      int auxM = 5L, double alpha = 100.0) {
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

  NumericMatrix rho( K, etaCols );
  NumericVector rhoSum( etaCols );

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
  for (int i=0; i<ndat; i++){
  
  int di2e = diToEta[ i ];

    if ( di2e >= 0L ){
      etaNpt[ di2e ]--; 
      // Rprintf("etaN=%e\n",      etaNpt[ di2e ]);
      if ( etaNpt[ di2e ] == 0.0 ){	
		decN++;
		if( etaM-- > di2e ){
		  std::copy(etaNpt+1L+di2e, etaNpt+etaM, etaNpt+di2e);
		  std::copy(etapt + J * (1L + di2e ),
			    etapt + J * ( etaM ),
			    etapt + J * di2e );
		  for (int idi=0; idi<ndat; idi++) 
		    if (diToEta[ idi ] >= di2e) diToEta[ idi ]--;
		}
      }
    }
    ;
    
    // sample auxM from prior
    for (int j = 0; j<auxM; j++){
      eta( _, j + etaM) = rdirich(J);
      etaNpt[ j+etaM ] = alpha/auxM;
    }
    
    // rho and logprob

    int newind =
      newIndex(logprob(tab(di[ i ], _),om,eta,etaN,etaM+auxM-1L));

    // update-eta>

    if (newind >= etaM){
      eta(_, etaM) = eta(_,newind);
      etaN[etaM] = 1L;incNnew++;
      diToEta[ i ] = etaM;
      etaM++;
    } else {
      etaN[ newind ]++;incNold++;
      diToEta[ i] = newind;
    }
    // Rprintf("diToEta=%d\n",diToEta[ i ]);
  }

  Rprintf("dec = %d new = %d old = %d\n", decN, incNnew, incNold);
// return an R list; this is achieved
// through an implicit call to Rcpp::wrap
    return List::create(_["eta"] = eta,
			_["etaN"] = etaN,
			_["dataToEta"] = diToEta,
			_["etaM"] = etaM);
}

