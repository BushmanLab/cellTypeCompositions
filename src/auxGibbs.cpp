
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
        16-060-2019
*/

#define MORECOLS 10L
#define MAXCOLS 1000L

// uniform dirichlet random numbers

inline vec rdirich( int n ){
  vec rg( n );
  for (int i = 0; i<n; i++) {
    rg(i)=Rf_rgamma(1.0,1.0);
  }
  return rg/sum(rg);
}


// log probability vector 
inline rowvec logprob(irowvec& tabrow, mat& om, mat& eta,
		      rowvec& etaN, int etaLast){
  rowvec logpr(etaLast);
  int J = eta.n_rows;
  int K = om.n_cols;
  for (int rc = 0; rc<etaLast; rc++){
    double rhosum = 0.0;
    int tabsum = 1L;
    double logprc = 0.0;
 
    for (int k = 0; k<K; k++){
      tabsum+=tabrow(k);
      double rhoelt = 0.0;
      for (int j=0;j<J;j++) rhoelt+=om(j,k)*eta(j,rc);
      logprc+= tabrow(k) * log(rhoelt);
      rhosum+= rhoelt;
    }
    logprc-= (double) tabsum * log(rhosum);
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
              int etaM = 0L,
              int auxM = 5L, double alpha = 100.0,
              int ijvals = 0L,
              int verbose = 0L) {
  // we get a list from R
  // pull std::vector<double> from R list
  // this is achieved through an implicit
  // call to Rcpp::as
  
  int etaCols = eta.n_cols;
  int J = om.n_rows;
  imat tab = wtab["tab"];
  ivec di = wtab["data.index"];
  di = di - 1L;
  int ndat = di.size();
  
  int decN = 0L;
  int incNnew = 0L;
  int incNold = 0L;
  for (int i=ijvals; i<ndat && etaM+auxM <= etaCols; i++){
  
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
    
    

    // sample auxM from prior
    for (int j = 0; j < auxM-etaN1; j++){
      eta.col(j + etaM ) = rdirich(J);
      etaN( j+etaM ) = alpha/auxM;
    }
    
    // rho and logprob
    
    irowvec tr = tab.row(di( i ));
    int newind =
      newIndex( logprob( tr, om, eta, etaN,  etaM + auxM - (int) etaN1));
    
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
  }
  
  if (etaM+auxM > etaCols) warning("etaM+auxM = %d would have exceeded %d (maximum)", etaM+auxM, MAXCOLS);
  if (verbose)  Rprintf("delete = %d add = %d use existing = %d\n",
			decN, incNnew, incNold);
  // return an R list; this is achieved
  // through an implicit call to Rcpp::wrap
  return List::create(_["eta"] = eta,
                      _["etaN"] = etaN,
                      _["dataToEta"] = diToEta,
                      _["etaM"] = etaM);
}

