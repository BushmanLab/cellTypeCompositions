
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

inline void rmultnm(int n, double* prob, int k, int* rn){
    double prsum = 0.0;
    for (int i=0; i<k; i++) prsum += prob[i];
    for (int i=0; i<k; i++) prob[i] /= prsum;
    Rf_rmultinom(n, prob,k,rn);
}



// [[Rcpp::export]]
List sampleParms(
		 List wtab,
		 arma::mat& om,
		 arma::ivec& dataToEta,
		 arma::ivec& dataToLambda,
		 arma::mat& eta, int etaM, 
		 arma::vec& lambda, arma::vec& lambdaN, int lambdaM,
		 double dprior,
		 double dpriorLambda, int verbose=0L){
  if (verbose) Rprintf("starting....\n");
  int J = om.n_rows;
  vec omsum = sum(om, 1L);
  imat tab = wtab["tab"];
  ivec r = sum( tab, 1L); // rowSums
  ivec di = wtab["data.index"];
  di = di - 1L;
  int ndat = di.size();
  dataToEta -= 1L;
  dataToLambda -= 1L;
    if (verbose) Rprintf("inits\n");
// eta.by.ct
  imat eta_by_ct( etaM, J , fill::zeros );
// eta.by.lambda.by.r
  imat eta_by_lambda_by_r(etaM, lambdaM, fill::zeros );

  for (int i = 0L; i<ndat; i++){
    try {
      eta_by_ct.row( dataToEta(i)) = eta_by_ct.row( dataToEta(i) ) + tab.row( di(i) );
      eta_by_lambda_by_r( dataToEta(i) , dataToLambda(i) ) += r( di(i) );
    }
    catch (...) {
      Rprintf("i = %d dataToEta(i)) = %d di(i) = %d dataToLambda(i)=%d\n", i,
	      dataToEta(i), di(i), dataToLambda(i));
      i = ndat;
    }
  }
  if (verbose) Rprintf("eta.by.lambda.by.r\n");
  
// rho.vec 
// rho.tilde 

  vec rhocomp = 1.0 - vectorise(trans( omsum ) * eta.head_cols( etaM ));
  mat rhoTilde = 1.0 - ( rhocomp * lambda.head( lambdaM).t() );
  if (verbose) Rprintf("rhoTilde\n");
// R_minus_r

  mat R_minus_r( etaM, lambdaM );
  for (int i = 0L; i<etaM; i++)
    for (int j = 0L; j<lambdaM; j++)
      R_minus_r(i, j) = (eta_by_lambda_by_r(i,j) == 0L) ? 0.0 :
	Rf_rnbinom(( double) eta_by_lambda_by_r(i,j), rhoTilde(i,j));
  if (verbose) Rprintf("R_minus_r\n");
  
// R_by_eta
  vec R_by_eta = sum(R_minus_r,1L);
  if (verbose) Rprintf("R_by_eta\n");
// dropped
  imat dropped(J, etaM, fill::zeros);
  
  for (int i = 0L; i<etaM; i++)
    if (R_by_eta(i) != 0.0){
      vec eta_out = eta.col(i) % (1.0 - omsum);
      rmultnm(R_by_eta(i),eta_out.memptr(),J,dropped.colptr(i));
    }
    if (verbose) Rprintf("dropped\n");
  
// R
  mat R = R_minus_r + eta_by_lambda_by_r;
    if (verbose) Rprintf("R\n");
// eta.by.ct.all
  imat eta_by_ct_all = eta_by_ct + dropped.t();
    if (verbose) Rprintf("eta_by_ct_all\n");
// new.eta
  for (int i = 0L; i< etaM; i++)
    for (int j = 0L; j < J; j++)
      eta(j,i) = Rf_rgamma((double) eta_by_ct_all(i,j) + dprior, 1); 
  if (verbose) Rprintf("gammas\n");
  eta.head_cols(etaM).each_row() /=  sum( eta.head_cols(etaM), 0L); 
  if (verbose) Rprintf("new_eta\n");

// rho.vec 
  rhocomp = 1.0 - vectorise(trans(omsum) * eta.head_cols(etaM));
// rho.tilde
  rhoTilde = 1.0 - ( rhocomp * lambda.head( lambdaM).t() );
  if (verbose) Rprintf("rhoTilde\n");

  // R_minus_r

  for (int i = 0L; i<etaM; i++)
    for (int j = 0L; j<lambdaM; j++)
      R_minus_r(i, j) = (eta_by_lambda_by_r(i,j) == 0L) ? 0.0 :
	Rf_rnbinom(( double) eta_by_lambda_by_r(i,j), rhoTilde(i,j));
  
// R

  R = R_minus_r + eta_by_lambda_by_r;
// 
// update lambda

  for (int i = 0L; i<lambdaM; i++)
    lambda(i) = Rf_rbeta(1.0 + sum(R.col(i)), 1.0 + lambda(i));
  
//

  return List::create(_["eta"] = eta,
                      _["dataToEta"] = dataToEta,
                      _["etaM"] = etaM,
		      _["lambda"] = lambda,
                      _["lambdaN"] = lambdaN,
                      _["dataToLambda"] = dataToLambda,
                      _["lambdaM"] = lambdaM
		      );
  
}
