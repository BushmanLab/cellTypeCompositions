
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

// assume di, dataTo[ Eta | Lambda ] are zero based indexes

// [[Rcpp::export]]
List sampleParms(
		 arma::imat& tab,
		 arma::ivec& di,
		 arma::mat& om,
		 arma::ivec dataToEta,
		 arma::ivec dataToLambda,
		 arma::mat eta, int etaM, 
		 arma::vec lambda, arma::vec& lambdaN, int lambdaM,
		 double dprior,
		 double dpriorLambda, int verbose=0L){
  if (verbose) Rprintf("starting....\n");
  int J = om.n_rows;
  vec omsum = sum(om, 1L);
  ivec r = sum( tab, 1L); // rowSums
  int ndat = di.size();
    if (verbose) Rprintf("inits\n");
// eta.by.ct
  imat eta_by_ct( etaM, J , fill::zeros );
// eta.by.lambda.by.r
  imat eta_by_lambda_by_r(etaM, lambdaM, fill::zeros );

  for (int i = 0L; i<ndat; i++){
    eta_by_ct.row( dataToEta(i)) = eta_by_ct.row( dataToEta(i) ) + tab.row( di(i) );
    eta_by_lambda_by_r( dataToEta(i) , dataToLambda(i) ) += r( di(i) );
  }
  if (verbose>1L) Rprintf("eta.by.lambda.by.r\n");
  
// rho.vec 
// rho.tilde 

  vec rhocomp = 1.0 - vectorise(trans( omsum ) * eta.head_cols( etaM ));
  mat rhoTilde = 1.0 - ( rhocomp * lambda.head( lambdaM).t() );
  if (verbose>1L) Rprintf("rhoTilde\n");
// R_minus_r

  mat R_minus_r( etaM, lambdaM );
  for (int i = 0L; i<etaM; i++)
    for (int j = 0L; j<lambdaM; j++)
      R_minus_r(i, j) = (eta_by_lambda_by_r(i,j) == 0L) ? 0.0 :
	Rf_rnbinom(( double) eta_by_lambda_by_r(i,j), rhoTilde(i,j));
  if (verbose>1L) Rprintf("R_minus_r\n");
  
// R_by_eta
  vec R_by_eta = sum(R_minus_r,1L);
  if (verbose>1L) Rprintf("R_by_eta\n");
// dropped
  imat dropped(J, etaM, fill::zeros);
  
  for (int i = 0L; i<etaM; i++)
    if (R_by_eta(i) != 0.0){
      vec eta_out = eta.col(i) % (1.0 - omsum);
      rmultnm(R_by_eta(i),eta_out.memptr(),J,dropped.colptr(i));
    }
    if (verbose>1L) Rprintf("dropped\n");
  
// R
  mat R = R_minus_r + eta_by_lambda_by_r;
    if (verbose>1L) Rprintf("R\n");
// eta.by.ct.all
  imat eta_by_ct_all = eta_by_ct + dropped.t();
    if (verbose>1L) Rprintf("eta_by_ct_all\n");
// new.eta
  for (int i = 0L; i< etaM; i++)
    for (int j = 0L; j < J; j++)
      eta(j,i) = Rf_rgamma((double) eta_by_ct_all(i,j) + dprior, 1); 
  if (verbose>1L) Rprintf("gammas\n");
  eta.head_cols(etaM).each_row() /=  sum( eta.head_cols(etaM), 0L); 
  if (verbose>1L) Rprintf("new_eta\n");

// rho.vec 
  rhocomp = 1.0 - vectorise(trans(omsum) * eta.head_cols(etaM));
// rho.tilde
  rhoTilde = 1.0 - ( rhocomp * lambda.head( lambdaM).t() );
  if (verbose>1L) Rprintf("rhoTilde\n");

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
    lambda(i) = Rf_rbeta(dpriorLambda + sum(R.col(i)), 1.0 + lambdaN(i));
  
//

  return List::create(_["eta"] = eta,
		      _["lambda"] = lambda);
  
}
