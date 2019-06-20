
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

/* ctcomps.cpp

   Gibbs Sampler for a Dirichlet multinomial model of cell type
   proportions with vague prior for cell counts

   Author: Charles C. Berry
   Date: 13-03-2019
         16-06-2018
*/


/* rmultnm assumes that:
    - elts of prob sum to 1.0
    - elts of rn have been initialized to zero
   On exit rn should have a multinomial sample (or zeroes if n==0)
*/

#define NDEBUG


/* do what FixupProb in main/random.c does */

inline void rmultnm(int n, double* prob, int k, int* rn){
    double prsum = 0.0;
    for (int i=0; i<k; i++) prsum += prob[i];
    for (int i=0; i<k; i++) prob[i] /= prsum;
    Rf_rmultinom(n, prob,k,rn);
}

// [[Rcpp::export]]
NumericMatrix samplePI(arma::ivec gi, arma::mat& om, arma::rowvec pi0,
		       int nkeep, int nthin = 1L, int nburn = 0L, double dprior = 1.0){
  int nsamps = (nkeep - 1L) * nthin + nburn + 1L;
  int isamp = 0L;
  int nrom = om.n_rows;
  int ncom = om.n_cols;
  vec rwom = sum(om,1L);
  int r = sum( gi );

  mat gammat(nrom, nkeep);

  for (int i=0L; i < nsamps; i++){
    double rho_obs =  dot( pi0, rwom ); // observation prob
    // sample unseen cells
    int R_minus_r = rnbinom(1L,r,rho_obs)[0L];
    vec pi_miss = trans(pi0) % (1.0 - rwom) / (1-rho_obs);
    ivec dropped_f(nrom);
    rmultnm(R_minus_r, pi_miss.memptr(),nrom,dropped_f.memptr());

    // sample true types of observed cells
    mat pi_given_gi = om.each_col() % trans(pi0); // columnwise multiply
    pi_given_gi.each_row() /= sum( pi_given_gi, 0L); // rowwise divide

    ivec shuffled_f(nrom, fill::zeros);
    for (int j=0L; j<ncom; j++){
      ivec ftmp(nrom);
      rmultnm( gi[j], pi_given_gi.memptr() + j*nrom, nrom, ftmp.memptr());
      shuffled_f = shuffled_f + ftmp;
    }
    
    vec dirichlet_parm =  dprior + conv_to<vec>::from(shuffled_f +  dropped_f);
    
    rowvec gamma_vals(nrom);
    for (int j = 0; j<nrom; j++) gamma_vals[j] = Rf_rgamma(dirichlet_parm[j],1.0);
    
    pi0 = gamma_vals / sum(gamma_vals);
    
    if ( (i >= nburn) && ((i-nburn)%nthin == 0L) )  
      gammat.col(isamp++) = trans(pi0) * Rf_rgamma((double)(R_minus_r + r),1.0);

  }

  return wrap(gammat);

}
