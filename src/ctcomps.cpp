
#include <Rcpp.h>
using namespace Rcpp;

/* ctcomps.cpp

   Gibbs Sampler for a Dirichlet multinomial model of cell type
   proportions with vague prior for cell counts

   Author: Charles C. Berry
   Date: 13-03-2019
*/




/* rmultnm assumes that:
    - elts of prob sum to 1.0
    - elts of rn have been initialized to zero
   On exit rn should have a multinomial sample (or zeroes if n==0)
*/

/* do what FixupProb in main/random.c does */

inline void rmultnm(int n, double* prob, int k, int* rn){
    double prsum = 0.0;
    for (int i=0; i<k; i++) prsum += prob[i];
    for (int i=0; i<k; i++) prob[i] /= prsum;
    Rf_rmultinom(n, prob,k,rn);
}

// [[Rcpp::export]]
NumericMatrix samplePI(IntegerVector gi, NumericMatrix om, NumericVector pi0,
		       int nkeep, int nthin = 1L, int nburn = 0L, double dprior = 1.0){
  int nsamps = (nkeep - 1L) * nthin + nburn + 1L;
  int isamp = 0L;
  int nrom = om.nrow();
  int ncom = om.ncol();
  NumericVector rwom(nrom);
  for (int i = 0; i<nrom; i++){
    double xx = 0.0;
    for (int j = 0L; j < ncom; j++) xx += om(i,j);
    rwom(i) = xx;
  }
  int r = std::accumulate(gi.begin(), gi.end(), 0L);

  NumericMatrix gammat(nrom, nkeep);

  for (int i=0L; i < nsamps; i++){
    double rho_obs=std::inner_product(pi0.begin(),pi0.end(),rwom.begin(),0.0);
    int R_minus_r = rnbinom(1L,r,rho_obs)[0L];
    NumericVector pi_miss = pi0 * (1.0 - rwom) / (1-rho_obs);
    IntegerVector dropped_f(nrom);
    rmultnm(R_minus_r, REAL(pi_miss),nrom,INTEGER(dropped_f));
    NumericMatrix pi_given_gi(nrom,ncom); 
    for (int j = 0L; j<ncom; j++) {
      pi_given_gi( _, j ) = om( _, j ) * pi0;
      double colsm = std::accumulate(pi_given_gi( _,j).begin(),pi_given_gi( _, j).end(),0.0);
      pi_given_gi( _, j) = pi_given_gi( _, j)/colsm;
    }
    IntegerVector shuffled_f(nrom);
    for (int j=0L; j<ncom; j++){
      IntegerVector ftmp(nrom);
      rmultnm( gi[j], REAL(pi_given_gi)+j*nrom, nrom, INTEGER(ftmp));
      shuffled_f = shuffled_f + ftmp;
    }
    NumericVector dirichlet_parm(nrom);
    for (int j = 0; j<nrom; j++){ 
      dirichlet_parm[j] = dprior + (double) (shuffled_f[j] + dropped_f[j]);
    }
    NumericVector gamma_vals(nrom);
    for (int j = 0; j<nrom; j++) gamma_vals[j] = Rf_rgamma(dirichlet_parm[j],1.0);

    pi0 = gamma_vals / std::accumulate(gamma_vals.begin(),gamma_vals.end(),0.0);
    if ( (i >= nburn) && ((i-nburn)%nthin == 0L) )  
      gammat(_, isamp++) = pi0 * Rf_rgamma((double)(R_minus_r + r),1.0);
  }

  return gammat;

}
