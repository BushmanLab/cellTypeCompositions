#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


#define LOGFACTORIAL(x) 0.0 
// uniform dirichlet random numbers

inline vec rdirich( int n, double dprior=1.0 ){
  vec rg( n );
  for (int i = 0; i<n; i++) {
    rg(i)=Rf_rgamma(dprior, 1.0);
  }
  return rg/sum(rg);
}



// log probability of r given p and lambda

// log likelihood kernel sans constants in data r1!, r2!, ...
inline double llkr( int rsum, int n, double p, double lambda ){
  double llk;
  if (ISNA(lambda)) {
    llk = 0.0;
  }
  else
    {
      if (lambda<DOUBLE_XMIN) lambda = DOUBLE_XMIN;
      llk = (double) rsum * (log(lambda)+log(p)) -
	(double) n*(lambda*p + log1p(-exp(-lambda*p)));
    }
  return llk;
}

// lambda vectorized log probability of r given p and lambda

inline arma::rowvec logprob_p( int r, double p, arma::rowvec lambda ){
  rowvec result(lambda.size());
  for (int i = 0L; i<lambda.size(); i++) 
    result[i] = llkr( r, 1L, p, lambda[i] );
  return result;
}    

// p vectorized log probability of r given p and lambda

 inline arma::rowvec logprob_l( int r, arma::rowvec p,  double lambda ){
  rowvec result(p.size());
  for (int i = 0L; i<p.size(); i++) 
    result[i] = llkr( r, 1L, p[i], lambda );
  return result;
}    

#define LOWLIM 8



inline double log_pi_Y_approx(double y, double lowlim, double rho, double alpha, double beta){
  return log(beta) * alpha - lgamma(alpha) + log(rho)*(y-1.0)
    - LOGFACTORIAL(y) + lgamma( alpha + y ) - log (alpha + y - 1.0)
    - log( beta + rho * (lowlim + 1.0)) * (alpha + y - 1.0);
}
inline double log_pi_ky(int k, double y,  double rho, double alpha, double beta){
  return log(beta)*alpha - lgamma(alpha) +
    log(rho) * y - LOGFACTORIAL(y) + lgamma(alpha+y) -
    log(beta+rho*(((double) k + 1.0))) * (alpha+y);
}


inline double log_piy(double y, double rho, double alpha, double beta){
  double log_pi0 = log_pi_ky( 0, y, rho, alpha, beta);
  double pitotal = 1.0 +
    exp(log_pi_Y_approx(y, (double) LOWLIM - 0.5,  rho, alpha, beta) - log_pi0);
  for (int i = 1; i<LOWLIM; i++){
    pitotal +=
      exp( log_pi_ky( i, y, rho, alpha, beta) - log_pi0 );
  }
  return log_pi0 + log(pitotal);
}

// log probability vector (sans multinomial coefficient)

inline arma::rowvec logprob(arma::irowvec& tabrow, arma::mat& om,
			    arma::mat& eta, arma::rowvec& etaN,
			    int etaLast, double lambda){
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
    logprc+= llkr( tabsum, 1L, rhosum, lambda );
    logprc+= log( etaN(rc));
    logpr(rc) = logprc;
  }
  return logpr;
}


// sample one index
inline int newIndex(arma::rowvec logpr){
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

#define INDEF_INT( x )   -pow((beta+rho*(1+x)),-(alpha+r-1))/(rho*(alpha+r-1))

inline double rkReject (double rho,double alpha,double beta, int r,int z){
  vec kpz(z+1L, fill::zeros);
  double kpztot = 0.0;
  int k;
  for (k = 0L; k < z; k++){
    kpz[k] = pow(beta+rho*(1+k),-(alpha+r));
    kpztot += kpz[k];
  }
  kpz[ z ] = -INDEF_INT(z-0.5);
  kpztot += kpz[z];
  bool  reject = true;
  while( reject ){
    double kx = Rf_runif(0.0,kpztot);
    double kpsum=kpz[ 0L ];
    for (k = 0L; k<z && kpsum < kx; k++)
      kpsum += kpz[k+1L];
    if (k==z){
      double x = Rf_runif(0.0,1.0);
      double kcand =
	(beta/rho + z + 0.5)/pow(x, 1.0/(alpha+r-1.0)) - beta/rho - 0.5;
      kcand = Rf_ftrunc(kcand);
      double v = INDEF_INT(kcand+0.5)-INDEF_INT(kcand-0.5);
      double y = Rf_runif(0.0,v);
      reject = y >= pow(beta+rho*(1+kcand),-(alpha+r));
      k = (int) kcand;
    } else {
      reject = false;

    }
  }
  return (double) k;
}


#define KMAX 5L

inline double rlamGivenR(double rho,double alpha,double beta,int r){

  double indx = rkReject(rho, alpha, beta, r, KMAX);

  double rg = Rf_rgamma( alpha+r, 1.0/(beta+rho*(1+indx)));

  return rg;

}

/*
  update a matrix and associated N vector in place
  Author: Charles C. Berry
  Date: 12-04-2020

*/

#define MORESIZE 10L
#define MAXSIZE 500L


int updateXX( int newind,
	       int i,
	       arma::mat& XX,
	       arma::rowvec& XXN,
	       arma::ivec& diToXX,
	       int XXN1,	// singleton flag
	       int& XXM,
	       int& decXXN,  
	       int& incXXNnew,
	       int& incXXNold,
	       int auxXXM,
	       int verbose)
{

  int ndat = diToXX.n_elem;
  int di2XX = diToXX[ i ];
  int XXSize = XX.n_cols;

  if (XXN1){
    //singleton case
    if (verbose>2L) Rprintf("singleton\n");
    if (newind == di2XX)
      {
	// retain di2XX
	XXN( di2XX ) = 1;
      }
    else if (newind < XXM)
      {
	// use existing element in place of this one and move di2e
	XXN( newind )++;
	diToXX( i ) = newind;
	// shift left
	XXN.shed_col(di2XX);
	XX.shed_col(di2XX);
	if (verbose>2L)
	  Rprintf("XXN.n_cols=%d XX.n_cols=%d\n", XXN.n_cols, XX.n_cols);
	XXM--;
	XXSize--;
	decXXN++;
	for (int idi=0; idi<ndat; idi++) 
	  if (diToXX[ idi ] >= di2XX) diToXX[ idi ]--;
	if (verbose > 2L) Rprintf("diToXX=%d\n",diToXX[ i ]);
      }
    else
      {
	// use new element
	// copy to di2XX
	XXN( di2XX ) = 1; incXXNnew++;
	XX.col( di2XX ) = XX.col( newind );
      } 
  }
  else
    // initial run or XXN[ di2XX ] >= 2
    {
      if (verbose>2L) Rprintf("initial run or N>=2\n");
      if (newind >= XXM)
	{
	  //use new element
	  XXN( XXM ) =1;
	  if (newind>XXM) XX.col(XXM) = XX.col(newind);
	  if (verbose>2L) Rprintf("XX(0,XXM)=%f\n", XX(0,XXM));
	  diToXX( i ) = XXM;
	  XXM++;
	  incXXNnew++;
	  if (verbose>2L) Rprintf("XXM=%d\n", XXM);
	}
      else
	{
	  // use existing element
	  XXN( newind )++;
	  diToXX( i ) = newind;
	  incXXNold++;
	}

    }
  // check size and pad as needed
  if (XXM+auxXXM > XXSize){
    if (verbose) {
      Rprintf("XX has %d Elts ", XX.n_elem);
      Rprintf("XXSize = %d XXM = %d auxXXM = %d\n",
	      XXSize, XXM, auxXXM);
    }
    int addSize = MORESIZE;
    if (XXSize + addSize <= MAXSIZE){
      XXSize += addSize;
      XX.resize( XX.n_rows, XXSize ); 
      XXN.resize( XXSize );
    } else {
      Rcpp::stop("Cannot resize XX");
    }
    if (verbose) Rprintf("XX has %d Elts\n", XX.n_elem);
  }
  return XXSize;
}

/*
  auxGibbs.cpp

  Auxiliary Gibbs Sampler for negative multinomial sampler of cell
  type proportions.

  Author: Charles C. Berry
  Date:
  28-07-2020
  24-01-2020
  10-01-2020
  22-04-2019
  16-06-2019
*/

/* assume
   imat tab = wtab["tab"];
   ivec di = wtab["data.index"];
   di = di - 1L;
   di, dataTo[Eta|Lambda] are zero based
*/

// [[Rcpp::export]]
List auxGibbs(arma::imat& tab, arma::ivec& di, arma::mat& om, 
	      arma::mat eta,
	      arma::rowvec etaN,
	      arma::ivec diToEta,
	      arma::rowvec lambda,
	      arma::rowvec lambdaN,
	      arma::ivec diToLambda,
	      int etaM = 0L,
	      int auxM = 5L, double alpha = 100.0,
	      int lambdaM = 0L,
	      int auxLambdaM = 1L, double alphaLambda = 5.0,
	      int ijvals = 0L,
	      int verbose = 0L,
	      double dprior=1.0,
	      double lambdaShape=1.0,
	      double lambdaRate=0.01) {
  // we get a list from R
  // pull std::vector<double> from R list
  // this is achieved through an implicit
  // call to Rcpp::as

  int etaCols = eta.n_cols;
  int lambdaSize = lambda.size();
  int J = om.n_rows;
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

    if (verbose>1L) Rprintf("i = %d",i);

    int di2e = diToEta[ i ];
    int etaN1; // singletons need one less 
    if (di2e >= 0L && etaN( di2e ) == 1.0){
      etaN1 = 1;
      etaN( di2e ) = alpha/auxM;
    }
    else
      { etaN1 = 0; 
	if (di2e>=0L) etaN( di2e )--;
      }

    int di2lam = diToLambda[ i ];
    int lambdaN1; // singletons need one less 
    if (di2lam >= 0 && lambdaN( di2lam ) == 1.0){
      lambdaN1 = 1;
      lambdaN( di2lam ) = alphaLambda;
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
			etaM + auxM - etaN1, lambdaVal));

    // update-eta

    etaCols = updateXX(  newind, i, eta, etaN, diToEta, etaN1, etaM,
	       decN, incNnew, incNold, auxM, verbose);

    // update-lambda

    // sample lambdaM from posterior
    double rhosum = (double) accu( trans(eta.col(newind))*om );
    int tabsum = arma::sum( tr );
    int lambdaElts = (auxLambdaM==0 || lambdaN1) ? lambdaM : lambdaM + 1;
    arma::rowvec lambdaProbs =
      logprob_p( tabsum, rhosum, lambda.head(lambdaElts) ) +
      log( lambdaN.head( lambdaElts));
    int intFGindx = lambdaN1 ? di2lam : lambdaM;
    if (lambdaN1 || auxLambdaM){
      lambdaProbs(intFGindx) =
	log_piy((double) tabsum , rhosum, lambdaShape, lambdaRate) + log(alphaLambda);
    }

    newind = newIndex(lambdaProbs);

    if (newind == intFGindx)
	lambda(newind) = rlamGivenR( rhosum, lambdaShape, lambdaRate, tabsum);

    lambdaSize = updateXX(  newind, i, lambda, lambdaN, diToLambda, lambdaN1, lambdaM,
	       decLambdaN, incLambdaNnew, incLambdaNold, auxLambdaM, verbose);
    if (verbose>1L) Rprintf(" done\n");
  }
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


inline void rmultnm(int n, double* prob, int k, int* rn){
  double prsum = 0.0;
  for (int i=0; i<k; i++) prsum += prob[i];
  for (int i=0; i<k; i++) prob[i] /= prsum;
  Rf_rmultinom(n, prob,k,rn);
}
inline int rN0GivenN1( int N1, double lambda, double rhoplus ){
  return  Rf_rnbinom( (double) N1, -expm1(-lambda * rhoplus ) );
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
		 arma::rowvec lambda, int lambdaM,
		 double dprior,
		 double lambdaAlpha, double lambdaBeta,
		 int niter = 5L, int verbose=0L){
  if (verbose>1L) Rprintf("starting....\n");
  int J = om.n_rows;
  int ndat = di.size();
  ivec rplus = sum( tab, 1L);
  imat etaLambdaTab(etaM, lambdaM, fill::zeros );
  ivec lambdaPostTab(lambdaM, fill::zeros );
  imat etaPostTab(etaM, J, fill::zeros );

  for (int idat = 0; idat < ndat; idat++){
    etaLambdaTab( dataToEta(idat) , dataToLambda(idat))++;
    lambdaPostTab( dataToLambda(idat) ) += rplus(di(idat));
    etaPostTab.row( dataToEta(idat) ) += tab.row(di(idat));
  }

  for (int iter = 0L; iter<niter; iter++){
    vec lambdaN0( etaM , fill::zeros);

    for (int ieta = 0L; ieta<etaM; ieta++)
      for (int ilambda = 0L; ilambda < lambdaM; ilambda++)
	if (etaLambdaTab(ieta, ilambda) != 0L){
	  lambdaN0(ieta) +=
	    lambda(ilambda) * ((double)
			       rN0GivenN1(etaLambdaTab(ieta,ilambda),
					  lambda(ilambda),
					  sum(trans(eta.col(ieta))*om))
			       + (double) etaLambdaTab(ieta,ilambda));
	}

    mat lambdaRhoComp = trans((1.0 - sum(om,1L)) % eta.head_cols(etaM).each_col());
    lambdaRhoComp.each_col() %= lambdaN0;

    mat etaTab( J, etaM, fill::zeros);
    for (int ieta = 0L; ieta<etaM; ieta++){
      mat XY = diagmat(eta.col(ieta)) * om;
      for (int iY = 0; iY<J; iY++){
	ivec X(J, fill::zeros);
	rmultnm(etaPostTab(ieta,iY), XY.colptr(iY), J, X.memptr());
	etaTab.col(ieta) += conv_to<vec>::from(X);
      }
      for (int iX = 0L; iX<J; iX++){
	etaTab(iX,ieta) += Rf_rpois(lambdaRhoComp(ieta, iX));
	etaTab(iX,ieta) = Rf_rgamma(etaTab(iX,ieta)+dprior, 1.0);
      }
    }

    rowvec etaColSum = sum(etaTab, 0L);
    etaTab.each_row() /= etaColSum;

    eta.head_cols(etaM) = etaTab;

    mat lambdaN1(etaM, lambdaM, fill::zeros );

    for (int ieta = 0L; ieta < etaM; ieta++){
      double rhoplus =  sum(trans(eta.col(ieta))*om);
      for (int ilambda = 0L; ilambda <lambdaM; ilambda++)
	if (etaLambdaTab(ieta,ilambda) != 0L)
	  lambdaN1(ieta, ilambda) =
	    rhoplus * (
		       (double) rN0GivenN1(etaLambdaTab(ieta,ilambda),
					   lambda(ilambda), rhoplus) +
		       (double) etaLambdaTab(ieta,ilambda));
    }

    rowvec lambdaUpdate(lambdaM);
    for (int ilambda = 0L; ilambda<lambdaM; ilambda++)
      lambdaUpdate(ilambda) =
	Rf_rgamma(lambdaAlpha + (double) lambdaPostTab(ilambda),
		  1.0/(lambdaBeta + sum(lambdaN1.col(ilambda))));

    lambda.head(lambdaM) = lambdaUpdate;
  }

  double dev = as_scalar(log(lambda.head(lambdaM)) * lambdaPostTab);
  mat rho = trans(eta.head_cols(etaM))*om;
  for (int ieta = 0L; ieta<etaM; ieta++){
    double rhoplus =  sum(rho.row(ieta));
    for (int ilambda = 0L; ilambda<lambdaM; ilambda++)
      dev -=
	etaLambdaTab(ieta, ilambda) *
	(lambda(ilambda)*rhoplus +
	 log( 1.0 - exp( - lambda(ilambda)*rhoplus)));
    for (int j = 0L; j<J; j++) dev += etaPostTab(ieta,j)*log(rho(ieta,j));
  }

  return List::create(_["eta"] = eta,
		      _["lambda"] = lambda,
		      _["etaM"] = etaM,
		      _["lambdaM"] = lambdaM,
		      _["logLik"] = dev);
}
