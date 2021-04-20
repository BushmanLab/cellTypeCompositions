// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// auxGibbs
List auxGibbs(arma::imat& tab, arma::ivec& di, arma::mat& om, arma::mat eta, arma::rowvec etaN, arma::ivec diToEta, arma::rowvec lambda, arma::rowvec lambdaN, arma::ivec diToLambda, int etaM, int auxM, double alpha, int lambdaM, int auxLambdaM, double alphaLambda, int ijvals, int verbose, double dprior, double lambdaShape, double lambdaRate);
RcppExport SEXP _cellTypeCompositions_auxGibbs(SEXP tabSEXP, SEXP diSEXP, SEXP omSEXP, SEXP etaSEXP, SEXP etaNSEXP, SEXP diToEtaSEXP, SEXP lambdaSEXP, SEXP lambdaNSEXP, SEXP diToLambdaSEXP, SEXP etaMSEXP, SEXP auxMSEXP, SEXP alphaSEXP, SEXP lambdaMSEXP, SEXP auxLambdaMSEXP, SEXP alphaLambdaSEXP, SEXP ijvalsSEXP, SEXP verboseSEXP, SEXP dpriorSEXP, SEXP lambdaShapeSEXP, SEXP lambdaRateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type tab(tabSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type di(diSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type om(omSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type etaN(etaNSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type diToEta(diToEtaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type lambdaN(lambdaNSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type diToLambda(diToLambdaSEXP);
    Rcpp::traits::input_parameter< int >::type etaM(etaMSEXP);
    Rcpp::traits::input_parameter< int >::type auxM(auxMSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type lambdaM(lambdaMSEXP);
    Rcpp::traits::input_parameter< int >::type auxLambdaM(auxLambdaMSEXP);
    Rcpp::traits::input_parameter< double >::type alphaLambda(alphaLambdaSEXP);
    Rcpp::traits::input_parameter< int >::type ijvals(ijvalsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type dprior(dpriorSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaShape(lambdaShapeSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaRate(lambdaRateSEXP);
    rcpp_result_gen = Rcpp::wrap(auxGibbs(tab, di, om, eta, etaN, diToEta, lambda, lambdaN, diToLambda, etaM, auxM, alpha, lambdaM, auxLambdaM, alphaLambda, ijvals, verbose, dprior, lambdaShape, lambdaRate));
    return rcpp_result_gen;
END_RCPP
}
// sampleParms
List sampleParms(arma::imat& tab, arma::ivec& di, arma::mat& om, arma::ivec dataToEta, arma::ivec dataToLambda, arma::mat eta, int etaM, arma::rowvec lambda, int lambdaM, double dprior, double lambdaAlpha, double lambdaBeta, int niter, int verbose);
RcppExport SEXP _cellTypeCompositions_sampleParms(SEXP tabSEXP, SEXP diSEXP, SEXP omSEXP, SEXP dataToEtaSEXP, SEXP dataToLambdaSEXP, SEXP etaSEXP, SEXP etaMSEXP, SEXP lambdaSEXP, SEXP lambdaMSEXP, SEXP dpriorSEXP, SEXP lambdaAlphaSEXP, SEXP lambdaBetaSEXP, SEXP niterSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type tab(tabSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type di(diSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type om(omSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type dataToEta(dataToEtaSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type dataToLambda(dataToLambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< int >::type etaM(etaMSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type lambdaM(lambdaMSEXP);
    Rcpp::traits::input_parameter< double >::type dprior(dpriorSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaAlpha(lambdaAlphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaBeta(lambdaBetaSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleParms(tab, di, om, dataToEta, dataToLambda, eta, etaM, lambda, lambdaM, dprior, lambdaAlpha, lambdaBeta, niter, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cellTypeCompositions_auxGibbs", (DL_FUNC) &_cellTypeCompositions_auxGibbs, 20},
    {"_cellTypeCompositions_sampleParms", (DL_FUNC) &_cellTypeCompositions_sampleParms, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_cellTypeCompositions(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
