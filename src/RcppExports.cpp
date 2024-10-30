// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_Laplace
List Rcpp_Laplace(NumericVector lat0, double v, int K, Function PDcorrect, Function Dloglik, Function D2loglik, Function Qv);
RcppExport SEXP _mixcurelps_Rcpp_Laplace(SEXP lat0SEXP, SEXP vSEXP, SEXP KSEXP, SEXP PDcorrectSEXP, SEXP DloglikSEXP, SEXP D2loglikSEXP, SEXP QvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lat0(lat0SEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Function >::type PDcorrect(PDcorrectSEXP);
    Rcpp::traits::input_parameter< Function >::type Dloglik(DloglikSEXP);
    Rcpp::traits::input_parameter< Function >::type D2loglik(D2loglikSEXP);
    Rcpp::traits::input_parameter< Function >::type Qv(QvSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_Laplace(lat0, v, K, PDcorrect, Dloglik, D2loglik, Qv));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_Laplace2
List Rcpp_Laplace2(NumericVector lat0, double v, int K, Function Dloglik, Function D2loglik, Function Qv);
RcppExport SEXP _mixcurelps_Rcpp_Laplace2(SEXP lat0SEXP, SEXP vSEXP, SEXP KSEXP, SEXP DloglikSEXP, SEXP D2loglikSEXP, SEXP QvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lat0(lat0SEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Function >::type Dloglik(DloglikSEXP);
    Rcpp::traits::input_parameter< Function >::type D2loglik(D2loglikSEXP);
    Rcpp::traits::input_parameter< Function >::type Qv(QvSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_Laplace2(lat0, v, K, Dloglik, D2loglik, Qv));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_cubicBspline
NumericMatrix Rcpp_cubicBspline(NumericVector x, double lower, double upper, int K);
RcppExport SEXP _mixcurelps_Rcpp_cubicBspline(SEXP xSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_cubicBspline(x, lower, upper, K));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mixcurelps_Rcpp_Laplace", (DL_FUNC) &_mixcurelps_Rcpp_Laplace, 7},
    {"_mixcurelps_Rcpp_Laplace2", (DL_FUNC) &_mixcurelps_Rcpp_Laplace2, 6},
    {"_mixcurelps_Rcpp_cubicBspline", (DL_FUNC) &_mixcurelps_Rcpp_cubicBspline, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mixcurelps(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
