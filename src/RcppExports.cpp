// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logSumExp_cpp
double logSumExp_cpp(NumericVector x);
RcppExport SEXP _SpatMet_logSumExp_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logSumExp_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// state_dist_cpp
NumericVector state_dist_cpp(double G12, double G21);
RcppExport SEXP _SpatMet_state_dist_cpp(SEXP G12SEXP, SEXP G21SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type G12(G12SEXP);
    Rcpp::traits::input_parameter< double >::type G21(G21SEXP);
    rcpp_result_gen = Rcpp::wrap(state_dist_cpp(G12, G21));
    return rcpp_result_gen;
END_RCPP
}
// GeneralLoglikelihood_cpp
double GeneralLoglikelihood_cpp(NumericMatrix y, NumericVector r, NumericVector s, NumericVector u, NumericMatrix Gamma, NumericMatrix e_it, NumericVector B, int model, NumericMatrix z_it, NumericMatrix z_it2);
RcppExport SEXP _SpatMet_GeneralLoglikelihood_cpp(SEXP ySEXP, SEXP rSEXP, SEXP sSEXP, SEXP uSEXP, SEXP GammaSEXP, SEXP e_itSEXP, SEXP BSEXP, SEXP modelSEXP, SEXP z_itSEXP, SEXP z_it2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type e_it(e_itSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z_it(z_itSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z_it2(z_it2SEXP);
    rcpp_result_gen = Rcpp::wrap(GeneralLoglikelihood_cpp(y, r, s, u, Gamma, e_it, B, model, z_it, z_it2));
    return rcpp_result_gen;
END_RCPP
}
// GeneralLoglikelihood_cpp2
double GeneralLoglikelihood_cpp2(NumericMatrix y, NumericVector r, NumericVector s, NumericVector u, NumericMatrix Gamma, NumericMatrix e_it, NumericVector B, int model, NumericMatrix z_it, NumericMatrix z_it2);
RcppExport SEXP _SpatMet_GeneralLoglikelihood_cpp2(SEXP ySEXP, SEXP rSEXP, SEXP sSEXP, SEXP uSEXP, SEXP GammaSEXP, SEXP e_itSEXP, SEXP BSEXP, SEXP modelSEXP, SEXP z_itSEXP, SEXP z_it2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type e_it(e_itSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z_it(z_itSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z_it2(z_it2SEXP);
    rcpp_result_gen = Rcpp::wrap(GeneralLoglikelihood_cpp2(y, r, s, u, Gamma, e_it, B, model, z_it, z_it2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpatMet_logSumExp_cpp", (DL_FUNC) &_SpatMet_logSumExp_cpp, 1},
    {"_SpatMet_state_dist_cpp", (DL_FUNC) &_SpatMet_state_dist_cpp, 2},
    {"_SpatMet_GeneralLoglikelihood_cpp", (DL_FUNC) &_SpatMet_GeneralLoglikelihood_cpp, 10},
    {"_SpatMet_GeneralLoglikelihood_cpp2", (DL_FUNC) &_SpatMet_GeneralLoglikelihood_cpp2, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpatMet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
