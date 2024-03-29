// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// computeLinear
NumericVector computeLinear(double curValue, double target, const NumericVector& x, const NumericVector& w, double boundLinear);
RcppExport SEXP _surveysd_computeLinear(SEXP curValueSEXP, SEXP targetSEXP, SEXP xSEXP, SEXP wSEXP, SEXP boundLinearSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type curValue(curValueSEXP);
    Rcpp::traits::input_parameter< double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type boundLinear(boundLinearSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLinear(curValue, target, x, w, boundLinear));
    return rcpp_result_gen;
END_RCPP
}
// computeLinearG1_old
NumericVector computeLinearG1_old(double curValue, double target, const NumericVector& x, const NumericVector& w, double boundLinear);
RcppExport SEXP _surveysd_computeLinearG1_old(SEXP curValueSEXP, SEXP targetSEXP, SEXP xSEXP, SEXP wSEXP, SEXP boundLinearSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type curValue(curValueSEXP);
    Rcpp::traits::input_parameter< double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type boundLinear(boundLinearSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLinearG1_old(curValue, target, x, w, boundLinear));
    return rcpp_result_gen;
END_RCPP
}
// computeLinearG1
NumericVector computeLinearG1(double curValue, double target, const NumericVector& x, const NumericVector& w, double boundLinear);
RcppExport SEXP _surveysd_computeLinearG1(SEXP curValueSEXP, SEXP targetSEXP, SEXP xSEXP, SEXP wSEXP, SEXP boundLinearSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type curValue(curValueSEXP);
    Rcpp::traits::input_parameter< double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type boundLinear(boundLinearSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLinearG1(curValue, target, x, w, boundLinear));
    return rcpp_result_gen;
END_RCPP
}
// geometric_mean_reference
void geometric_mean_reference(NumericVector& w, const IntegerVector& classes);
RcppExport SEXP _surveysd_geometric_mean_reference(SEXP wSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type classes(classesSEXP);
    geometric_mean_reference(w, classes);
    return R_NilValue;
END_RCPP
}
// geometric_mean
NumericVector geometric_mean(const NumericVector& w, const IntegerVector& classes);
RcppExport SEXP _surveysd_geometric_mean(SEXP wSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type classes(classesSEXP);
    rcpp_result_gen = Rcpp::wrap(geometric_mean(w, classes));
    return rcpp_result_gen;
END_RCPP
}
// arithmetic_mean
NumericVector arithmetic_mean(const NumericVector& w, const IntegerVector& classes);
RcppExport SEXP _surveysd_arithmetic_mean(SEXP wSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type classes(classesSEXP);
    rcpp_result_gen = Rcpp::wrap(arithmetic_mean(w, classes));
    return rcpp_result_gen;
END_RCPP
}
// rollMeanC
NumericVector rollMeanC(NumericVector x, int k, char type);
RcppExport SEXP _surveysd_rollMeanC(SEXP xSEXP, SEXP kSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< char >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(rollMeanC(x, k, type));
    return rcpp_result_gen;
END_RCPP
}
// rollSumC
NumericVector rollSumC(NumericVector x, int k, char type);
RcppExport SEXP _surveysd_rollSumC(SEXP xSEXP, SEXP kSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< char >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(rollSumC(x, k, type));
    return rcpp_result_gen;
END_RCPP
}
// weightedRatio
double weightedRatio(NumericVector x, NumericVector w);
RcppExport SEXP _surveysd_weightedRatio(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(weightedRatio(x, w));
    return rcpp_result_gen;
END_RCPP
}
// weightedSum
double weightedSum(NumericVector x, NumericVector w);
RcppExport SEXP _surveysd_weightedSum(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(weightedSum(x, w));
    return rcpp_result_gen;
END_RCPP
}
// ipf_step_ref
void ipf_step_ref(NumericVector w, IntegerVector classes, NumericVector targets);
RcppExport SEXP _surveysd_ipf_step_ref(SEXP wSEXP, SEXP classesSEXP, SEXP targetsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type classes(classesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type targets(targetsSEXP);
    ipf_step_ref(w, classes, targets);
    return R_NilValue;
END_RCPP
}
// ipf_step
NumericVector ipf_step(NumericVector w, IntegerVector classes, NumericVector targets);
RcppExport SEXP _surveysd_ipf_step(SEXP wSEXP, SEXP classesSEXP, SEXP targetsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type classes(classesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type targets(targetsSEXP);
    rcpp_result_gen = Rcpp::wrap(ipf_step(w, classes, targets));
    return rcpp_result_gen;
END_RCPP
}
// ipf_step_f
NumericVector ipf_step_f(NumericVector w, IntegerVector classes, NumericVector targets);
RcppExport SEXP _surveysd_ipf_step_f(SEXP wSEXP, SEXP classesSEXP, SEXP targetsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type classes(classesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type targets(targetsSEXP);
    rcpp_result_gen = Rcpp::wrap(ipf_step_f(w, classes, targets));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_surveysd_computeLinear", (DL_FUNC) &_surveysd_computeLinear, 5},
    {"_surveysd_computeLinearG1_old", (DL_FUNC) &_surveysd_computeLinearG1_old, 5},
    {"_surveysd_computeLinearG1", (DL_FUNC) &_surveysd_computeLinearG1, 5},
    {"_surveysd_geometric_mean_reference", (DL_FUNC) &_surveysd_geometric_mean_reference, 2},
    {"_surveysd_geometric_mean", (DL_FUNC) &_surveysd_geometric_mean, 2},
    {"_surveysd_arithmetic_mean", (DL_FUNC) &_surveysd_arithmetic_mean, 2},
    {"_surveysd_rollMeanC", (DL_FUNC) &_surveysd_rollMeanC, 3},
    {"_surveysd_rollSumC", (DL_FUNC) &_surveysd_rollSumC, 3},
    {"_surveysd_weightedRatio", (DL_FUNC) &_surveysd_weightedRatio, 2},
    {"_surveysd_weightedSum", (DL_FUNC) &_surveysd_weightedSum, 2},
    {"_surveysd_ipf_step_ref", (DL_FUNC) &_surveysd_ipf_step_ref, 3},
    {"_surveysd_ipf_step", (DL_FUNC) &_surveysd_ipf_step, 3},
    {"_surveysd_ipf_step_f", (DL_FUNC) &_surveysd_ipf_step_f, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_surveysd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
