// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ordered
IntegerVector ordered(NumericVector vec);
RcppExport SEXP _StatComp22003_ordered(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(ordered(vec));
    return rcpp_result_gen;
END_RCPP
}
// tree_predict
int tree_predict(List node_info, NumericVector x_test);
RcppExport SEXP _StatComp22003_tree_predict(SEXP node_infoSEXP, SEXP x_testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type node_info(node_infoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_test(x_testSEXP);
    rcpp_result_gen = Rcpp::wrap(tree_predict(node_info, x_test));
    return rcpp_result_gen;
END_RCPP
}
// tree_create
List tree_create(NumericMatrix x, IntegerVector y, int q, int s_min, bool random_split);
RcppExport SEXP _StatComp22003_tree_create(SEXP xSEXP, SEXP ySEXP, SEXP qSEXP, SEXP s_minSEXP, SEXP random_splitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type s_min(s_minSEXP);
    Rcpp::traits::input_parameter< bool >::type random_split(random_splitSEXP);
    rcpp_result_gen = Rcpp::wrap(tree_create(x, y, q, s_min, random_split));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp22003_ordered", (DL_FUNC) &_StatComp22003_ordered, 1},
    {"_StatComp22003_tree_predict", (DL_FUNC) &_StatComp22003_tree_predict, 2},
    {"_StatComp22003_tree_create", (DL_FUNC) &_StatComp22003_tree_create, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp22003(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}