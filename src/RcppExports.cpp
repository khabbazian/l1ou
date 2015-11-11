// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// add_configuration_score_to_db
void add_configuration_score_to_db(std::string str_key, double value);
RcppExport SEXP l1ou_add_configuration_score_to_db(SEXP str_keySEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type str_key(str_keySEXP);
    Rcpp::traits::input_parameter< double >::type value(valueSEXP);
    add_configuration_score_to_db(str_key, value);
    return R_NilValue;
END_RCPP
}
// get_stored_config_score
Rcpp::List get_stored_config_score();
RcppExport SEXP l1ou_get_stored_config_score() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(get_stored_config_score());
    return __result;
END_RCPP
}
// erase_configuration_score_db
void erase_configuration_score_db();
RcppExport SEXP l1ou_erase_configuration_score_db() {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    erase_configuration_score_db();
    return R_NilValue;
END_RCPP
}
// get_score_of_configuration
Rcpp::List get_score_of_configuration(std::string str_key);
RcppExport SEXP l1ou_get_score_of_configuration(SEXP str_keySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type str_key(str_keySEXP);
    __result = Rcpp::wrap(get_score_of_configuration(str_key));
    return __result;
END_RCPP
}
// cmp_sqrt_OU_covariance
Rcpp::List cmp_sqrt_OU_covariance(Rcpp::NumericMatrix edgeList, int nTips);
RcppExport SEXP l1ou_cmp_sqrt_OU_covariance(SEXP edgeListSEXP, SEXP nTipsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type edgeList(edgeListSEXP);
    Rcpp::traits::input_parameter< int >::type nTips(nTipsSEXP);
    __result = Rcpp::wrap(cmp_sqrt_OU_covariance(edgeList, nTips));
    return __result;
END_RCPP
}
