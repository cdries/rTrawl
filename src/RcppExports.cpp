// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// levy_alpha2nu
arma::mat levy_alpha2nu(arma::mat levy_alpha, double b, double beta_0);
RcppExport SEXP _rTrawl_levy_alpha2nu(SEXP levy_alphaSEXP, SEXP bSEXP, SEXP beta_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type levy_alpha(levy_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type beta_0(beta_0SEXP);
    rcpp_result_gen = Rcpp::wrap(levy_alpha2nu(levy_alpha, b, beta_0));
    return rcpp_result_gen;
END_RCPP
}
// levy_alpha_beta
List levy_alpha_beta(arma::vec p_grid, double T0, double TT);
RcppExport SEXP _rTrawl_levy_alpha_beta(SEXP p_gridSEXP, SEXP T0SEXP, SEXP TTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p_grid(p_gridSEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    rcpp_result_gen = Rcpp::wrap(levy_alpha_beta(p_grid, T0, TT));
    return rcpp_result_gen;
END_RCPP
}
// simulate_trawl_uv
List simulate_trawl_uv(std::string levy_seed, arma::vec levy_par, std::string trawl, arma::vec trawl_par, double T0, double TT, double observed_freq, double b);
RcppExport SEXP _rTrawl_simulate_trawl_uv(SEXP levy_seedSEXP, SEXP levy_parSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP, SEXP T0SEXP, SEXP TTSEXP, SEXP observed_freqSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type levy_seed(levy_seedSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type levy_par(levy_parSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par(trawl_parSEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< double >::type observed_freq(observed_freqSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_trawl_uv(levy_seed, levy_par, trawl, trawl_par, T0, TT, observed_freq, b));
    return rcpp_result_gen;
END_RCPP
}
// simulate_trawl_mv
List simulate_trawl_mv(std::string levy_seed, arma::mat levy_par, List trawl, arma::mat trawl_par, arma::mat design_matrix, double T0, double TT, double observed_freq, arma::vec b);
RcppExport SEXP _rTrawl_simulate_trawl_mv(SEXP levy_seedSEXP, SEXP levy_parSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP, SEXP design_matrixSEXP, SEXP T0SEXP, SEXP TTSEXP, SEXP observed_freqSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type levy_seed(levy_seedSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type levy_par(levy_parSEXP);
    Rcpp::traits::input_parameter< List >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type trawl_par(trawl_parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type design_matrix(design_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< double >::type observed_freq(observed_freqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_trawl_mv(levy_seed, levy_par, trawl, trawl_par, design_matrix, T0, TT, observed_freq, b));
    return rcpp_result_gen;
END_RCPP
}
// survival_GIG
arma::vec survival_GIG(arma::vec unif_seed, double gamma, double delta, double nu, double Tmax, double b, double observed_freq);
RcppExport SEXP _rTrawl_survival_GIG(SEXP unif_seedSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP nuSEXP, SEXP TmaxSEXP, SEXP bSEXP, SEXP observed_freqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type unif_seed(unif_seedSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type Tmax(TmaxSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type observed_freq(observed_freqSEXP);
    rcpp_result_gen = Rcpp::wrap(survival_GIG(unif_seed, gamma, delta, nu, Tmax, b, observed_freq));
    return rcpp_result_gen;
END_RCPP
}
// number_parameters_trawl
int number_parameters_trawl(std::string trawl);
RcppExport SEXP _rTrawl_number_parameters_trawl(SEXP trawlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    rcpp_result_gen = Rcpp::wrap(number_parameters_trawl(trawl));
    return rcpp_result_gen;
END_RCPP
}
// trawl_bounds
List trawl_bounds(std::string trawl);
RcppExport SEXP _rTrawl_trawl_bounds(SEXP trawlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    rcpp_result_gen = Rcpp::wrap(trawl_bounds(trawl));
    return rcpp_result_gen;
END_RCPP
}
// trawl_x0
arma::vec trawl_x0(std::string trawl);
RcppExport SEXP _rTrawl_trawl_x0(SEXP trawlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    rcpp_result_gen = Rcpp::wrap(trawl_x0(trawl));
    return rcpp_result_gen;
END_RCPP
}
// vs_sample
arma::vec vs_sample(arma::vec h, arma::vec x_grid, arma::vec p_grid, double T0, double TT, int multi);
RcppExport SEXP _rTrawl_vs_sample(SEXP hSEXP, SEXP x_gridSEXP, SEXP p_gridSEXP, SEXP T0SEXP, SEXP TTSEXP, SEXP multiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_grid(p_gridSEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< int >::type multi(multiSEXP);
    rcpp_result_gen = Rcpp::wrap(vs_sample(h, x_grid, p_grid, T0, TT, multi));
    return rcpp_result_gen;
END_RCPP
}
// vs_SY
List vs_SY(arma::vec h, std::string trawl, arma::vec trawl_par, double beta_0, arma::mat levy_alpha, bool include_cum1, double b, bool include_b);
RcppExport SEXP _rTrawl_vs_SY(SEXP hSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP, SEXP beta_0SEXP, SEXP levy_alphaSEXP, SEXP include_cum1SEXP, SEXP bSEXP, SEXP include_bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par(trawl_parSEXP);
    Rcpp::traits::input_parameter< double >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type levy_alpha(levy_alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type include_cum1(include_cum1SEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type include_b(include_bSEXP);
    rcpp_result_gen = Rcpp::wrap(vs_SY(h, trawl, trawl_par, beta_0, levy_alpha, include_cum1, b, include_b));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rTrawl_levy_alpha2nu", (DL_FUNC) &_rTrawl_levy_alpha2nu, 3},
    {"_rTrawl_levy_alpha_beta", (DL_FUNC) &_rTrawl_levy_alpha_beta, 3},
    {"_rTrawl_simulate_trawl_uv", (DL_FUNC) &_rTrawl_simulate_trawl_uv, 8},
    {"_rTrawl_simulate_trawl_mv", (DL_FUNC) &_rTrawl_simulate_trawl_mv, 9},
    {"_rTrawl_survival_GIG", (DL_FUNC) &_rTrawl_survival_GIG, 7},
    {"_rTrawl_number_parameters_trawl", (DL_FUNC) &_rTrawl_number_parameters_trawl, 1},
    {"_rTrawl_trawl_bounds", (DL_FUNC) &_rTrawl_trawl_bounds, 1},
    {"_rTrawl_trawl_x0", (DL_FUNC) &_rTrawl_trawl_x0, 1},
    {"_rTrawl_vs_sample", (DL_FUNC) &_rTrawl_vs_sample, 6},
    {"_rTrawl_vs_SY", (DL_FUNC) &_rTrawl_vs_SY, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_rTrawl(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
