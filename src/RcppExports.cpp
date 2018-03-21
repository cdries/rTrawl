// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// acf_sample_p
arma::vec acf_sample_p(double h, arma::vec x_grid, arma::vec p_grid, double TT, int lag_max);
RcppExport SEXP _rTrawl_acf_sample_p(SEXP hSEXP, SEXP x_gridSEXP, SEXP p_gridSEXP, SEXP TTSEXP, SEXP lag_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_grid(p_gridSEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< int >::type lag_max(lag_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(acf_sample_p(h, x_grid, p_grid, TT, lag_max));
    return rcpp_result_gen;
END_RCPP
}
// acf_sample_dp
arma::vec acf_sample_dp(double h, arma::vec x_grid, arma::vec p_grid, double T0, double TT, int lag_max, int multi);
RcppExport SEXP _rTrawl_acf_sample_dp(SEXP hSEXP, SEXP x_gridSEXP, SEXP p_gridSEXP, SEXP T0SEXP, SEXP TTSEXP, SEXP lag_maxSEXP, SEXP multiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_grid(p_gridSEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< int >::type lag_max(lag_maxSEXP);
    Rcpp::traits::input_parameter< int >::type multi(multiSEXP);
    rcpp_result_gen = Rcpp::wrap(acf_sample_dp(h, x_grid, p_grid, T0, TT, lag_max, multi));
    return rcpp_result_gen;
END_RCPP
}
// acf_trawl_p
arma::vec acf_trawl_p(double h, std::string trawl, arma::vec trawl_par, int lag_max);
RcppExport SEXP _rTrawl_acf_trawl_p(SEXP hSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP, SEXP lag_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par(trawl_parSEXP);
    Rcpp::traits::input_parameter< int >::type lag_max(lag_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(acf_trawl_p(h, trawl, trawl_par, lag_max));
    return rcpp_result_gen;
END_RCPP
}
// acf_trawl_dp
arma::vec acf_trawl_dp(double h, std::string trawl, arma::vec trawl_par, double b, int lag_max);
RcppExport SEXP _rTrawl_acf_trawl_dp(SEXP hSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP, SEXP bSEXP, SEXP lag_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par(trawl_parSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type lag_max(lag_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(acf_trawl_dp(h, trawl, trawl_par, b, lag_max));
    return rcpp_result_gen;
END_RCPP
}
// acf_BN_V
List acf_BN_V(double h, std::string trawl, arma::vec trawl_par, int lag_max);
RcppExport SEXP _rTrawl_acf_BN_V(SEXP hSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP, SEXP lag_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par(trawl_parSEXP);
    Rcpp::traits::input_parameter< int >::type lag_max(lag_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(acf_BN_V(h, trawl, trawl_par, lag_max));
    return rcpp_result_gen;
END_RCPP
}
// ccf_sample_p
arma::vec ccf_sample_p(double h, arma::vec x_grid1, arma::vec p_grid1, arma::vec x_grid2, arma::vec p_grid2, double TT, int lag_max);
RcppExport SEXP _rTrawl_ccf_sample_p(SEXP hSEXP, SEXP x_grid1SEXP, SEXP p_grid1SEXP, SEXP x_grid2SEXP, SEXP p_grid2SEXP, SEXP TTSEXP, SEXP lag_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid1(x_grid1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_grid1(p_grid1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid2(x_grid2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_grid2(p_grid2SEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< int >::type lag_max(lag_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(ccf_sample_p(h, x_grid1, p_grid1, x_grid2, p_grid2, TT, lag_max));
    return rcpp_result_gen;
END_RCPP
}
// ccf_sample_dp
arma::vec ccf_sample_dp(double h, arma::vec x_grid1, arma::vec p_grid1, arma::vec x_grid2, arma::vec p_grid2, double T0, double TT, int lag_max, int multi);
RcppExport SEXP _rTrawl_ccf_sample_dp(SEXP hSEXP, SEXP x_grid1SEXP, SEXP p_grid1SEXP, SEXP x_grid2SEXP, SEXP p_grid2SEXP, SEXP T0SEXP, SEXP TTSEXP, SEXP lag_maxSEXP, SEXP multiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid1(x_grid1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_grid1(p_grid1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid2(x_grid2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_grid2(p_grid2SEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< int >::type lag_max(lag_maxSEXP);
    Rcpp::traits::input_parameter< int >::type multi(multiSEXP);
    rcpp_result_gen = Rcpp::wrap(ccf_sample_dp(h, x_grid1, p_grid1, x_grid2, p_grid2, T0, TT, lag_max, multi));
    return rcpp_result_gen;
END_RCPP
}
// ccf_trawl_p
arma::vec ccf_trawl_p(double h, std::string trawl1, arma::vec trawl_par1, std::string trawl2, arma::vec trawl_par2, std::string levy_seed, arma::mat levy_par, arma::mat design_matrix, int lag_max);
RcppExport SEXP _rTrawl_ccf_trawl_p(SEXP hSEXP, SEXP trawl1SEXP, SEXP trawl_par1SEXP, SEXP trawl2SEXP, SEXP trawl_par2SEXP, SEXP levy_seedSEXP, SEXP levy_parSEXP, SEXP design_matrixSEXP, SEXP lag_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl1(trawl1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par1(trawl_par1SEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl2(trawl2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par2(trawl_par2SEXP);
    Rcpp::traits::input_parameter< std::string >::type levy_seed(levy_seedSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type levy_par(levy_parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type design_matrix(design_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type lag_max(lag_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(ccf_trawl_p(h, trawl1, trawl_par1, trawl2, trawl_par2, levy_seed, levy_par, design_matrix, lag_max));
    return rcpp_result_gen;
END_RCPP
}
// cum_sample
double cum_sample(int ord, arma::vec x_grid, arma::vec p_grid, double TT);
RcppExport SEXP _rTrawl_cum_sample(SEXP ordSEXP, SEXP x_gridSEXP, SEXP p_gridSEXP, SEXP TTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ord(ordSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_grid(p_gridSEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    rcpp_result_gen = Rcpp::wrap(cum_sample(ord, x_grid, p_grid, TT));
    return rcpp_result_gen;
END_RCPP
}
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
// levy_cum_fit
arma::vec levy_cum_fit(std::string levy_seed, double k1_sample, double k2_sample);
RcppExport SEXP _rTrawl_levy_cum_fit(SEXP levy_seedSEXP, SEXP k1_sampleSEXP, SEXP k2_sampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type levy_seed(levy_seedSEXP);
    Rcpp::traits::input_parameter< double >::type k1_sample(k1_sampleSEXP);
    Rcpp::traits::input_parameter< double >::type k2_sample(k2_sampleSEXP);
    rcpp_result_gen = Rcpp::wrap(levy_cum_fit(levy_seed, k1_sample, k2_sample));
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
List simulate_trawl_mv(std::string levy_seed, arma::mat levy_par, List trawl, List trawl_par, arma::mat design_matrix, double T0, double TT, double observed_freq, arma::vec b);
RcppExport SEXP _rTrawl_simulate_trawl_mv(SEXP levy_seedSEXP, SEXP levy_parSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP, SEXP design_matrixSEXP, SEXP T0SEXP, SEXP TTSEXP, SEXP observed_freqSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type levy_seed(levy_seedSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type levy_par(levy_parSEXP);
    Rcpp::traits::input_parameter< List >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< List >::type trawl_par(trawl_parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type design_matrix(design_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type T0(T0SEXP);
    Rcpp::traits::input_parameter< double >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< double >::type observed_freq(observed_freqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_trawl_mv(levy_seed, levy_par, trawl, trawl_par, design_matrix, T0, TT, observed_freq, b));
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
// leb_AtA
arma::vec leb_AtA(arma::vec h, std::string trawl, arma::vec trawl_par);
RcppExport SEXP _rTrawl_leb_AtA(SEXP hSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par(trawl_parSEXP);
    rcpp_result_gen = Rcpp::wrap(leb_AtA(h, trawl, trawl_par));
    return rcpp_result_gen;
END_RCPP
}
// trawl_function
arma::vec trawl_function(arma::vec h, std::string trawl, arma::vec trawl_par);
RcppExport SEXP _rTrawl_trawl_function(SEXP hSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par(trawl_parSEXP);
    rcpp_result_gen = Rcpp::wrap(trawl_function(h, trawl, trawl_par));
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
// vs_C
List vs_C(arma::vec h, std::string trawl, arma::vec trawl_par, double omega, double xi, bool include_b, double eta, bool include_cum1);
RcppExport SEXP _rTrawl_vs_C(SEXP hSEXP, SEXP trawlSEXP, SEXP trawl_parSEXP, SEXP omegaSEXP, SEXP xiSEXP, SEXP include_bSEXP, SEXP etaSEXP, SEXP include_cum1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type trawl(trawlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trawl_par(trawl_parSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< bool >::type include_b(include_bSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< bool >::type include_cum1(include_cum1SEXP);
    rcpp_result_gen = Rcpp::wrap(vs_C(h, trawl, trawl_par, omega, xi, include_b, eta, include_cum1));
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
    {"_rTrawl_acf_sample_p", (DL_FUNC) &_rTrawl_acf_sample_p, 5},
    {"_rTrawl_acf_sample_dp", (DL_FUNC) &_rTrawl_acf_sample_dp, 7},
    {"_rTrawl_acf_trawl_p", (DL_FUNC) &_rTrawl_acf_trawl_p, 4},
    {"_rTrawl_acf_trawl_dp", (DL_FUNC) &_rTrawl_acf_trawl_dp, 5},
    {"_rTrawl_acf_BN_V", (DL_FUNC) &_rTrawl_acf_BN_V, 4},
    {"_rTrawl_ccf_sample_p", (DL_FUNC) &_rTrawl_ccf_sample_p, 7},
    {"_rTrawl_ccf_sample_dp", (DL_FUNC) &_rTrawl_ccf_sample_dp, 9},
    {"_rTrawl_ccf_trawl_p", (DL_FUNC) &_rTrawl_ccf_trawl_p, 9},
    {"_rTrawl_cum_sample", (DL_FUNC) &_rTrawl_cum_sample, 4},
    {"_rTrawl_levy_alpha2nu", (DL_FUNC) &_rTrawl_levy_alpha2nu, 3},
    {"_rTrawl_levy_alpha_beta", (DL_FUNC) &_rTrawl_levy_alpha_beta, 3},
    {"_rTrawl_levy_cum_fit", (DL_FUNC) &_rTrawl_levy_cum_fit, 3},
    {"_rTrawl_simulate_trawl_uv", (DL_FUNC) &_rTrawl_simulate_trawl_uv, 8},
    {"_rTrawl_simulate_trawl_mv", (DL_FUNC) &_rTrawl_simulate_trawl_mv, 9},
    {"_rTrawl_number_parameters_trawl", (DL_FUNC) &_rTrawl_number_parameters_trawl, 1},
    {"_rTrawl_trawl_bounds", (DL_FUNC) &_rTrawl_trawl_bounds, 1},
    {"_rTrawl_trawl_x0", (DL_FUNC) &_rTrawl_trawl_x0, 1},
    {"_rTrawl_leb_AtA", (DL_FUNC) &_rTrawl_leb_AtA, 3},
    {"_rTrawl_trawl_function", (DL_FUNC) &_rTrawl_trawl_function, 3},
    {"_rTrawl_vs_sample", (DL_FUNC) &_rTrawl_vs_sample, 6},
    {"_rTrawl_vs_C", (DL_FUNC) &_rTrawl_vs_C, 8},
    {"_rTrawl_vs_SY", (DL_FUNC) &_rTrawl_vs_SY, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_rTrawl(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
