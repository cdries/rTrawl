# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

acf_sample_p <- function(h, x_grid, p_grid, TT, lag_max) {
    .Call('_rTrawl_acf_sample_p', PACKAGE = 'rTrawl', h, x_grid, p_grid, TT, lag_max)
}

acf_sample_dp <- function(h, x_grid, p_grid, T0, TT, lag_max, multi) {
    .Call('_rTrawl_acf_sample_dp', PACKAGE = 'rTrawl', h, x_grid, p_grid, T0, TT, lag_max, multi)
}

acf_trawl_p <- function(h, trawl, trawl_par, lag_max) {
    .Call('_rTrawl_acf_trawl_p', PACKAGE = 'rTrawl', h, trawl, trawl_par, lag_max)
}

acf_trawl_dp <- function(h, trawl, trawl_par, b, lag_max) {
    .Call('_rTrawl_acf_trawl_dp', PACKAGE = 'rTrawl', h, trawl, trawl_par, b, lag_max)
}

acf_BN_V <- function(h, trawl, trawl_par, lag_max) {
    .Call('_rTrawl_acf_BN_V', PACKAGE = 'rTrawl', h, trawl, trawl_par, lag_max)
}

acov <- function(h, x_grid, p_grid, T0, TT, lag_max) {
    .Call('_rTrawl_acov', PACKAGE = 'rTrawl', h, x_grid, p_grid, T0, TT, lag_max)
}

add_processes <- function(x_grid1, p_grid1, x_grid2, p_grid2, T0, TT, observed_freq) {
    .Call('_rTrawl_add_processes', PACKAGE = 'rTrawl', x_grid1, p_grid1, x_grid2, p_grid2, T0, TT, observed_freq)
}

ccf_sample_p <- function(h, x_grid1, p_grid1, x_grid2, p_grid2, TT, lag_max) {
    .Call('_rTrawl_ccf_sample_p', PACKAGE = 'rTrawl', h, x_grid1, p_grid1, x_grid2, p_grid2, TT, lag_max)
}

ccf_sample_dp <- function(h, x_grid1, p_grid1, x_grid2, p_grid2, T0, TT, lag_max, multi) {
    .Call('_rTrawl_ccf_sample_dp', PACKAGE = 'rTrawl', h, x_grid1, p_grid1, x_grid2, p_grid2, T0, TT, lag_max, multi)
}

ccf_trawl_p <- function(h, trawl1, trawl_par1, trawl2, trawl_par2, levy_seed, levy_par, design_matrix, lag_max) {
    .Call('_rTrawl_ccf_trawl_p', PACKAGE = 'rTrawl', h, trawl1, trawl_par1, trawl2, trawl_par2, levy_seed, levy_par, design_matrix, lag_max)
}

ccf_trawl_dp <- function(h, trawl1, trawl_par1, trawl2, trawl_par2, b, levy_seed, levy_par, design_matrix, lag_max) {
    .Call('_rTrawl_ccf_trawl_dp', PACKAGE = 'rTrawl', h, trawl1, trawl_par1, trawl2, trawl_par2, b, levy_seed, levy_par, design_matrix, lag_max)
}

cum_sample <- function(ord, x_grid, p_grid, TT) {
    .Call('_rTrawl_cum_sample', PACKAGE = 'rTrawl', ord, x_grid, p_grid, TT)
}

levy_cum_mv2fit <- function(T0, TT, x_grid, p_grid, trawl, trawl_par, p) {
    .Call('_rTrawl_levy_cum_mv2fit', PACKAGE = 'rTrawl', T0, TT, x_grid, p_grid, trawl, trawl_par, p)
}

levy_alpha_beta <- function(x_grid, p_grid, T0, TT) {
    .Call('_rTrawl_levy_alpha_beta', PACKAGE = 'rTrawl', x_grid, p_grid, T0, TT)
}

levy_alpha2nu <- function(levy_alpha, b, beta_0) {
    .Call('_rTrawl_levy_alpha2nu', PACKAGE = 'rTrawl', levy_alpha, b, beta_0)
}

levy_cum_fit <- function(levy_seed, k1_sample, k2_sample) {
    .Call('_rTrawl_levy_cum_fit', PACKAGE = 'rTrawl', levy_seed, k1_sample, k2_sample)
}

levy_varcovar <- function(levy_seed, levy_par, design_matrix) {
    .Call('_rTrawl_levy_varcovar', PACKAGE = 'rTrawl', levy_seed, levy_par, design_matrix)
}

observe_process <- function(x_grid_latent, p_grid_latent, T0, TT, observed_freq) {
    .Call('_rTrawl_observe_process', PACKAGE = 'rTrawl', x_grid_latent, p_grid_latent, T0, TT, observed_freq)
}

simulate_trawl_uv <- function(levy_seed, levy_par, trawl, trawl_par, T0, TT, observed_freq, b) {
    .Call('_rTrawl_simulate_trawl_uv', PACKAGE = 'rTrawl', levy_seed, levy_par, trawl, trawl_par, T0, TT, observed_freq, b)
}

simulate_trawl_mv_Poisson <- function(levy_par, trawl, trawl_par, design_matrix, T0, TT, observed_freq, b) {
    .Call('_rTrawl_simulate_trawl_mv_Poisson', PACKAGE = 'rTrawl', levy_par, trawl, trawl_par, design_matrix, T0, TT, observed_freq, b)
}

simulate_trawl_mv_negBin <- function(levy_par, trawl, trawl_par, design_matrix, T0, TT, observed_freq, b) {
    .Call('_rTrawl_simulate_trawl_mv_negBin', PACKAGE = 'rTrawl', levy_par, trawl, trawl_par, design_matrix, T0, TT, observed_freq, b)
}

number_parameters_trawl <- function(trawl) {
    .Call('_rTrawl_number_parameters_trawl', PACKAGE = 'rTrawl', trawl)
}

trawl_bounds <- function(trawl) {
    .Call('_rTrawl_trawl_bounds', PACKAGE = 'rTrawl', trawl)
}

trawl_x0 <- function(trawl) {
    .Call('_rTrawl_trawl_x0', PACKAGE = 'rTrawl', trawl)
}

leb_AtA <- function(h, trawl, trawl_par) {
    .Call('_rTrawl_leb_AtA', PACKAGE = 'rTrawl', h, trawl, trawl_par)
}

trawl_function <- function(h, trawl, trawl_par) {
    .Call('_rTrawl_trawl_function', PACKAGE = 'rTrawl', h, trawl, trawl_par)
}

vs_sample <- function(h, x_grid, p_grid, T0, TT, multi) {
    .Call('_rTrawl_vs_sample', PACKAGE = 'rTrawl', h, x_grid, p_grid, T0, TT, multi)
}

vs_C <- function(h, trawl, trawl_par, omega, xi, include_b, eta, include_cum1) {
    .Call('_rTrawl_vs_C', PACKAGE = 'rTrawl', h, trawl, trawl_par, omega, xi, include_b, eta, include_cum1)
}

vs_SY <- function(h, trawl, trawl_par, beta_0, levy_alpha, include_cum1, b, include_b) {
    .Call('_rTrawl_vs_SY', PACKAGE = 'rTrawl', h, trawl, trawl_par, beta_0, levy_alpha, include_cum1, b, include_b)
}

