// simulate MNHMMs
#include "config.h"
#include "nhmm.h"

// [[Rcpp::export]]
Rcpp::List Rcpp_simulate_nhmm(
    const arma::ucube& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::cube& X_A,
    const arma::field<arma::cube>& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& gamma_pi,
    const arma::cube& gamma_A,
    const arma::field<arma::cube>& gamma_B) {
  
  nhmm model(
      obs, Ti, M, X_pi, X_A, X_B, 
      icpt_only_pi, icpt_only_A, icpt_only_B,
      iv_A, iv_B, tv_A, tv_B, gamma_pi, gamma_A, gamma_B
  );

  return model.simulate();
}
