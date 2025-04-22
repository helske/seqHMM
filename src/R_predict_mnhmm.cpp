// forward predictions for NHMM
#include "config.h"
#include "mnhmm.h"
#include "list_to_field.h"

// [[Rcpp::export]]
Rcpp::List Rcpp_predict_mnhmm(
    const arma::field<arma::umat>& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::field<arma::mat>& X_A,
    const Rcpp::List& X_B,
    const arma::mat& X_omega,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool icpt_only_omega,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::field<arma::mat>& gamma_pi,
    const arma::field<arma::cube>& gamma_A,
    const Rcpp::List& gamma_B,
    const arma::mat& gamma_omega) {
  
  mnhmm model(
      obs, Ti, M, X_pi, X_A, matlist_to_2d_field(X_B), X_omega, 
      icpt_only_pi, icpt_only_A, icpt_only_B, icpt_only_omega,
      iv_A, iv_B, tv_A, tv_B, gamma_pi, gamma_A, cubelist_to_2d_field(gamma_B), 
      gamma_omega
  );
  return model.predict();
}
