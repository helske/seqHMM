// simulate MNHMMs
#include "config.h"
#include "fanhmm.h"
#include "list_to_field.h"

// [[Rcpp::export]]
Rcpp::List Rcpp_simulate_fanhmm(
    const arma::field<arma::umat>& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::field<arma::mat>& X_A,
    const Rcpp::List& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& gamma_pi,
    const arma::cube& gamma_A,
    const arma::field<arma::cube>& gamma_B,
    const arma::vec& prior_y,
    const Rcpp::List& W_X_B,
    const Rcpp::List& W_A, 
    const Rcpp::List& W_B) {
  
  fanhmm model(
      obs, Ti, M, X_pi, X_A, matlist_to_2d_field(X_B), 
      icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, gamma_pi, gamma_A, gamma_B, prior_y, W_X_B
  );

  return model.simulate(matlist_to_2d_field(W_A), matlist_to_3d_field(W_B));
}
