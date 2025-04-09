#include "config.h"

// [[Rcpp::export]]
arma::mat cost_matrix(
    const arma::mat& gamma_pi_est, const arma::mat& gamma_pi_ref, 
    const arma::cube& gamma_A_est, const arma::cube& gamma_A_ref,
    const arma::field<arma::cube>& gamma_B_est, const arma::field<arma::cube>& gamma_B_ref) {
  
  arma::uword S = gamma_A_ref.n_slices;
  arma::uword C = gamma_B_est.n_elem;
  arma::mat costs(S, S);
  
  for (arma::uword j = 0; j < S; ++j) {
    for (arma::uword k = 0; k < S; ++k) {
      double cost_pi = arma::norm(gamma_pi_est.row(k) - gamma_pi_ref.row(j));
      double cost_A = arma::norm(arma::vectorise(gamma_A_est.slice(k) - gamma_A_ref.slice(j)));
      double cost_B = 0;
      for (arma::uword c = 0; c < C; ++c){
        cost_B += arma::norm(arma::vectorise(gamma_B_est(c).slice(k) - gamma_B_ref(c).slice(j)));
      }
      costs(k, j) = cost_pi + cost_A + cost_B;
    }
  }
  return costs;
}
// [[Rcpp::export]]
arma::mat cost_matrix_clusters(
    const arma::mat& pcp_est, const arma::mat& pcp_mle) {
  
  arma::uword D = pcp_est.n_cols;
  arma::mat costs(D, D);
  
  for (arma::uword j = 0; j < D; ++j) {
    for (arma::uword k = 0; k < D; ++k) {
      costs(k, j) = arma::norm(pcp_est.col(k) - pcp_mle.col(j));
    }
  }
  return costs;
}
