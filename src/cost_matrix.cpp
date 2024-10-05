#include<RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat cost_matrix_singlechannel(
    const arma::mat& gamma_pi_est, const arma::mat& gamma_pi_ref, 
    const arma::cube& gamma_A_est, const arma::cube& gamma_A_ref,
    const arma::cube& gamma_B_est, const arma::cube& gamma_B_ref) {
  
  unsigned int S = gamma_A_ref.n_slices;
  arma::mat costs(S, S);
  
  for (unsigned int j = 0; j < S; j++) {
    for (unsigned int k = 0; k < S; k++) {
      double cost_pi = arma::norm(gamma_pi_est.row(j) - gamma_pi_ref.row(k));
      double cost_A = arma::norm(arma::vectorise(gamma_A_est.slice(j) - gamma_A_ref.slice(k)));
      double cost_B = arma::norm(arma::vectorise(gamma_B_est.slice(j) - gamma_B_ref.slice(k)));
      costs(k, j) = cost_pi + cost_A + cost_B;
    }
  }
  return costs.t();
}
// [[Rcpp::export]]
arma::mat cost_matrix_multichannel(
    const arma::mat& gamma_pi_est, const arma::mat& gamma_pi_ref, 
    const arma::cube& gamma_A_est, const arma::cube& gamma_A_ref,
    const arma::field<arma::cube>& gamma_B_est, arma::field<arma::cube>& gamma_B_ref) {
  
  unsigned int S = gamma_A_ref.n_slices;
  unsigned int C = gamma_B_est.n_elem;
  arma::mat costs(S, S);
  
  for (unsigned int j = 0; j < S; j++) {
    for (unsigned int k = 0; k < S; k++) {
      double cost_pi = arma::norm(gamma_pi_est.row(j) - gamma_pi_ref.row(k));
      double cost_A = arma::norm(arma::vectorise(gamma_A_est.slice(j) - gamma_A_ref.slice(k)));
      double cost_B = 0;
      for (unsigned int c = 0; c < C; c++){
        cost_B += arma::norm(arma::vectorise(gamma_B_est(c).slice(j) - gamma_B_ref(c).slice(k)));
      }
      costs(k, j) = cost_pi + cost_A + cost_B;
    }
  }
  return costs.t();
}