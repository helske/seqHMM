#include "config.h"

// [[Rcpp::export]]
arma::mat cost_matrix(
    const arma::mat& gamma_pi_est, const arma::mat& gamma_pi_ref, 
    const arma::cube& gamma_A_est, const arma::cube& gamma_A_ref,
    const arma::field<arma::cube>& gamma_B_est, const arma::field<arma::cube>& gamma_B_ref) {
  
  arma::uword S = gamma_A_ref.n_slices;
  arma::uword C = gamma_B_est.n_elem;
  arma::mat costs(S, S);
  arma::uword n_pi = gamma_pi_est.n_cols;
  arma::uword n_A = gamma_A_est.n_rows * gamma_A_est.n_cols;
  arma::uvec n_Bc(C);
  for (arma::uword c = 0; c < C; ++c){
    n_Bc(c) = gamma_B_est(c).n_rows * gamma_B_est(c).n_cols;
  }
  arma::uword n = n_pi + n_A + arma::accu(n_Bc);
  arma::vec est_vec(n);
  arma::vec ref_vec(n);
  arma::uword i;
  for (arma::uword j = 0; j < S; ++j) {
    i = 0;
    ref_vec.subvec(i, i + n_pi - 1) = gamma_pi_ref.row(j).t();
    i += n_pi;
    ref_vec.subvec(i, i + n_A - 1) = arma::vectorise(gamma_A_ref.slice(j));
    i += n_A;
    for (arma::uword c = 0; c < C; ++c) {
      ref_vec.subvec(i, i + n_Bc(c) - 1) = arma::vectorise(gamma_B_ref(c).slice(j));
      i += n_Bc(c);
    }
    for (arma::uword k = 0; k < S; ++k) {
      i = 0;
      est_vec.subvec(i, i + n_pi - 1) = gamma_pi_est.row(k).t();
      i += n_pi;
      est_vec.subvec(i, i + n_A - 1) = arma::vectorise(gamma_A_est.slice(k));
      i += n_A;
      for (arma::uword c = 0; c < C; ++c) {
        est_vec.subvec(i, i + n_Bc(c) - 1) = arma::vectorise(gamma_B_est(c).slice(k));
        i += n_Bc(c);
      }
      costs(k, j) = arma::mean(arma::square(est_vec - ref_vec));
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
      costs(k, j) = arma::mean(arma::square(pcp_est.col(k) - pcp_mle.col(j)));
    }
  }
  return costs;
}
