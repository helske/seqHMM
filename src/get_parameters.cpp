// log-likelihood of HMM
#include <RcppArmadillo.h>
#include "softmax.h"

// gamma_omega_raw is (D - 1) x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_omega(const arma::mat& gamma_omega_raw, const arma::vec X, const int logspace) {
  arma::mat gamma_omega = arma::join_cols(arma::zeros<arma::rowvec>(gamma_omega_raw.n_cols), gamma_omega_raw);
  return softmax(gamma_omega * X, logspace);
}

// gamma_raw is (S - 1) x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_pi(const arma::mat& gamma_raw, const arma::vec X, const int logspace) {
  arma::mat beta = arma::join_cols(arma::zeros<arma::rowvec>(gamma_raw.n_cols), gamma_raw);
  return softmax(beta * X, logspace);
}
// gamma_raw is (S - 1) x K x S (transition to, covariates, transition from)
// X is K x T matrix (covariates, time points)
// [[Rcpp::export]]
arma::cube get_A(const arma::cube& gamma_raw, const arma::mat& X, 
                 const int logspace) {
  unsigned int S = gamma_raw.n_slices;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  arma::cube beta(S, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(arma::zeros<arma::rowvec>(K), gamma_raw.slice(i));
  }
  arma::cube A(S, S, T);
  arma::mat Atmp(S, S);
  for (unsigned int t = 0; t < T; t++) { // time
    for (unsigned int j = 0; j < S; j ++) { // from states
      Atmp.col(j) = softmax(beta.slice(j) * X.col(t), logspace);
    }
    A.slice(t) = Atmp.t();
  }
  return A;
}
// gamma_raw is (M - 1) x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
arma::cube get_B(const arma::cube& gamma_raw, const arma::mat& X, 
                 const int logspace, const int add_missing) {
  unsigned int S = gamma_raw.n_slices;
  unsigned int M = gamma_raw.n_rows + 1;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  arma::cube beta(M, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(arma::zeros<arma::rowvec>(K), gamma_raw.slice(i));
  }
  arma::cube B(S, M + add_missing, T);
  arma::mat Btmp(M + add_missing, S);
  if (add_missing) {
    Btmp.row(M).fill(1.0 - logspace);
  }
  for (unsigned int t = 0; t < T; t++) { // time
    for (unsigned int j = 0; j < S; j ++) { // from states
      Btmp.col(j).rows(0, M - 1) = softmax(
        beta.slice(j) * X.col(t), logspace
      );
    }
    B.slice(t) = Btmp.t();
  }
  return B;
}
// gamma_raw is a a field of (M_c - 1) x K x S cubes
// X is K x T (covariates, time point)
// [[Rcpp::export]]
arma::field<arma::cube> get_B(
    const arma::field<arma::cube>& gamma_raw, 
    const arma::mat& X, const arma::uvec& M, 
    const int logspace, const int add_missing) {
  unsigned int C = M.n_elem;
  arma::field<arma::cube> B(C); // C field of cubes, each S x M_c x T
  for (unsigned int c = 0; c < C; c++) {
    B(c) = get_B(gamma_raw(c), X, logspace, add_missing);
  }
  return B;
}
