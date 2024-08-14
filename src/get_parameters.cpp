// log-likelihood of HMM
#include <RcppArmadillo.h>
#include "logsumexp.h"

// beta_raw is (S - 1) x K (start from, covariates)
// X is N x K (sequences, covariates)
// [[Rcpp::export]]
arma::mat get_pi(const arma::mat& beta_raw, const arma::mat& X) {
  arma::mat beta = arma::join_cols(arma::zeros<arma::rowvec>(beta_raw.n_cols), beta_raw);
  arma::mat Pi = beta * X.t();
  for (unsigned int i = 0; i < Pi.n_cols; i++) {
    Pi.col(i) = softmax(Pi.col(i));
  }
  return Pi;
}
// beta_raw is (S - 1) x K x S (transition to, covariates, transition from)
// X is K x N x T (covariates, sequences, time points)
// [[Rcpp::export]]
arma::field<arma::cube> get_A(const arma::cube& beta_raw, const arma::cube& X) {
  unsigned int S = beta_raw.n_slices;
  unsigned int K = X.n_rows;
  unsigned int N = X.n_cols;
  unsigned int T = X.n_slices;
  arma::cube beta(S, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(arma::zeros<arma::rowvec>(K), beta_raw.slice(i));
  }
  arma::field<arma::cube> A(T); // field of T cubes, each S x S x N
  A.fill(arma::zeros<arma::cube>(S, S, N));
  for (unsigned int t = 0; t < T; t++) { // time
    for (unsigned int i = 0; i < N; i++) { // sequences
      arma::mat Atmp(S, S);
      for (unsigned int j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(beta.slice(j) * X.slice(t).col(i));
      }
      A(t).slice(i) = Atmp.t();
    }
  }
  return A;
}
// beta_raw is (M - 1) x K x S (symbols, covariates, transition from)
// X is K x N x T (covariates, sequences, time points)
// [[Rcpp::export]]
arma::field<arma::cube> get_B(const arma::cube& beta_raw, const arma::cube& X) {
  unsigned int S = beta_raw.n_slices;
  unsigned int M = beta_raw.n_rows + 1;
  unsigned int K = X.n_rows;
  unsigned int N = X.n_cols;
  unsigned int T = X.n_slices;
  arma::cube beta(M, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(arma::zeros<arma::rowvec>(K), beta_raw.slice(i));
  }
  arma::field<arma::cube> B(T); // field of T cubes, each S x M x N
  B.fill(arma::zeros<arma::cube>(S, M, N));
  for (unsigned int t = 0; t < T; t++) { // time
    for (unsigned int i = 0; i < N; i++) { // sequences
      arma::mat Btmp(M, S);
      for (unsigned int j = 0; j < S; j ++) { // from states
        Btmp.col(j) = softmax(beta.slice(j) * X.slice(t).col(i));
      }
      B(t).slice(i) = Btmp.t();
    }
  }
  return B;
}
