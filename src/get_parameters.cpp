// log-likelihood of HMM
#include <RcppArmadillo.h>
#include "logsumexp.h"

// beta_raw is (S - 1) x K (start from, covariates)
// X is K x N (covariates, sequences)
// [[Rcpp::export]]
arma::mat get_pi(const arma::mat& beta_raw, const arma::mat& X, const int logspace) {
  arma::mat beta = arma::join_cols(arma::zeros<arma::rowvec>(beta_raw.n_cols), beta_raw);
  arma::mat Pi = beta * X;
  for (unsigned int i = 0; i < Pi.n_cols; i++) {
    Pi.col(i) = softmax(Pi.col(i), logspace);
  }
  return Pi;
}
// beta_raw is (S - 1) x K x S (transition to, covariates, transition from)
// X is K x T x N (covariates, time points, sequences)
// [[Rcpp::export]]
arma::field<arma::cube> get_A(const arma::cube& beta_raw, const arma::cube& X, 
                              const int logspace) {
  unsigned int S = beta_raw.n_slices;
  unsigned int K = X.n_rows;
  unsigned int N = X.n_slices;
  unsigned int T = X.n_cols;
  arma::cube beta(S, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(arma::zeros<arma::rowvec>(K), beta_raw.slice(i));
  }
  arma::field<arma::cube> A(N); // field of N cubes, each S x S x T
  A.fill(arma::zeros<arma::cube>(S, S, T));
  arma::mat Atmp(S, S);
  for (unsigned int i = 0; i < N; i++) { // sequences
    for (unsigned int t = 0; t < T; t++) { // time
      for (unsigned int j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(beta.slice(j) * X.slice(i).col(t), logspace);
      }
      A(i).slice(t) = Atmp.t();
    }
  }
  return A;
}
// beta_raw is (M - 1) x K x S (symbols, covariates, transition from)
// X is K x T x N (covariates, time points, sequences)
// [[Rcpp::export]]
arma::field<arma::cube> get_B(const arma::cube& beta_raw, const arma::cube& X, 
                              const int logspace) {
  unsigned int S = beta_raw.n_slices;
  unsigned int M = beta_raw.n_rows + 1;
  unsigned int K = X.n_rows;
  unsigned int N = X.n_slices;
  unsigned int T = X.n_cols;
  arma::cube beta(M, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(arma::zeros<arma::rowvec>(K), beta_raw.slice(i));
  }
  arma::field<arma::cube> B(N); // field of N cubes, each S x M x T
  B.fill(arma::zeros<arma::cube>(S, M, T));
  for (unsigned int i = 0; i < N; i++) { // sequences
    for (unsigned int t = 0; t < T; t++) { // time
      arma::mat Btmp(M, S);
      for (unsigned int j = 0; j < S; j ++) { // from states
        Btmp.col(j) = softmax(beta.slice(j) * X.slice(i).col(t), logspace);
      }
      B(i).slice(t) = Btmp.t();
    }
  }
  return B;
}
// beta_raw is a vector of raw coefficients
// X is K x T x N (covariates, time points, sequences)
// [[Rcpp::export]]
arma::field<arma::cube> get_multichannel_B(
    const arma::vec& beta_raw, 
    const arma::cube& X, unsigned int S, unsigned int C,
    const arma::uvec& M, const int logspace) {
  unsigned int K = X.n_rows;
  unsigned int N = X.n_slices;
  unsigned int T = X.n_cols;
  
  arma::field<arma::cube> beta(C);
  int idx = 0;
  for (unsigned int c = 0; c < C; c++) {
    beta(c) = arma::zeros<arma::cube>(M(c), K, S);
    for (unsigned int s = 0; s < S; s++) {
      arma::mat mat_slice = beta(c).slice(s);
      for (unsigned int m = 1; m < M(c); m++) {
        beta(c).slice(s).row(m - 1) = beta_raw.subvec(idx, idx + K - 1).t();
        idx += K;
      }
    }
  }
  arma::field<arma::cube> B(N, C); // NxC field of cubes, each S x M_c x T
  
  for (unsigned int i = 0; i < N; i++) { // sequences
    for (unsigned int c = 0; c < C; c++) { // channels  
      B(i, c) = arma::zeros<arma::cube>(S, M(c), T);
      for (unsigned int t = 0; t < T; t++) { // time
        arma::mat Btmp(M(c), S);
        for (unsigned int j = 0; j < S; j ++) { // from states
          Btmp.col(j) = softmax(beta(c).slice(j) * X.slice(i).col(t), logspace);
        }
        B(i, c).slice(t) = Btmp.t();
      }
    }
  }
  return B;
}
