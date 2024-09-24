// log-likelihood of HMM
#include <RcppArmadillo.h>
#include "softmax.h"

// gamma_omega_raw is (D - 1) x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_omega(const arma::mat& gamma_raw, const arma::vec& X, const bool logspace) {
  arma::mat beta = arma::join_cols(arma::zeros<arma::rowvec>(gamma_raw.n_cols), gamma_raw);
  return softmax(beta * X, logspace);
}
// [[Rcpp::export]]
arma::mat get_omega_all(const arma::mat& gamma_raw, const arma::mat& X, const bool logspace) {
  arma::mat beta = arma::join_cols(arma::zeros<arma::rowvec>(gamma_raw.n_cols), gamma_raw);
  arma::mat omega(beta.n_rows, X.n_cols);
  for (unsigned int i = 0; i < X.n_cols; i++) {
    omega.col(i) = softmax(beta * X.col(i), logspace);
  }
  return omega;
}
// gamma_raw is (S - 1) x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_pi(const arma::mat& gamma_raw, const arma::vec& X, const bool logspace) {
  arma::mat beta = arma::join_cols(arma::zeros<arma::rowvec>(gamma_raw.n_cols), gamma_raw);
  return softmax(beta * X, logspace);
}
// gamma_raw is (S - 1) x K (start from, covariates)
// X a K x N matrix
// [[Rcpp::export]]
arma::mat get_pi_all(const arma::mat& gamma_raw, const arma::mat& X, const bool logspace) {
  arma::mat beta = arma::join_cols(arma::zeros<arma::rowvec>(gamma_raw.n_cols), gamma_raw);
  arma::mat pi(beta.n_rows, X.n_cols);
  for (unsigned int i = 0; i < X.n_cols; i++) {
    pi.col(i) = softmax(beta * X.col(i), logspace);
  }
  return pi;
}

// gamma_raw is (S - 1) x K x S (transition to, covariates, transition from)
// X is K x T matrix (covariates, time points)
// [[Rcpp::export]]
arma::cube get_A(const arma::cube& gamma_raw, const arma::mat& X, 
                 const bool logspace, const bool tv) {
  unsigned int S = gamma_raw.n_slices;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  arma::cube beta(S, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(arma::zeros<arma::rowvec>(K), gamma_raw.slice(i));
  }
  arma::cube A(S, S, T);
  arma::mat Atmp(S, S);
  if (tv) {
    for (unsigned int t = 0; t < T; t++) { // time
      for (unsigned int j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(beta.slice(j) * X.col(t), logspace);
      }
      A.slice(t) = Atmp.t();
    }
  } else {
    for (unsigned int j = 0; j < S; j ++) { // from states
      Atmp.col(j) = softmax(beta.slice(j) * X.col(0), logspace);
    }
    A.each_slice() = Atmp.t();
  }
  return A;
}
// gamma_raw is (S - 1) x K x S (transition to, covariates, transition from)
// X is K x T x N cube (covariates, time points, sequences)
// [[Rcpp::export]]
arma::field<arma::cube> get_A_all(const arma::cube& gamma_raw, 
                                  const arma::cube& X, const bool logspace, 
                                  const bool tv) {
  unsigned int S = gamma_raw.n_slices;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  unsigned int N = X.n_slices;
  arma::cube beta(S, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(
      arma::zeros<arma::rowvec>(K), gamma_raw.slice(i)
    );
  }
  arma::field<arma::cube> A(N);
  arma::mat Atmp(S, S);
  if (tv) {
    for (unsigned int i = 0; i < N; i++) {
      A(i) = arma::cube(S, S, T);
      for (unsigned int t = 0; t < T; t++) { // time
        for (unsigned int j = 0; j < S; j ++) { // from states
          Atmp.col(j) = softmax(beta.slice(j) * X.slice(i).col(t), logspace);
        }
        A(i).slice(t) = Atmp.t();
      }
    }
  } else {
    for (unsigned int i = 0; i < N; i++) {
      A(i) = arma::cube(S, S, T);
      for (unsigned int j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(beta.slice(j) * X.slice(i).col(0), logspace);
      }
      A(i).each_slice() = Atmp.t();
    }
  }
  return A;
}
// gamma_raw is (M - 1) x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
arma::cube get_B(const arma::cube& gamma_raw, const arma::mat& X, 
                 const bool logspace, const bool add_missing, const bool tv) {
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
  if (tv) {
    for (unsigned int t = 0; t < T; t++) { // time
      for (unsigned int j = 0; j < S; j ++) { // from states
        Btmp.col(j).rows(0, M - 1) = softmax(
          beta.slice(j) * X.col(t), logspace
        );
      }
      B.slice(t) = Btmp.t();
    }
  } else {
    for (unsigned int j = 0; j < S; j ++) { // from states
      Btmp.col(j).rows(0, M - 1) = softmax(
        beta.slice(j) * X.col(0), logspace
      );
    }
    B.each_slice() = Btmp.t();
  }
  return B;
}
// gamma_raw is a a field of (M_c - 1) x K x S cubes
// X is K x T (covariates, time point)
arma::field<arma::cube> get_B(
    const arma::field<arma::cube>& gamma_raw, 
    const arma::mat& X, const arma::uvec& M, 
    const bool logspace, const bool add_missing, const bool tv) {
  unsigned int C = M.n_elem;
  arma::field<arma::cube> B(C); // C field of cubes, each S x M_c x T
  for (unsigned int c = 0; c < C; c++) {
    B(c) = get_B(gamma_raw(c), X, logspace, add_missing, tv);
  }
  return B;
}

// gamma_raw is (M - 1) x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
arma::field<arma::cube> get_B_all(
    const arma::cube& gamma_raw,  const arma::cube& X, const bool logspace, 
    const bool add_missing, const bool tv) {
  unsigned int S = gamma_raw.n_slices;
  unsigned int M = gamma_raw.n_rows + 1;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  unsigned int N = X.n_slices;
  arma::cube beta(M, K, S);
  for (unsigned int i = 0; i < S; i++) {
    beta.slice(i) = arma::join_cols(
      arma::zeros<arma::rowvec>(K), gamma_raw.slice(i)
    );
  }
  arma::field<arma::cube> B(N);
  arma::mat Btmp(M + add_missing, S);
  if (add_missing) {
    Btmp.row(M).fill(1.0 - logspace);
  }
  if (tv) {
    for (unsigned int i = 0; i < N; i++) {
      B(i) = arma::cube(S, M + add_missing, T);
      for (unsigned int t = 0; t < T; t++) { // time
        for (unsigned int j = 0; j < S; j ++) { // from states
          Btmp.col(j).rows(0, M - 1) = softmax(
            beta.slice(j) * X.slice(i).col(t), logspace
          );
        }
        B(i).slice(t) = Btmp.t();
      }
    }
  } else {
    for (unsigned int i = 0; i < N; i++) {
      B(i) = arma::cube(S, M + add_missing, T);
      for (unsigned int j = 0; j < S; j ++) { // from states
        Btmp.col(j).rows(0, M - 1) = softmax(
          beta.slice(j) * X.slice(i).col(0), logspace
        );
      }
      B(i).each_slice() = Btmp.t();
    }
  }
  return B;
}