#include "get_parameters.h"
#include "sum_to_zero.h"

// eta_omega is (D - 1) x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_omega(const arma::mat& gamma, const arma::vec& X) {
  return softmax(gamma * X);
}
// eta_omega is (D - 1) x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_log_omega(const arma::mat& gamma, const arma::vec& X) {
  return arma::log(softmax(gamma * X));
}
// [[Rcpp::export]]
arma::mat get_omega_all(const arma::mat& gamma, const arma::mat& X) {
  arma::mat omega(gamma.n_rows, X.n_cols);
  for (unsigned int i = 0; i < X.n_cols; i++) {
    omega.col(i) = softmax(gamma * X.col(i));
  }
  return omega;
}

// gamma is (S - 1) x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_pi(const arma::mat& gamma, const arma::vec& X) {
  return softmax(gamma * X);
}
// gamma is (S - 1) x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_log_pi(const arma::mat& gamma, const arma::vec& X) {
  return arma::log(softmax(gamma * X));
}
// gamma is (S - 1) x K (start from, covariates)
// X a K x N matrix
// [[Rcpp::export]]
arma::mat get_pi_all(const arma::mat& gamma, const arma::mat& X) {
  arma::mat pi(gamma.n_rows, X.n_cols);
  for (unsigned int i = 0; i < X.n_cols; i++) {
    pi.col(i) = softmax(gamma * X.col(i));
  }
  return pi;
}
// gamma is (S - 1) x K x B (start from, covariates, sample)
// X a K x N matrix
// [[Rcpp::export]]
arma::cube get_pi_boot(const arma::cube& gamma, const arma::mat& X) {
  unsigned int B = gamma.n_slices;
  arma::cube pi(gamma.n_rows, X.n_cols, B);
  for (unsigned int b = 0; b < B; b++) {
    for (unsigned int i = 0; i < X.n_cols; i++) {
      pi.slice(b).col(i) = softmax(gamma.slice(b) * X.col(i));
    }
  }
  return pi;
}
// gamma is (S - 1) x K x S (transition to, covariates, transition from)
// X is K x T matrix (covariates, time points)
// [[Rcpp::export]]
arma::cube get_A(const arma::cube& gamma, const arma::mat& X, 
                 const bool tv) {
  unsigned int S = gamma.n_slices;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  arma::cube A(S, S, T);
  arma::mat Atmp(S, S);
  if (tv) {
    for (unsigned int t = 0; t < T; t++) { // time
      for (unsigned int j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(gamma.slice(j) * X.col(t));
      }
      A.slice(t) = Atmp.t();
    }
  } else {
    for (unsigned int j = 0; j < S; j ++) { // from states
      Atmp.col(j) = softmax(gamma.slice(j) * X.col(0));
    }
    A.each_slice() = Atmp.t();
  }
  return A;
}
// gamma is (S - 1) x K x S (transition to, covariates, transition from)
// X is K x T matrix (covariates, time points)
// [[Rcpp::export]]
arma::cube get_log_A(const arma::cube& gamma, const arma::mat& X, 
                     const bool tv) {
  unsigned int S = gamma.n_slices;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  arma::cube A(S, S, T);
  arma::mat Atmp(S, S);
  if (tv) {
    for (unsigned int t = 0; t < T; t++) { // time
      for (unsigned int j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(gamma.slice(j) * X.col(t));
      }
      A.slice(t) = Atmp.t();
    }
  } else {
    for (unsigned int j = 0; j < S; j ++) { // from states
      Atmp.col(j) = softmax(gamma.slice(j) * X.col(0));
    }
    A.each_slice() = Atmp.t();
  }
  return arma::log(A);
}
// gamma is (S - 1) x K x S (transition to, covariates, transition from)
// X is K x T x N cube (covariates, time points, sequences)
// [[Rcpp::export]]
arma::field<arma::cube> get_A_all(const arma::cube& gamma, const arma::cube& X, 
                                  const bool tv) {
  unsigned int S = gamma.n_slices;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  unsigned int N = X.n_slices;
  arma::field<arma::cube> A(N);
  arma::mat Atmp(S, S);
  if (tv) {
    for (unsigned int i = 0; i < N; i++) {
      A(i) = arma::cube(S, S, T);
      for (unsigned int t = 0; t < T; t++) { // time
        for (unsigned int j = 0; j < S; j ++) { // from states
          Atmp.col(j) = softmax(gamma.slice(j) * X.slice(i).col(t));
        }
        A(i).slice(t) = Atmp.t();
      }
    }
  } else {
    for (unsigned int i = 0; i < N; i++) {
      A(i) = arma::cube(S, S, T);
      for (unsigned int j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(gamma.slice(j) * X.slice(i).col(0));
      }
      A(i).each_slice() = Atmp.t();
    }
  }
  return A;
}
// gamma is (M - 1) x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
arma::cube get_B(const arma::cube& gamma, const arma::mat& X, 
                 const bool add_missing, const bool tv) {
  unsigned int S = gamma.n_slices;
  unsigned int M = gamma.n_rows;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  arma::cube B(S, M + add_missing, T);
  arma::mat Btmp(M + add_missing, S);
  if (add_missing) {
    Btmp.row(M).fill(1.0);
  }
  if (tv) {
    for (unsigned int t = 0; t < T; t++) { // time
      for (unsigned int j = 0; j < S; j ++) { // from states
        Btmp.col(j).rows(0, M - 1) = softmax(gamma.slice(j) * X.col(t));
      }
      B.slice(t) = Btmp.t();
    }
  } else {
    for (unsigned int j = 0; j < S; j ++) { // from states
      Btmp.col(j).rows(0, M - 1) = softmax(
        gamma.slice(j) * X.col(0)
      );
    }
    B.each_slice() = Btmp.t();
  }
  return B;
}
// gamma is a a field of (M_c - 1) x K x S cubes
// X is K x T (covariates, time point)
arma::field<arma::cube> get_B(
    const arma::field<arma::cube>& gamma, 
    const arma::mat& X, const arma::uvec& M, 
    const bool add_missing, const bool tv) {
  unsigned int C = M.n_elem;
  arma::field<arma::cube> B(C); // C field of cubes, each S x M_c x T
  for (unsigned int c = 0; c < C; c++) {
    B(c) = get_B(gamma(c), X, add_missing, tv);
  }
  return B;
}
// gamma is (M - 1) x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
arma::cube get_log_B(const arma::cube& gamma, const arma::mat& X, 
                     const bool add_missing, const bool tv) {
  unsigned int S = gamma.n_slices;
  unsigned int M = gamma.n_rows;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  arma::cube B(S, M + add_missing, T);
  arma::mat Btmp(M + add_missing, S);
  if (add_missing) {
    Btmp.row(M).fill(1.0);
  }
  if (tv) {
    for (unsigned int t = 0; t < T; t++) { // time
      for (unsigned int j = 0; j < S; j ++) { // from states
        Btmp.col(j).rows(0, M - 1) = softmax(gamma.slice(j) * X.col(t));
      }
      B.slice(t) = Btmp.t();
    }
  } else {
    for (unsigned int j = 0; j < S; j ++) { // from states
      Btmp.col(j).rows(0, M - 1) = softmax(
        gamma.slice(j) * X.col(0)
      );
    }
    B.each_slice() = Btmp.t();
  }
  return arma::log(B);
}
// gamma is a a field of (M_c - 1) x K x S cubes
// X is K x T (covariates, time point)
arma::field<arma::cube> get_log_B(
    const arma::field<arma::cube>& gamma, const arma::mat& X, 
    const arma::uvec& M, const bool add_missing, const bool tv) {
  unsigned int C = M.n_elem;
  arma::field<arma::cube> log_B(C); // C field of cubes, each S x M_c x T
  for (unsigned int c = 0; c < C; c++) {
    log_B(c) = get_log_B(gamma(c), X, add_missing, tv);
  }
  return log_B;
}
// gamma is (M - 1) x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
arma::field<arma::cube> get_B_all(
    const arma::cube& gamma,  const arma::cube& X, 
    const bool add_missing, const bool tv) {
  unsigned int S = gamma.n_slices;
  unsigned int M = gamma.n_rows;
  unsigned int K = X.n_rows;
  unsigned int T = X.n_cols;
  unsigned int N = X.n_slices;
  arma::field<arma::cube> B(N);
  arma::mat Btmp(M + add_missing, S);
  if (add_missing) {
    Btmp.row(M).fill(1.0);
  }
  if (tv) {
    for (unsigned int i = 0; i < N; i++) {
      B(i) = arma::cube(S, M + add_missing, T);
      for (unsigned int t = 0; t < T; t++) { // time
        for (unsigned int j = 0; j < S; j ++) { // from states
          Btmp.col(j).rows(0, M - 1) = softmax(
            gamma.slice(j) * X.slice(i).col(t)
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
          gamma.slice(j) * X.slice(i).col(0)
        );
      }
      B(i).each_slice() = Btmp.t();
    }
  }
  return B;
}
