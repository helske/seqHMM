#include "get_parameters.h"
#include "list_to_field.h"

arma::cube get_A(const arma::cube& gamma, const arma::mat& X, 
                 const bool tv) {
  arma::uword S = gamma.n_slices;
  arma::uword T = X.n_cols;
  arma::cube A(S, S, T);
  arma::mat Atmp(S, S);
  if (tv) {
    for (arma::uword t = 0; t < T; ++t) { // time
      for (arma::uword j = 0; j < S; ++j) { // from states
        Atmp.col(j) = softmax(gamma.slice(j) * X.col(t));
      }
      A.slice(t) = Atmp.t();
    }
  } else {
    for (arma::uword j = 0; j < S; ++j) { // from states
      Atmp.col(j) = softmax(gamma.slice(j) * X.col(0));
    }
    A.each_slice() = Atmp.t();
  }
  return A;
}

arma::cube get_B(const arma::cube& gamma, const arma::mat& X, 
                 const bool tv, const bool add_missing) {
  arma::uword S = gamma.n_slices;
  arma::uword M = gamma.n_rows;
  arma::uword T = X.n_cols;
  arma::cube B(S, M + add_missing, T);
  arma::mat Btmp(M + add_missing, S);
  if (add_missing) {
    Btmp.row(M).fill(1.0);
  }
  if (tv) {
    for (arma::uword t = 0; t < T; ++t) { // time
      for (arma::uword j = 0; j < S; ++j) { // from states
        Btmp.col(j).rows(0, M - 1) = softmax(gamma.slice(j) * X.col(t));
      }
      B.slice(t) = Btmp.t();
    }
  } else {
    for (arma::uword j = 0; j < S; ++j) { // from states
      Btmp.col(j).rows(0, M - 1) = softmax(
        gamma.slice(j) * X.col(0)
      );
    }
    B.each_slice() = Btmp.t();
  }
  return B;
}

// [[Rcpp::export]]
arma::mat get_omega_all(const arma::mat& gamma, const arma::mat& X) {
  arma::mat omega(gamma.n_rows, X.n_cols);
  for (arma::uword i = 0; i < X.n_cols; ++i) {
    omega.col(i) = softmax(gamma * X.col(i));
  }
  return omega;
}
// gamma is S x K (start from, covariates)
// X a K x N matrix
// [[Rcpp::export]]
arma::mat get_pi_all(const arma::mat& gamma, const arma::mat& X) {
  arma::mat pi(gamma.n_rows, X.n_cols);
  for (arma::uword i = 0; i < X.n_cols; ++i) {
    pi.col(i) = softmax(gamma * X.col(i));
  }
  return pi;
}
// gamma is S x K x S (transition to, covariates, transition from)
// X is K x T x N cube (covariates, time points, sequences)
// [[Rcpp::export]]
arma::field<arma::cube> get_A_all(
    const arma::cube& gamma, const arma::field<arma::mat>& X, const bool tv) {
  arma::uword N = X.n_elem;
  arma::field<arma::cube> A(N);
  for (arma::uword i = 0; i < N; ++i) {
    A(i) = get_A(gamma, X(i), tv);
  }
  return A;
}
// gamma is M x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
arma::field<arma::cube> get_B_all(
    const arma::cube& gamma,  const arma::field<arma::mat>& X, const bool tv) {
  arma::uword N = X.n_elem;
  arma::field<arma::cube> B(N);
  for (arma::uword i = 0; i < N; ++i) {
    B(i) = get_B(gamma, X(i), tv);
  }
  return B;
}

// gamma is M x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
arma::field<arma::cube> get_B1(
    const arma::field<arma::cube>& gamma, const arma::uvec M, 
    const Rcpp::List& W_X_B, const arma::vec& prior_y) {
  arma::uword S = gamma(0).n_slices;
  arma::uword C = M.n_elem;
  arma::field<arma::vec> X = veclist_to_3d_field(W_X_B);
  arma::uword N = X.n_slices;
  arma::field<arma::cube> B(C);
  for (arma::uword c = 0; c < C; ++c) {
    arma::mat Btmp(M(c), S);
    B(c) = arma::cube(S, M(c), N);
    for (arma::uword i = 0; i < N; ++i) {
      Btmp.zeros();
      for (arma::uword j = 0; j < prior_y.n_elem; ++j) {
        for (arma::uword s = 0; s < S; ++s) { // from states
          Btmp.col(s) += softmax(
            gamma(c).slice(s) * X(j, c, i)
          ) * prior_y(j);
        }
      }
      B(c).slice(i) = Btmp.t();
    }
  }
  return B;
}
