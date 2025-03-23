#include "get_parameters.h"

arma::vec get_omega(const arma::mat& gamma, const arma::vec& X) {
  return softmax(gamma * X);
}
arma::vec get_log_omega(const arma::mat& gamma, const arma::vec& X) {
  return arma::log(softmax(gamma * X));
}


arma::vec get_pi(const arma::mat& gamma, const arma::vec& X) {
  return softmax(gamma * X);
}
arma::vec get_log_pi(const arma::mat& gamma, const arma::vec& X) {
  return arma::log(softmax(gamma * X));
}
arma::cube get_A(const arma::cube& gamma, const arma::mat& X, 
                 const bool tv) {
  arma::uword S = gamma.n_slices;
  arma::uword T = X.n_cols;
  arma::cube A(S, S, T);
  arma::mat Atmp(S, S);
  if (tv) {
    for (arma::uword t = 0; t < T; t++) { // time
      for (arma::uword j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(gamma.slice(j) * X.col(t));
      }
      A.slice(t) = Atmp.t();
    }
  } else {
    for (arma::uword j = 0; j < S; j ++) { // from states
      Atmp.col(j) = softmax(gamma.slice(j) * X.col(0));
    }
    A.each_slice() = Atmp.t();
  }
  return A;
}
arma::cube get_log_A(const arma::cube& gamma, const arma::mat& X, 
                     const bool tv) {
  arma::uword S = gamma.n_slices;
  arma::uword T = X.n_cols;
  arma::cube A(S, S, T);
  arma::mat Atmp(S, S);
  if (tv) {
    for (arma::uword t = 0; t < T; t++) { // time
      for (arma::uword j = 0; j < S; j ++) { // from states
        Atmp.col(j) = softmax(gamma.slice(j) * X.col(t));
      }
      A.slice(t) = Atmp.t();
    }
  } else {
    for (arma::uword j = 0; j < S; j ++) { // from states
      Atmp.col(j) = softmax(gamma.slice(j) * X.col(0));
    }
    A.each_slice() = Atmp.t();
  }
  return arma::log(A);
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
    for (arma::uword t = 0; t < T; t++) { // time
      for (arma::uword j = 0; j < S; j ++) { // from states
        Btmp.col(j).rows(0, M - 1) = softmax(gamma.slice(j) * X.col(t));
      }
      B.slice(t) = Btmp.t();
    }
  } else {
    for (arma::uword j = 0; j < S; j ++) { // from states
      Btmp.col(j).rows(0, M - 1) = softmax(
        gamma.slice(j) * X.col(0)
      );
    }
    B.each_slice() = Btmp.t();
  }
  return B;
}
arma::field<arma::cube> get_B(
    const arma::field<arma::cube>& gamma, const arma::mat& X, 
    const arma::uvec& M, const bool tv, const bool add_missing) {
  arma::uword C = M.n_elem;
  arma::field<arma::cube> B(C); // C field of cubes, each S x M_c x T
  for (arma::uword c = 0; c < C; c++) {
    B(c) = get_B(gamma(c), X, tv, add_missing);
  }
  return B;
}
arma::cube get_log_B(const arma::cube& gamma, const arma::mat& X, 
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
    for (arma::uword t = 0; t < T; t++) { // time
      for (arma::uword j = 0; j < S; j ++) { // from states
        Btmp.col(j).rows(0, M - 1) = softmax(gamma.slice(j) * X.col(t));
      }
      B.slice(t) = Btmp.t();
    }
  } else {
    for (arma::uword j = 0; j < S; j ++) { // from states
      Btmp.col(j).rows(0, M - 1) = softmax(
        gamma.slice(j) * X.col(0)
      );
    }
    B.each_slice() = Btmp.t();
  }
  return arma::log(B);
}
arma::field<arma::cube> get_log_B(
    const arma::field<arma::cube>& gamma, const arma::mat& X, 
    const arma::uvec& M, const bool tv, const bool add_missing) {
  arma::uword C = M.n_elem;
  arma::field<arma::cube> log_B(C); // C field of cubes, each S x M_c x T
  for (arma::uword c = 0; c < C; c++) {
    log_B(c) = get_log_B(gamma(c), X, tv, add_missing);
  }
  return log_B;
}

// [[Rcpp::export]]
arma::mat get_omega_all(const arma::mat& gamma, const arma::mat& X) {
  arma::mat omega(gamma.n_rows, X.n_cols);
  for (arma::uword i = 0; i < X.n_cols; i++) {
    omega.col(i) = softmax(gamma * X.col(i));
  }
  return omega;
}
// gamma is S x K (start from, covariates)
// X a K x N matrix
// [[Rcpp::export]]
arma::mat get_pi_all(const arma::mat& gamma, const arma::mat& X) {
  arma::mat pi(gamma.n_rows, X.n_cols);
  for (arma::uword i = 0; i < X.n_cols; i++) {
    pi.col(i) = softmax(gamma * X.col(i));
  }
  return pi;
}
// gamma is S x K x S (transition to, covariates, transition from)
// X is K x T x N cube (covariates, time points, sequences)
// [[Rcpp::export]]
arma::field<arma::cube> get_A_all(const arma::cube& gamma, const arma::cube& X, 
                                  const bool tv) {
  arma::uword N = X.n_slices;
  arma::field<arma::cube> A(N);
  for (arma::uword i = 0; i < N; i++) {
    A(i) = get_A(gamma, X.slice(i), tv);
  }
  return A;
}
// gamma is M x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
arma::field<arma::cube> get_B_all(
    const arma::cube& gamma,  const arma::cube& X, const bool tv) {
  arma::uword N = X.n_slices;
  arma::field<arma::cube> B(N);
  for (arma::uword i = 0; i < N; i++) {
    B(i) = get_B(gamma, X.slice(i), tv);
  }
  return B;
}

// gamma is S x K x L (start from, covariates, sample)
// X a K x N matrix
// probs is a vector of length P
// output is S x N x P cube
// [[Rcpp::export]]
arma::mat get_pi_qs(const arma::field<arma::mat>& gamma, const arma::mat& X, 
                    const arma::vec& probs) {
  arma::uword S = gamma(0).n_rows;
  arma::uword L = gamma.n_elem;
  arma::uword N = X.n_cols;
  arma::uword P = probs.n_elem;
  arma::mat pi(S, L);
  arma::mat qs(S * N, P);
  for (arma::uword i = 0; i < N; i++) {
    for (arma::uword l = 0; l < L; l++) {
      pi.col(l) = get_pi(gamma(l), X.col(i));
    }
    qs.rows(i * S, (i + 1) * S - 1) = arma::quantile(pi, probs, 1);
  }
  return qs;
}
// gamma is L field of S x K x S cubes
// X is K x T x N cube (covariates, time points, sequences)
// probs is a vector of length P
// output is SST * N x P matrix
// [[Rcpp::export]]
arma::mat get_A_qs(const arma::field<arma::cube>& gamma, 
                   const arma::cube& X, const bool tv,
                   const arma::vec& probs, const arma::uvec& Ti) {
  
  arma::uword S = gamma(0).n_rows;
  arma::uword L = gamma.n_elem;
  arma::uword N = X.n_slices;
  arma::uword T = X.n_cols;
  arma::uword P = probs.n_elem;
  arma::uword SST = S * S * T;
  arma::mat A(SST, L, arma::fill::zeros);
  arma::mat qs(SST * N, P);
  for (arma::uword i = 0; i < N; i++) {
    for (arma::uword l = 0; l < L; l++) {
      A.col(l) = arma::vectorise(get_A(gamma(l), X.slice(i), tv));
    }
    //qs.rows(i * SST, (i + 1) * SST - 1) = arma::quantile(A, probs, 1);
    qs.rows(i * SST, i * SST + S * S * Ti(i) - 1) = arma::quantile(A.rows(0, S * S * Ti(i) - 1), probs, 1);
  }
  return qs;
}

// gamma is B field of M x K x S cubes
// X is K x T x N cube (covariates, time points, sequences)
//  output is T x N field of (S x M x T) cubes
// [[Rcpp::export]]
arma::mat get_B_qs(const arma::field<arma::cube>& gamma, 
                   const arma::cube& X, const bool tv,
                   const arma::vec& probs, const arma::uvec& Ti) {
  arma::uword M = gamma(0).n_rows;
  arma::uword S = gamma(0).n_slices;
  arma::uword L = gamma.n_elem;
  arma::uword N = X.n_slices;
  arma::uword T = X.n_cols;
  arma::uword P = probs.n_elem;
  arma::uword SMT = S * M * T;
  arma::mat B(SMT, L);
  arma::mat qs(SMT * N, P, arma::fill::zeros);
  for (arma::uword i = 0; i < N; i++) {
    for (arma::uword l = 0; l < L; l++) {
      B.col(l) = arma::vectorise(get_B(gamma(l), X.slice(i), tv));
    }
    //qs.rows(i * SMT, (i + 1) * SMT - 1) = arma::quantile(B, probs, 1);
    qs.rows(i * SMT, i * SMT + S * M * Ti(i) - 1) = arma::quantile(B.rows(0, S * M * Ti(i) - 1), probs, 1);
  }
  return qs;
}
// gamma is D x K x B (cluster, covariates, sample)
// X a K x N matrix
// output is D x N x B cube
// [[Rcpp::export]]
arma::mat get_omega_qs(const arma::field<arma::mat>& gamma, const arma::mat& X,
                       const arma::vec& probs) {
  arma::uword D = gamma(0).n_rows;
  arma::uword L = gamma.n_elem;
  arma::uword N = X.n_cols;
  arma::uword P = probs.n_elem;
  arma::mat omega(D, L);
  arma::mat qs(D * N, P);
  for (arma::uword i = 0; i < N; i++) {
    for (arma::uword l = 0; l < L; l++) {
      omega.col(l) = get_omega(gamma(l), X.col(i));
    }
    qs.rows(i * D, (i + 1) * D - 1) = arma::quantile(omega, probs, 1);
  }
  return qs;
}
