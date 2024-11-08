#include "get_parameters.h"

// gamma_omega is D x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_omega(const arma::mat& gamma, const arma::vec& X) {
  return softmax(gamma * X);
}
// eta_omega is D x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_log_omega(const arma::mat& gamma, const arma::vec& X) {
  return arma::log(softmax(gamma * X));
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
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_pi(const arma::mat& gamma, const arma::vec& X) {
  return softmax(gamma * X);
}
// gamma is S x K (start from, covariates)
// X a vector of length K
// [[Rcpp::export]]
arma::vec get_log_pi(const arma::mat& gamma, const arma::vec& X) {
  return arma::log(softmax(gamma * X));
}
// gamma is S x K x S (transition to, covariates, transition from)
// X is K x T matrix (covariates, time points)
// [[Rcpp::export]]
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
// gamma is S x K x S (transition to, covariates, transition from)
// X is K x T matrix (covariates, time points)
// [[Rcpp::export]]
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
// gamma is M x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
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
// gamma is a a field of M_c x K x S cubes
// X is K x T (covariates, time point)
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
// gamma is M x K x S (symbols, covariates, transition from)
// X is K x T (covariates, time points)
// [[Rcpp::export]]
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
// gamma is a a field of M_c x K x S cubes
// X is K x T (covariates, time point)
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
// output is P x N field of (S x S x T) cubes
// [[Rcpp::export]]
arma::mat get_A_qs(const arma::field<arma::cube>& gamma, 
                   const arma::cube& X, const bool tv,
                   const arma::vec& probs) {
  
  arma::uword S = gamma(0).n_rows;
  arma::uword L = gamma.n_elem;
  arma::uword N = X.n_slices;
  arma::uword T = X.n_cols;
  arma::uword P = probs.n_elem;
  arma::uword SST = S * S * T;
  arma::mat A(SST, L);
  arma::mat qs(SST * N, P);
  for (arma::uword i = 0; i < N; i++) {
    for (arma::uword l = 0; l < L; l++) {
      A.col(l) = arma::vectorise(get_A(gamma(l), X.slice(i), tv));
    }
    qs.rows(i * SST, (i + 1) * SST - 1) = arma::quantile(A, probs, 1);
  }
  return qs;
}

// gamma is B field of M x K x S cubes
// X is K x T x N cube (covariates, time points, sequences)
//  output is T x N field of (S x M x T) cubes
// [[Rcpp::export]]
arma::mat get_B_qs(const arma::field<arma::cube>& gamma, 
                   const arma::cube& X, const bool tv,
                   const arma::vec& probs) {
  arma::uword M = gamma(0).n_rows;
  arma::uword S = gamma(0).n_slices;
  arma::uword L = gamma.n_elem;
  arma::uword N = X.n_slices;
  arma::uword T = X.n_cols;
  arma::uword P = probs.n_elem;
  arma::uword SMT = S * M * T;
  arma::mat B(SMT, L);
  arma::mat qs(SMT * N, P);
  for (arma::uword i = 0; i < N; i++) {
    for (arma::uword l = 0; l < L; l++) {
      B.col(l) = arma::vectorise(get_B(gamma(l), X.slice(i), tv));
    }
    qs.rows(i * SMT, (i + 1) * SMT - 1) = arma::quantile(B, probs, 1);
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


// gamma is S x K x L (start from, covariates, sample)
// X a K x N matrix
// probs is a vector of length P
// [[Rcpp::export]]
arma::mat get_pi_ame(const arma::field<arma::mat>& gamma, 
                     const arma::mat& X1, const arma::mat& X2, 
                     const arma::vec& probs) {
  arma::uword S = gamma(0).n_rows;
  arma::uword L = gamma.n_elem;
  arma::uword N = X1.n_cols;
  double invN = 1.0 / N;
  arma::mat pi(S, L, arma::fill::zeros);
  for (arma::uword l = 0; l < L; l++) {
    for (arma::uword i = 0; i < N; i++) {
      pi.col(l) += invN * (
        get_pi(gamma(l), X1.col(i)) - get_pi(gamma(l), X2.col(i))
      );
    }
  }
  return arma::quantile(pi, probs, 1);;
}
// gamma is L field of S x K x S cubes
// X is K x T x N cube (covariates, time points, sequences)
// probs is a vector of length P
// [[Rcpp::export]]
arma::mat get_A_ame(const arma::field<arma::cube>& gamma, 
                    const arma::cube& X1, const arma::cube& X2, const bool tv,
                    const arma::vec& probs) {
  
  arma::uword S = gamma(0).n_rows;
  arma::uword L = gamma.n_elem;
  arma::uword N = X1.n_slices;
  arma::uword T = X1.n_cols;
  double invN = 1.0 / N;
  arma::uword SST = S * S * T;
  arma::mat A(SST, L, arma::fill::zeros);
  for (arma::uword l = 0; l < L; l++) {
    for (arma::uword i = 0; i < N; i++) {
      A.col(l) += invN * (
        arma::vectorise(get_A(gamma(l), X1.slice(i), tv)) - 
          arma::vectorise(get_A(gamma(l), X2.slice(i), tv))
      );
    }
  }
  return arma::quantile(A, probs, 1);
}

// gamma is L field of M x K x S cubes
// X is K x T x N cube (covariates, time points, sequences)
// [[Rcpp::export]]
arma::mat get_B_ame(const arma::field<arma::cube>& gamma, 
                    const arma::cube& X1, const arma::cube& X2, const bool tv,
                    const arma::vec& probs) {
  arma::uword M = gamma(0).n_rows;
  arma::uword S = gamma(0).n_slices;
  arma::uword L = gamma.n_elem;
  arma::uword N = X1.n_slices;
  arma::uword T = X1.n_cols;
  arma::uword SMT = S * M * T;
  double invN = 1.0 / N;
  arma::mat B(SMT, L, arma::fill::zeros);
  for (arma::uword l = 0; l < L; l++) {
    for (arma::uword i = 0; i < N; i++) {
      B.col(l) += invN * (
        arma::vectorise(get_B(gamma(l), X1.slice(i), tv)) - 
          arma::vectorise(get_B(gamma(l), X2.slice(i), tv))
      );
    }
  }
  return arma::quantile(B, probs, 1);
}
// gamma is D x K x L (cluster, covariates, sample)
// X a K x N matrix
// [[Rcpp::export]]
arma::mat get_omega_ame(const arma::field<arma::mat>& gamma, 
                     const arma::mat& X1, const arma::mat& X2, 
                     const arma::vec& probs) {
  arma::uword D = gamma(0).n_rows;
  arma::uword L = gamma.n_elem;
  arma::uword N = X1.n_cols;
  double invN = 1.0 / N;
  arma::mat omega(D, L, arma::fill::zeros);
  for (arma::uword l = 0; l < L; l++) {
    for (arma::uword i = 0; i < N; i++) {
      omega.col(l) += invN * (
        get_omega(gamma(l), X1.col(i)) - get_omega(gamma(l), X2.col(i))
      );
    }
  }
  return arma::quantile(omega, probs, 1);
}
// 
// // Compute row of A for EM algorithm
// // gamma is S x K (transition to, covariates)
// // X is K x T matrix (covariates, time points)
// // [[Rcpp::export]]
// arma::mat get_A_em(const arma::mat& gamma, const arma::mat& X, 
//                  const bool tv) {
//   arma::uword S = gamma.rows;
//   arma::uword T = X.n_cols;
//   arma::mat A(S, T);
//   arma::mat Atmp(S, S);
//   if (tv) {
//     for (arma::uword t = 0; t < T; t++) { // time
//       for (arma::uword j = 0; j < S; j ++) { // from states
//         Atmp.col(j) = softmax(gamma.slice(j) * X.col(t));
//       }
//       A.slice(t) = Atmp.t();
//     }
//   } else {
//     for (arma::uword j = 0; j < S; j ++) { // from states
//       Atmp.col(j) = softmax(gamma.slice(j) * X.col(0));
//     }
//     A.each_slice() = Atmp.t();
//   }
//   return A;
// }
