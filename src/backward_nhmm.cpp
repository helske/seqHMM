// backward algorithm for NHMM
#include "backward_nhmm.h"
#include "get_parameters.h"
#include "eta_to_gamma.h"

// [[Rcpp::export]]
arma::cube backward_nhmm_singlechannel(
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::cube& eta_B, const arma::cube& X_o,
    const arma::umat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int M = eta_B.n_rows + 1;
  arma::cube log_beta(S, T, N);
  arma::mat log_py(S, T);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::cube gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_A = get_log_A(gamma_A, X_s.slice(i));
    log_B = get_log_B(gamma_B, X_o.slice(i), true);
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = log_B(s, obs(t, i), t);
      }
    }
    univariate_backward_nhmm(log_beta.slice(i), log_A, log_py);
  }
  return log_beta;
}

// [[Rcpp::export]]
arma::cube backward_nhmm_multichannel(
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::ucube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int C = obs.n_rows;
  arma::cube log_beta(S, T, N);
  arma::mat log_py(S, T);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::field gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_py.zeros();
    log_A = get_log_A(gamma_A, X_s.slice(i));
    log_B = get_log_B(gamma_B, X_o.slice(i), M, true);
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int c = 0; c < C; c++) {
        log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
      }
    }
    univariate_backward_nhmm(log_beta.slice(i), log_A, log_py);
  }
  return log_beta;
}

// [[Rcpp::export]]
arma::cube backward_mnhmm_singlechannel(
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::umat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_A.n_elem;
  unsigned int M = eta_B(0).n_rows + 1;
  arma::cube log_beta(S * D, T, N);
  arma::mat log_py(S, T);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  arma::field gamma_A = eta_to_gamma(eta_A);
  arma::field gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int d = 0; d < D; d++) {
      log_A = get_log_A(gamma_A(d), X_s.slice(i));
      log_B = get_log_B(gamma_B(d), X_o.slice(i), true);
      for (unsigned int t = 0; t < T; t++) {
        log_py.col(t) = log_B.slice(t).col(obs(t, i));
      }
      arma::subview<double> submat = log_beta.slice(i).rows(d * S, (d + 1) * S - 1);
      univariate_backward_nhmm(submat, log_A, log_py);
    }
  }
  return log_beta;
}
// [[Rcpp::export]]
arma::cube backward_mnhmm_multichannel(
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::ucube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_A.n_elem;
  unsigned int C = obs.n_rows;
  arma::cube log_beta(S * D, T, N);
  arma::mat log_py(S, T);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  arma::field gamma_A = eta_to_gamma(eta_A);
  arma::field gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int d = 0; d < D; d++) {
      log_py.zeros();
      log_A = get_log_A(gamma_A(d), X_s.slice(i));
      log_B = get_log_B(
        gamma_B.rows(d * C, (d + 1) * C - 1), X_o.slice(i), M, true
      );
      for (unsigned int t = 0; t < T; t++) {
        for (unsigned int c = 0; c < C; c++) {
          log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
        }
      }
      arma::subview<double> submat = log_beta.slice(i).rows(d * S, (d + 1) * S - 1);
      univariate_backward_nhmm(submat, log_A, log_py);
    }
  }
  return log_beta;
}
