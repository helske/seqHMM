// forward algorithm for NHMM
#include "forward_nhmm.h"
#include "get_parameters.h"
#include "eta_to_gamma.h"

// [[Rcpp::export]]
arma::cube forward_nhmm_singlechannel(
    const arma::mat& eta_pi, const arma::mat& X_i,
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::cube& eta_B, const arma::cube& X_o,
    const arma::umat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int M = eta_B.n_rows + 1;
  arma::cube log_alpha(S, T, N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  arma::mat gamma_pi = eta_to_gamma(eta_pi);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::cube gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_Pi = get_log_pi(gamma_pi, X_i.col(i));
    log_A = get_log_A(gamma_A, X_s.slice(i));
    log_B = get_log_B(gamma_B, X_o.slice(i), true);
    for (unsigned int t = 0; t < T; t++) {
      log_py.col(t) = log_B.slice(t).col(obs(t, i));
    }
    univariate_forward_nhmm(log_alpha.slice(i), log_Pi, log_A, log_py);
  }
  return log_alpha;
}

// [[Rcpp::export]]
arma::cube forward_nhmm_multichannel(
    const arma::mat& eta_pi, const arma::mat& X_i,
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::ucube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int C = obs.n_rows;
  arma::cube log_alpha(S, T, N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  arma::mat gamma_pi = eta_to_gamma(eta_pi);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_py.zeros();
    log_Pi = get_log_pi(gamma_pi, X_i.col(i));
    log_A = get_log_A(gamma_A, X_s.slice(i));
    log_B = get_log_B(gamma_B, X_o.slice(i), M, true);
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int c = 0; c < C; c++) {
        log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
      }
    }
    univariate_forward_nhmm(log_alpha.slice(i), log_Pi, log_A, log_py);
  }
  return log_alpha;
}

// [[Rcpp::export]]
arma::cube forward_mnhmm_singlechannel(
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_i,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d,
    const arma::umat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int M = eta_B(0).n_rows + 1;
  arma::cube log_alpha(S * D, T, N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  arma::vec log_omega(D);
  arma::mat gamma_omega = eta_to_gamma(eta_omega);
  arma::field<arma::mat> gamma_pi = eta_to_gamma(eta_pi);
  arma::field<arma::cube> gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_omega = get_log_omega(gamma_omega, X_d.col(i));
    for (unsigned int d = 0; d < D; d++) {
      log_Pi = get_log_pi(gamma_pi(d), X_i.col(i));
      log_A = get_log_A(gamma_A(d), X_s.slice(i));
      log_B = get_log_B(gamma_B(d), X_o.slice(i), true);
      for (unsigned int t = 0; t < T; t++) {
        log_py.col(t) = log_B.slice(t).col(obs(t, i));
      }
      log_alpha.slice(i).rows(d * S, (d + 1) * S - 1).fill(log_omega(d)); 
      arma::subview<double> submat = log_alpha.slice(i).rows(d * S, (d + 1) * S - 1);
      univariate_forward_nhmm(submat, log_Pi, log_A, log_py);
    }
  }
  return log_alpha;
}
// [[Rcpp::export]]
arma::cube forward_mnhmm_multichannel(
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_i,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d,
    const arma::ucube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int C = obs.n_rows;
  arma::cube log_alpha(S * D, T, N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  arma::vec log_omega(D);
  arma::mat gamma_omega = eta_to_gamma(eta_omega);
  arma::field<arma::mat> gamma_pi = eta_to_gamma(eta_pi);
  arma::field<arma::cube> gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_omega = get_log_omega(gamma_omega, X_d.col(i));
    for (unsigned int d = 0; d < D; d++) {
      log_py.zeros();
      log_Pi = get_log_pi(gamma_pi(d), X_i.col(i));
      log_A = get_log_A(gamma_A(d), X_s.slice(i));
      log_B = get_log_B(
        gamma_B.rows(d * C, (d + 1) * C - 1), X_o.slice(i), M, true
      );
      for (unsigned int t = 0; t < T; t++) {
        for (unsigned int c = 0; c < C; c++) {
          log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
        }
      }
      log_alpha.slice(i).rows(d * S, (d + 1) * S - 1).fill(log_omega(d)); 
      arma::subview<double> submat = log_alpha.slice(i).rows(d * S, (d + 1) * S - 1);
      univariate_forward_nhmm(submat, log_Pi, log_A, log_py);
    }
  }
  return log_alpha;
}
