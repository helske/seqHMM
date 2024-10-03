// backward algorithm for NHMM
#include "backward_nhmm.h"
#include "get_parameters.h"
#include "logsumexp.h"

arma::mat univariate_backward_nhmm(
    const arma::cube& log_transition, 
    const arma::mat& log_py) {
  
  unsigned int S = log_py.n_rows;
  unsigned int T = log_py.n_cols;
  
  arma::mat log_beta(S, T);
  log_beta.col(T - 1).zeros();
  for (int t = (T - 2); t >= 0; t--) {
    arma::vec tmpbeta(S);
    for (unsigned int i = 0; i < S; i++) {
      log_beta(i, t) = logSumExp(
        log_beta.col(t + 1) + log_transition.slice(t).row(i).t() + log_py.col(t + 1)
      );
    }
  }
  return log_beta;
}

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
  for (unsigned int i = 0; i < N; i++) {
    log_A = get_log_A(eta_A, X_s.slice(i));
    log_B = get_log_B(eta_B, X_o.slice(i), true);
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = log_B(s, obs(t, i), t);
      }
    }
    log_beta.slice(i) = univariate_backward_nhmm(log_A, log_py);
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
  for (unsigned int i = 0; i < N; i++) {
    log_py.zeros();
    log_A = get_log_A(eta_A, X_s.slice(i));
    log_B = get_log_B(eta_B, X_o.slice(i), M, true);
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int c = 0; c < C; c++) {
        log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
      }
    }
    log_beta.slice(i) = univariate_backward_nhmm(log_A, log_py);
  }
  return log_beta;
}

// [[Rcpp::export]]
arma::cube backward_mnhmm_singlechannel(
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d,
    const arma::umat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int M = eta_B(0).n_rows + 1;
  arma::cube log_beta(S * D, T, N);
  arma::mat log_py(S, T);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int d = 0; d < D; d++) {
      log_A = get_log_A(eta_A(d), X_s.slice(i));
      log_B = get_log_B(eta_B(d), X_o.slice(i), true);
      for (unsigned int t = 0; t < T; t++) {
        log_py.col(t) = log_B.slice(t).col(obs(t, i));
      }
      log_beta.slice(i).rows(d * S, (d + 1) * S - 1) = univariate_backward_nhmm(log_A, log_py);
    }
  }
  return log_beta;
}
// [[Rcpp::export]]
arma::cube backward_mnhmm_multichannel(
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d,
    const arma::ucube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int C = obs.n_rows;
  arma::cube log_beta(S * D, T, N);
  arma::mat log_py(S, T);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int d = 0; d < D; d++) {
      log_py.zeros();
      log_A = get_log_A(eta_A(d), X_s.slice(i));
      log_B = get_log_B(
        eta_B.rows(d * C, (d + 1) * C - 1), X_o.slice(i), M, true
      );
      for (unsigned int t = 0; t < T; t++) {
        for (unsigned int c = 0; c < C; c++) {
          log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
        }
      }
      log_beta.slice(i).rows(d * S, (d + 1) * S - 1) = univariate_backward_nhmm(log_A, log_py);
    }
  }
  return log_beta;
}
