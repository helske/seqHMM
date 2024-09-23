// forward algorithm for NHMM
#include "forward_nhmm.h"
#include "get_parameters.h"
#include "logsumexp.h"

arma::mat univariate_forward_nhmm(
    const arma::vec& log_init, 
    const arma::cube& log_transition, 
    const arma::mat& log_py) {
  
  unsigned int S = log_py.n_rows;
  unsigned int T = log_py.n_cols;
  
  arma::mat log_alpha(S, T);
  log_alpha.col(0) = log_init + log_py.col(0);
  for (unsigned int t = 1; t < T; t++) {
    for (unsigned int i = 0; i < S; i++) {
      log_alpha(i, t) = logSumExp(
        log_alpha.col(t - 1) + log_transition.slice(t - 1).col(i) + log_py(i, t)
      );
    }
  }
  return log_alpha;
}

// [[Rcpp::export]]
arma::cube forward_nhmm_singlechannel(
    const arma::mat& gamma_pi_raw, const arma::mat& X_i,
    const arma::cube& gamma_A_raw, const arma::cube& X_s,
    const arma::cube& gamma_B_raw, const arma::cube& X_o,
    const arma::mat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = gamma_A_raw.n_slices;
  unsigned int M = gamma_B_raw.n_rows + 1;
  arma::cube log_alpha(S, T, N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  for (unsigned int i = 0; i < N; i++) {
    log_Pi = get_pi(gamma_pi_raw, X_i.col(i), true);
    log_A = get_A(gamma_A_raw, X_s.slice(i), true);
    log_B = get_B(gamma_B_raw, X_o.slice(i), true, true);
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = log_B(s, obs(t, i), t);
      }
    }
    log_alpha.slice(i) = univariate_forward_nhmm(log_Pi, log_A, log_py);
  }
  return log_alpha;
}

// [[Rcpp::export]]
arma::cube forward_nhmm_multichannel(
    const arma::mat& gamma_pi_raw, const arma::mat& X_i,
    const arma::cube& gamma_A_raw, const arma::cube& X_s,
    const arma::field<arma::cube>& gamma_B_raw, const arma::cube& X_o,
    const arma::cube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = gamma_A_raw.n_slices;
  unsigned int C = obs.n_rows;
  arma::cube log_alpha(S, T, N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  for (unsigned int i = 0; i < N; i++) {
    log_Pi = get_pi(gamma_pi_raw, X_i.col(i), true);
    log_A = get_A(gamma_A_raw, X_s.slice(i), true);
    log_B = get_B(gamma_B_raw, X_o.slice(i), M, true, true);
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int s = 0; s < S; s++) {
        log_py(s, t) = 0;
        for (unsigned int c = 0; c < C; c++) {
          log_py(s, t) += log_B(c)(s, obs(c, t, i), t);
        }
      }
    }
    log_alpha.slice(i) = univariate_forward_nhmm(log_Pi, log_A, log_py);
  }
  return log_alpha;
}


// [[Rcpp::export]]
arma::cube forward_mnhmm_singlechannel(
    const arma::field<arma::mat>& gamma_pi_raw, const arma::mat& X_i,
    const arma::field<arma::cube>& gamma_A_raw, const arma::cube& X_s,
    const arma::field<arma::cube>& gamma_B_raw, const arma::cube& X_o,
    const arma::mat& gamma_omega_raw, const arma::mat& X_d,
    const arma::mat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = gamma_A_raw(0).n_slices;
  unsigned int D = gamma_omega_raw.n_rows + 1;
  unsigned int M = gamma_B_raw(0).n_rows + 1;
  arma::cube log_alpha(S * D, T, N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::cube log_B(S, M + 1, T);
  arma::vec log_omega(D);
  for (unsigned int i = 0; i < N; i++) {
    log_omega = get_omega(gamma_omega_raw, X_d.col(i), true);
    for (unsigned int d = 0; d < D; d++) {
      log_Pi = get_pi(gamma_pi_raw(d), X_i.col(i), true);
      log_A = get_A(gamma_A_raw(d), X_s.slice(i), true);
      log_B = get_B(gamma_B_raw(d), X_o.slice(i), true, true);
      for (unsigned int t = 0; t < T; t++) {
        for (unsigned int s = 0; s < S; s++) {
          log_py(s, t) = log_B(s, obs(t, i), t);
        }
      }
      log_alpha.slice(i).rows(d * S, (d + 1) * S - 1) = log_omega(d) + 
        univariate_forward_nhmm(log_Pi, log_A, log_py);
    }
  }
  return log_alpha;
}
// [[Rcpp::export]]
arma::cube forward_mnhmm_multichannel(
    const arma::field<arma::mat>& gamma_pi_raw, const arma::mat& X_i,
    const arma::field<arma::cube>& gamma_A_raw, const arma::cube& X_s,
    const arma::field<arma::cube>& gamma_B_raw, const arma::cube& X_o,
    const arma::mat& gamma_omega_raw, const arma::mat& X_d,
    const arma::cube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = gamma_A_raw(0).n_slices;
  unsigned int D = gamma_omega_raw.n_rows + 1;
  unsigned int C = obs.n_rows;
  arma::cube log_alpha(S * D, T, N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  arma::vec log_omega(D);
  for (unsigned int i = 0; i < N; i++) {
    log_omega = get_omega(gamma_omega_raw, X_d.col(i), true);
    for (unsigned int d = 0; d < D; d++) {
      log_Pi = get_pi(gamma_pi_raw(d), X_i.col(i), true);
      log_A = get_A(gamma_A_raw(d), X_s.slice(i), true);
      log_B = get_B(
        gamma_B_raw.rows(d * C, (d + 1) * C - 1), X_o.slice(i), M, true, true
      );
      for (unsigned int t = 0; t < T; t++) {
        for (unsigned int s = 0; s < S; s++) {
          log_py(s, t) = 0;
          for (unsigned int c = 0; c < C; c++) {
            log_py(s, t) += log_B(c)(s, obs(c, t, i), t);
          }
        }
      }
      log_alpha.slice(i).rows(d * S, (d + 1) * S - 1) = log_omega(d) + univariate_forward_nhmm(log_Pi, log_A, log_py);
    }
  }
  return log_alpha;
}
