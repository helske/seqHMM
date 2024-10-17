//Viterbi algorithm for NHMM and MHMM, single_sequence
#include "get_parameters.h"
#include "eta_to_gamma.h"
#include "viterbi_nhmm.h"

double univariate_viterbi_nhmm(
    const arma::vec& log_init, 
    const arma::cube& log_transition, 
    const arma::mat& log_py,
    arma::subview_col<unsigned int> q) {
  
  unsigned int S = log_py.n_rows;
  unsigned int T = log_py.n_cols;
  
  arma::mat delta(S, T);
  arma::umat phi(S, T);
  delta.col(0) = log_init + log_py.col(0);
  phi.col(0).zeros();
  for (unsigned int t = 1; t < T; t++) {
    for (unsigned int j = 0; j < S; j++) {
      phi(j, t) = (delta.col(t - 1) + log_transition.slice(t).col(j)).index_max();
      delta(j, t) = delta(phi(j, t), t - 1) + log_transition(phi(j, t), j, t) + log_py(j, t);
    }
  }
  q(T - 1) = delta.col(T - 1).index_max();
  for (int t = (T - 2); t >= 0; t--) {
    q(t) = phi(q(t + 1), t + 1);
  }
  return delta.col(T - 1).max();
}
// [[Rcpp::export]]
Rcpp::List viterbi_nhmm_singlechannel(
    const arma::mat& eta_pi, const arma::mat& X_i,
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::cube& eta_B, const arma::cube& X_o,
    const arma::umat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int M = eta_B.n_rows + 1;
  arma::umat q(T, N, arma::fill::zeros);
  arma::vec logp(N);
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
    log_B = get_log_B(gamma_B, X_o.slice(i), true, true);
    for (unsigned int t = 0; t < T; t++) {
      log_py.col(t) = log_B.slice(t).col(obs(t, i));
    }
    logp(i) = univariate_viterbi_nhmm(log_Pi, log_A, log_py, q.col(i));
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}

// [[Rcpp::export]]
Rcpp::List viterbi_nhmm_multichannel(
    const arma::mat& eta_pi, const arma::mat& X_i,
    const arma::cube& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::ucube& obs, const arma::uvec M) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A.n_slices;
  unsigned int C = obs.n_rows;
  arma::umat q(T, N, arma::fill::zeros);
  arma::vec logp(N);
  arma::mat log_py(S, T);
  arma::vec log_Pi(S);
  arma::cube log_A(S, S, T);
  arma::field<arma::cube> log_B(C);
  arma::mat gamma_pi = eta_to_gamma(eta_pi);
  arma::cube gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B= eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_py.zeros();
    log_Pi = get_log_pi(gamma_pi, X_i.col(i));
    log_A = get_log_A(gamma_A, X_s.slice(i));
    log_B = get_log_B(gamma_B, X_o.slice(i), M, true, true);
    for (unsigned int t = 0; t < T; t++) {
      for (unsigned int c = 0; c < C; c++) {
        log_py.col(t) += log_B(c).slice(t).col(obs(c, t, i));
      }
    }
    logp(i) = univariate_viterbi_nhmm(log_Pi, log_A, log_py, q.col(i));
  }
  
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}

// [[Rcpp::export]]
Rcpp::List viterbi_mnhmm_singlechannel(
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_i,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d,
    const arma::umat& obs) {
  
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int SD = S * D;
  unsigned int M = eta_B(0).n_rows + 1;
  arma::umat q(T, N, arma::fill::zeros);
  arma::vec logp(N);
  arma::mat log_py(SD, T);
  arma::vec log_Pi(SD);
  arma::cube log_A(SD, SD, T, arma::fill::value(-arma::datum::inf));
  arma::cube log_B(SD, M + 1, T);
  arma::vec log_omega(D);
  arma::mat gamma_omega = eta_to_gamma(eta_omega);
  arma::field<arma::mat> gamma_pi = eta_to_gamma(eta_pi);
  arma::field<arma::cube> gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_omega = get_log_omega(gamma_omega, X_d.col(i));
    for (unsigned int d = 0; d < D; d++) {
      log_Pi.rows(d * S, (d + 1) * S - 1) = log_omega(d) + 
        get_log_pi(gamma_pi(d), X_i.col(i));
      log_A.tube(d * S, d * S, (d + 1) * S - 1, (d + 1) * S - 1) = 
        get_log_A(gamma_A(d), X_s.slice(i));
      log_B.rows(d * S, (d + 1) * S - 1) = get_log_B(
        gamma_B(d), X_o.slice(i), true, true
      );
    }
    for (unsigned int t = 0; t < T; t++) {
      log_py.col(t) = log_B.slice(t).col(obs(t, i));
    }
    logp(i) = univariate_viterbi_nhmm(log_Pi, log_A, log_py, q.col(i));
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
// [[Rcpp::export]]
Rcpp::List viterbi_mnhmm_multichannel(
    const arma::field<arma::mat>& eta_pi, const arma::mat& X_i,
    const arma::field<arma::cube>& eta_A, const arma::cube& X_s,
    const arma::field<arma::cube>& eta_B, const arma::cube& X_o,
    const arma::mat& eta_omega, const arma::mat& X_d,
    const arma::ucube& obs, const arma::uvec M) {
  unsigned int N = X_s.n_slices;
  unsigned int T = X_s.n_cols;
  unsigned int S = eta_A(0).n_slices;
  unsigned int D = eta_omega.n_rows + 1;
  unsigned int SD = S * D;
  unsigned int C = obs.n_rows;
  arma::umat q(T, N, arma::fill::zeros);
  arma::vec logp(N);
  arma::mat log_py(SD, T);
  arma::vec log_Pi(SD);
  arma::cube log_A(SD, SD, T, arma::fill::value(-arma::datum::inf));
  arma::field<arma::cube> log_B(C);
  arma::vec log_omega(D);
  arma::mat gamma_omega = eta_to_gamma(eta_omega);
  arma::field<arma::mat> gamma_pi = eta_to_gamma(eta_pi);
  arma::field<arma::cube> gamma_A = eta_to_gamma(eta_A);
  arma::field<arma::cube> gamma_B = eta_to_gamma(eta_B);
  for (unsigned int i = 0; i < N; i++) {
    log_omega = get_log_omega(gamma_omega, X_d.col(i));
    for (unsigned int d = 0; d < D; d++) {
      log_Pi.rows(d * S, (d + 1) * S - 1) = log_omega(d) + 
        get_log_pi(gamma_pi(d), X_i.col(i));
      log_A.tube(d * S, d * S, (d + 1) * S - 1, (d + 1) * S - 1) = 
        get_log_A(gamma_A(d), X_s.slice(i));
      log_B = get_log_B(
        gamma_B.rows(d * C, (d + 1) * C - 1), X_o.slice(i), M, true, true
      );
      for (unsigned int t = 0; t < T; t++) {
        for (unsigned int s = 0; s < S; s++) {
          log_py(d * S + s, t) = 0;
          for (unsigned int c = 0; c < C; c++) {
            log_py(d * S + s, t) += log_B(c)(s, obs(c, t, i), t);
          }
        }
      }
    }
    logp(i) = univariate_viterbi_nhmm(log_Pi, log_A, log_py, q.col(i));
  }
  return Rcpp::List::create(
    Rcpp::Named("q") = Rcpp::wrap(q), 
    Rcpp::Named("logp") = Rcpp::wrap(logp)
  );
}
